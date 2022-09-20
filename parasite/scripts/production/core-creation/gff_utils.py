import pandas as pd 
import re
from Bio import SeqIO
from argparse import ArgumentParser
from datetime import datetime

# def get_args():
#     parser = ArgumentParser()
#     parser.add_argument("input_gff",
#         help = "Unprocessed .gff file which requires parsing.")
#     parser.add_argument("output_file",
#     	help = "Output .gff file to write processed GFF to.")
#     parser.add_argument("synonyms_file",
#     	help = "Synonyms.tsv file for mapping scaffold names.")
#     args = parser.parse_args()
#     return args 


# Removes provided prefixes from all columns in a .gff3 dataframe.
def remove_prefixes_from_column_values(gff_df, prefixes):
	for prefix in prefixes:
		# Look only at beginning of line, escape any special characters.
		prefix = "^"+re.escape(prefix)
		for col in gff_df:
			# Avoid using re.sub on NaN values as this turns them into str
			nan_mask = gff_df.loc[:,col].isna()
			gff_df.loc[~nan_mask, col] = [re.sub(prefix, "", str(x)) for x in gff_df.loc[~nan_mask, col]]
	return gff_df


# Takes a GFF dataframe and adds a column for each attribute specified
# in the attributes column (column 9). Does not remove the original
# attributes column. NaN inserted where features lack a value for an
# attribute.
def _parse_attributes_field(gff_df):

	# Iterate over rows as tuples (~3x faster than iterrows)
	# Column 0 is index, column 9 is attributes 
	for row in gff_df.itertuples():

		# .gff attributes are delimited by ';'
		attribute_fields = row[9].split(';')

		for attribute_field in attribute_fields:

			# Skip empty rows output by .split
			if attribute_field != "":

				# Format is "attribute"="value"
				attribute_name = attribute_field.split('=')[0]
				attribute_value = attribute_field.split('=')[1]
				gff_df.at[row[0], attribute_name] = attribute_value

	return gff_df


# Makes IDs unique by adding | N as a suffix, where N is the cumulative
# occurence of a non-unique ID. The first occurence is left unchanged.
def _make_ids_unique(gff_df):
	# Suffixes cumulative count to a value.
	def _convert_count_to_unique_suffix(count):
		# 0 indicates first instance of non-unique ID -> empty string.
		if count == 0:
			return ""
		else:
			return "|" + str(count)

	# Modify all non-unique IDs with their cumulative count.
	gff_df["ID"] += gff_df.groupby("ID").cumcount().apply(_convert_count_to_unique_suffix) 
	return gff_df


# Reads in a .gff file as a Pandas dataframe. Rows and values are not
# modified but the attributes column is parsed into separate columns,
# one for each identified attribute (e.g. ID, Name, Parent).
# Attributes not present in the original .gff file are NOT added.
def parse_gff(input_gff):
	col_names = ["scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

	gff_df = pd.read_csv(input_gff, sep = '\t', comment = '#', names = col_names, index_col = False)

	gff_df = _parse_attributes_field(gff_df)

	gff_df = _make_ids_unique(gff_df)

	return gff_df


# Replaces INSDCID scaffold names with community scaffold names using 
# synonyms mapping provided in a synonyms.tsv file.
# TODO: Spits out NaNs as scaffold names if the .gff is already using
# community scaffold names. Implement a check for this? Or rely on 
# correct calls on the user's part? 
def rename_scaffolds(gff_df, synonyms_file):
	# TODO: Infer synonyms headers from file rather than hardcoding as below.
	colnames = ["toplevel", "CommunityID", "INSDCID", "INSDC"]
	scaffold_df = pd.read_csv (synonyms_file, sep = '\t', names = colnames, usecols = ["CommunityID", "INSDCID"])
	scaffold_df.set_index("INSDCID", inplace = True)
	gff_df.loc[:, "scaffold"] = gff_df["scaffold"].map(scaffold_df["CommunityID"])
	return gff_df



# Copies features ID attributes into Name attributes
def infer_and_overwrite_name_attribute_from_id(gff_df):
	gff_df["Name"] = gff_df["ID"]
	return gff_df


# As infer_and_overwrite_name_attribute_from_id, but should leave non-empty
# Name fields unchanged (keeps original Name found in .gff3.)
def infer_name_attribute_from_id(gff_df):
	no_name_mask = gff_df.loc[:, "Name"].isna()

	gff_df.loc[no_name_mask, "Name"] = gff_df.loc[no_name_mask, "ID"]
	
	gff_df["Name"] = gff_df["ID"]
	return gff_df

# Creates mRNA and exon features for genes on a given scaffold where only
# gene and CDS features are provided. Returns df with fabricated features
def extrapolate_missing_transcripts_for_scaffold(gff_df, scaffold):
	# Grab all features on scaffold
	scaffold_mask = gff_df.loc[:,"scaffold"] == scaffold 
	scaffold_df = gff_df.loc[scaffold_mask]

	scaffold_genes = get_all_gene_features(scaffold_df)

	scaffold_cds = get_all_cds_features(scaffold_df)

	# Extrapolate transcript features from gene features
	scaffold_transcripts = scaffold_genes.copy(deep=True).reset_index()

	scaffold_transcripts["type"] = "mRNA"
	# Set transcript's parent to gene it is based on.
	scaffold_transcripts["Parent"] = scaffold_transcripts["ID"]

	scaffold_transcripts["ID"] += "_mRNA"
	scaffold_transcripts["Name"] = scaffold_transcripts["ID"]

	# Change CDS parent from gene to corresponding, newly created transcript.
	gene_transcript_dict = dict(zip(scaffold_transcripts["Parent"], scaffold_transcripts["ID"]))

	scaffold_cds["Parent"] = scaffold_cds["Parent"].map(gene_transcript_dict)

	# Extrapolate exon features from CDS features
	scaffold_exons = scaffold_cds.copy(deep=True).reset_index()

	scaffold_exons["type"] = "exon"

	# Name exon names after their parent transcript, appending a number.
	scaffold_exons["ID"] = scaffold_exons["Parent"]
	scaffold_exons["ID"] += "-" + (scaffold_exons.groupby("ID").cumcount()+1).astype(str)
	scaffold_exons["Name"] = scaffold_exons["ID"]
	# Exons don't have a phase, unlike CDS features.
	scaffold_exons["phase"] = "."

	scaffold_df = pd.concat([scaffold_genes, scaffold_transcripts, scaffold_cds, scaffold_exons]).reset_index()

	# Replace the original features for given scaffold in the .gff
	gff_df.drop(gff_df.loc[scaffold_mask].index, inplace = True)

	gff_df = pd.concat([gff_df, scaffold_df])

	return(gff_df)

# Checks if a Name attribute is present. Returns False if not.
def has_name_attribute(gff_df):
	return ("Name" in gff_df)

#TODO: Might be better to combine the below add_prefix functions or write
# a wrapper function that calls both


# Adds specified string prefix to the ID attribute of all features.
def add_prefix_to_id(gff_df, prefix):
	gff_df["ID"] = prefix + gff_df["ID"]
	return gff_df


# Adds specified string prefix to the Name attribute of all features.
def add_prefix_to_name(gff_df, prefix):
	gff_df["Name"] = prefix + gff_df["Name"]
	return gff_df


# Replaces names in source column with provided argument.
def rename_sources(gff_df, new_source):
	gff_df["source"] = new_source
	return gff_df


# Returns dataframe consisting only of gene features.
def get_all_gene_features(gff_df):
	gene_mask = gff_df.loc[:, "type"] == "gene"
	gene_df = gff_df.loc[gene_mask]
	return gene_df


# Returns dataframe consisting only of transcript features.
def get_all_transcript_features(gff_df):
	# TODO: Could refer to a global list of features for maintainability
	transcript_types = ["mRNA", "tRNA", "rRNA", "transcript"]
	transcript_mask = gff_df.loc[:, "type"].isin(transcript_types) 
	transcript_df = gff_df.loc[transcript_mask]
	return transcript_df


# Returns dataframe consisting only of exon features.
def get_all_exon_features(gff_df):
	exon_mask = gff_df.loc[:, "type"] == "exon"
	exon_df = gff_df.loc[exon_mask]
	return exon_df


# Returns dataframe consisting only of CDS features.
def get_all_cds_features(gff_df):
	cds_mask = gff_df.loc[:, "type"] == "CDS"
	cds_df = gff_df.loc[cds_mask]
	return cds_df


# Creates Parent Gene column for all gene features. This contains the ID
# the gene's own gene ID (quick copy). Used for sorting the final .gff.
def _get_parent_gene_for_genes(gene_df):
	# TODO: Still raises a chained assignment warning:
	gene_df.loc[:,"Parent_Gene"] = gene_df["ID"]
	return gene_df 


# Creates a Parent Gene column for all transcript features. This contains
# the ID of the transcripts' parent gene (quick copy from Parent col).
# Used for sorting the final .gff.
def _get_parent_gene_for_transcripts(transcript_df):
	# TODO: Still raises a chained assignment warning
	transcript_df.loc[:,["Parent_Gene"]] = transcript_df["Parent"]
	return transcript_df


# Creates Parent Gene column for all exon or CDS features. This contains
# the ID of the exon's/CDS's parent gene. It has to be inferred from its 
# parent transcript, requiring iterrows instead of vectorised copying. 
# Used for sorting the final .gff.
def _get_parent_gene_for_exons_or_cds(exon_or_cds_df, transcript_df):
	# TODO: May need dramatic speeding up due to reliance on iterrows.
	for index, exon in exon_or_cds_df.iterrows():
		# TODO: Still raises a chained assignment warning:
		parent_transcript_id = exon["Parent"]
		parent_gene_id = transcript_df.at[parent_transcript_id, "Parent_Gene"]
		exon_or_cds_df.at[index, "Parent_Gene"] = parent_gene_id
	return exon_or_cds_df


# Returns numpy array containing all unique type values in .gff
def _get_all_unique_types(gff_df):
	return pd.unique(gff_df["type"])


# Returns numpy array containing all unique source values in .gff
def _get_all_unique_sources(gff_df):
	return pd.unique(gff_df["source"])


# Prints all unique type values across all features to stdout.
# Points out any non-standard types.
def print_all_unique_feature_types(gff_df):
	print(".gff file contains the following feature types: ")
	for unique_type in _get_all_unique_types(gff_df):
		# Could be refactored so that standard types are not defined in function.
		if unique_type not in ["gene", "mRNA", "tRNA", "rRNA", "exon", "CDS"]:
			suffix = "\t<---- Non-standard type"
		else:
			suffix = ""
		print("\t" + str(unique_type) + suffix)


def print_all_unique_sources(gff_df):
	print(".gff file contains the following source values: ")
	for unique_source in _get_all_unique_sources(gff_df):
		print("\t" + str(unique_source))


# Recreates attribute field according to ParaSite specification.
def finalise_attributes_column(gff_df):
	gff_df["attributes"] = "ID=" + gff_df["ID"].fillna('') + ';' + "Name=" + gff_df["Name"].fillna('') + ';' + "Parent=" + gff_df["Parent"].fillna('')

	# Omit Parent attribute in features that lack it (e.g. genes)
	gff_df["attributes"] = gff_df["attributes"].mask(gff_df["Parent"].isna(), 
		"ID=" + gff_df["ID"].fillna('') + ';' + "Name=" + gff_df["Name"].fillna(''))
	return gff_df


# Drops any columns created during .gff processing but not required
# in the final .gff file. Only call right before writing .gff to file.
def _drop_superfluous_columns(gff_df):
	# TODO: Could refer to a global list of final columns for maintainability
	final_columns = ["scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
	gff_df = gff_df.drop(columns=[column for column in gff_df if column not in final_columns])
	return gff_df


# Drops any columns created during .gff processing that are not required
# for final processing steps (e.g. sorting features).
def _drop_useless_columns(gff_df):
	useful_columns = ["scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes", "ID", "Parent", "Name"]
	gff_df = gff_df.drop(columns=[column for column in gff_df if column not in useful_columns])
	return gff_df


# Outputs a finalised .gff dataframe that contains only gene, transcript,
# exon and CDS features. Orders the dataframe so that features exist in a
# hierarchy: gene, gene's transcripts, transcript's exons/CDSs
def reorder_gff_features(gff_df):
	gff_df.set_index("ID", drop = False, inplace = True)

	gene_df = get_all_gene_features(gff_df)
	transcript_df = get_all_transcript_features(gff_df)
	exon_df = get_all_exon_features(gff_df)

	cds_df = get_all_cds_features(gff_df)
	


	# Create a "Parent Gene" column for each feature to allow sorting.
	gene_df = _get_parent_gene_for_genes(gene_df)
	transcript_df  = _get_parent_gene_for_transcripts(transcript_df)

	exon_df = _get_parent_gene_for_exons_or_cds(exon_df, transcript_df)
	exon_df.to_csv("exond_df.tsv", sep = '\t')
	cds_df = _get_parent_gene_for_exons_or_cds(cds_df, transcript_df)
	cds_df.to_csv("cds_df.tsv", sep = '\t')


	# Stick all features together and sort them by parent gene.
	reordered_gff = pd.concat([gene_df, transcript_df, exon_df, cds_df])

	reordered_gff["type"] = pd.Categorical(reordered_gff["type"], ["gene", "mRNA", "tRNA", "rRNA", "transcript", "exon", "CDS"])

	# TODO: Sort is purely alphabetic (e.g. g10 comes after g1, should be g2)
	reordered_gff = reordered_gff.sort_values(by=["Parent_Gene", "type"])
	return reordered_gff


# Updates attributes field and drops unecessary columns, then writes
# processed .gff to a specified output file.
def write_output_gff(gff_df, output_file):
	gff_df = finalise_attributes_column(gff_df)
	gff_df = _drop_useless_columns(gff_df)
	gff_df = reorder_gff_features(gff_df)
	gff_df = _drop_superfluous_columns(gff_df)

	with open(output_file, 'w') as output_gff:
		output_gff.write("##gff-version 3\n")

	gff_df.to_csv(output_file, mode = 'a', header = False, sep = '\t', index = False)


