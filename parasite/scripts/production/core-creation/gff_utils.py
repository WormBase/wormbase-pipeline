import pandas as pd 
import re
from Bio import SeqIO
from argparse import ArgumentParser

def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_gff",
        help = "Unprocessed .gff file which requires parsing.")
    parser.add_argument("output_file",
    	help = "Output .gff file to write processed GFF to.")
    args = parser.parse_args()
    return args 


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


# Reads in a .gff file as a Pandas dataframe. Rows and values are not
# modified but the attributes column is parsed into separate columns,
# one for each identified attribute (e.g. ID, Name, Parent).
# Attributes not present in the original .gff file are NOT added.
def parse_gff(input_gff):
	col_names = ["scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]

	gff_df = pd.read_csv(input_gff, sep = '\t', comment = '#', names = col_names, index_col = False)

	gff_df = _parse_attributes_field(gff_df)

	return gff_df

# Copies features ID attributes into Name attributes
def infer_name_attribute_from_id(gff_df):
	gff_df["Name"] = gff_df["ID"]
	return gff_df


# Checks if a Name attribute is present. Returns False if not.
def has_name_attribute(gff_df):
	return ("Name" in gff_df)


# Adds specified string prefix to the ID attribute of all features.
# Might be better to combine with add_prefix_to_name, or call with
# wrapper function that calls both.
def add_prefix_to_id(gff_df, prefix):
	gff_df["ID"] = prefix + gff_df["ID"]
	return gff_df


# Adds specified string prefix to the Name attribute of all features.
# Might be better to combine with add_prefix_to_id, or call with
# wrapper function that calls both.
def add_prefix_to_name(gff_df, prefix):
	gff_df["Name"] = prefix + gff_df["Name"]
	return gff_df


def get_all_gene_features(gff_df):
	gene_mask = gff_df.loc[:, "type"] == "gene"
	gene_df = gff_df.loc[gene_mask]
	return gene_df


def get_all_transcript_features(gff_df):
	transcript_types = ["mRNA", "tRNA", "rRNA", "transcript"]
	transcript_mask = gff_df.loc[:, "type"].isin(transcript_types) 
	transcript_df = gff_df.loc[transcript_mask]
	return transcript_df


def get_all_exon_features(gff_df):
	exon_mask = gff_df.loc[:, "type"] == "exon"
	exon_df = gff_df.loc[exon_mask]
	return exon_df


def get_all_cds_features(gff_df):
	cds_mask = gff_df.loc[:, "type"] == "CDS"
	cds_df = gff_df.loc[cds_mask]
	return cds_df


# Returns numpy array containing all unique type values in .gff
def get_all_unique_types(gff_df):
	return pd.unique(gff_df["type"])

# Returns numpy array containing all unique source values in .gff
def get_all_unique_sources(gff_df):
	return pd.unique(gff_df["source"])

def print_all_unique_types(gff_df):
	print(".gff file contains the following type values: ")
	for unique_type in get_all_unique_types(gff_df):
		if unique_type not in ["gene", "mRNA", "tRNA", "rRNA", "exon", "CDS"]:
			suffix = "\t<---- Non-standard type"
		else:
			suffix = ""
		print("\t" + str(unique_type) + suffix)


def print_all_unique_sources(gff_df):
	print(".gff file contains the following source values: ")
	for unique_source in get_all_unique_sources(gff_df):
		print("\t" + str(unique_source))


# Recreates attribute field according to ParaSite specification.
# TODO: Pretty crude implementation. There might be a better way to 
# concatenate string columns with pandas
def finalise_attributes_column(gff_df):
	gff_df["attributes"] = "ID=" + gff_df["ID"] + ';' + "Name=" + gff_df["Name"] + ';' + "Parent=" + gff_df["Parent"]
	return gff_df

# Drops any columns created during .gff processing but not required
# in the final .gff file. Only call right before writing .gff to file.
# TODO: Could refer to a global list of final columns for maintainability
def drop_superfluous_columns(gff_df):
	final_columns = ["scaffold", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
	gff_df = gff_df.drop(columns=[column for column in gff_df if column not in final_columns])
	return gff_df

# Updates attributes field and drops unecessary columns, then writes
# processed .gff to a specified output file.
def output_gff(gff_df, output_file):
	gff_df = finalise_attributes_column(gff_df)
	gff_df = drop_superfluous_columns(gff_df)

	gff_df.to_csv(output_file, sep = '\t')

args = get_args()

gff_df = parse_gff(args.input_gff)

# Infer Name attribute from ID if absent
if not has_name_attribute(gff_df):
	gff_df = infer_name_attribute_from_id(gff_df)


print(gff_df.head(n=20))

print_all_unique_types(gff_df)
print_all_unique_sources(gff_df)

# gene_df = get_all_gene_features(gff_df)
# transcript_df = get_all_transcript_features(gff_df)
# exon_df = get_all_exon_features(gff_df)
# cds_df = get_all_cds_features(gff_df)

# output_gff(gff_df, args.output_file)




