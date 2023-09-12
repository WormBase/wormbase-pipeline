import pandas as pd
import numpy as np
import re
import sys
from Bio import SeqIO
from argparse import ArgumentParser
from datetime import datetime
from ProductionUtils import *
from reformat_gff_dependencies import *

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
    gff_df = gff_df.replace(prefixes, "", regex=True)
    return gff_df

# Removes provided prefixes from the 'Name' column in a .gff3 dataframe.
def remove_prefixes_from_name_column(gff_df, prefixes):
    gff_df['Name'] = gff_df['Name'].replace(prefixes, "", regex=True)
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

def parse_synonyms(synonyms_file):
    colnames = ["toplevel", "CommunityID", "INSDCID", "INSDC"]
    scaffold_df = pd.read_csv (synonyms_file, sep = '\t', names = colnames, usecols = ["CommunityID", "INSDCID"])
    # scaffold_df.set_index("INSDCID", inplace = True)
    return scaffold_df

def scaffolds_needs_renaming(gff_df, fasta):
    # If gff_df and fasta have the same scaffold names then no renaming is needed
    # Extract the unique scaffold values from the gff_df DataFrame
    unique_scaffolds_in_df = gff_df["scaffold"].unique()
    # Get the scaffold names from the fasta object
    scaffold_names_in_fasta = fasta.scaffold_names()
    # Check if all values from gff_df["scaffold"] exist in scaffold_names_in_fasta
    if all(scaffold in scaffold_names_in_fasta for scaffold in unique_scaffolds_in_df):
        return False
    elif all(scaffold not in scaffold_names_in_fasta for scaffold in unique_scaffolds_in_df):
        return True
    else:
        return "Some"

    

def rename_scaffolds_to(gff_df, synonyms_df):
    # if gff_df scaffolds match scaffold_df["CommunityID"] then rename them to scaffold_df["INSDC"] if
    # they match scaffold_df["INSDC"] then rename them to scaffold_df["CommunityID"]
    if gff_df["scaffold"].unique().tolist() == synonyms_df["CommunityID"].unique().tolist():
        return "INSDCID"
    elif gff_df["scaffold"].unique().tolist() == synonyms_df["INSDCID"].unique().tolist():
        return("CommunityID")
    else:
        exit_with_error("GFF scaffolds do not match CommunityID or INSDCID columns in synonyms file. \
                        Cannot rename scaffolds.")



# Replaces INSDCID scaffold names with community scaffold names using
# synonyms mapping provided in a synonyms.tsv file.
# TODO: Spits out NaNs as scaffold names if the .gff is already using
# community scaffold names. Implement a check for this? Or rely on 
# correct calls on the user's part? 
def rename_scaffolds(gff_df, synonyms_df):
    rename_to = rename_scaffolds_to(gff_df, synonyms_df)
    rename_from = [x for x in synonyms_df.columns.to_list() if x !=rename_to][0]
    print_info(f"Renaming GFF scaffolds from {rename_from} to {rename_to}")

    # Merge the DataFrames to replace "scaffold" values
    merged_df = gff_df.merge(synonyms_df, left_on="scaffold", right_on=rename_from, how="left")
    merged_df["scaffold"] = merged_df[rename_to].fillna(merged_df["scaffold"])

    # Drop the unnecessary columns
    merged_df.drop(columns=[rename_from, rename_to], inplace=True)

    if merged_df["scaffold"].isna().any():
        problematic_scaffolds = merged_df.loc[merged_df["scaffold"].isna(), "scaffold"].unique()
        # Exit with error stating the problematic scaffolds
        exit_with_error("Some scaffold names could not rename: " + str(problematic_scaffolds))
    if len(gff_df)!=len(merged_df):
        # Print mismatching rows between gff_df and merged_df
        exit_with_error(f"Mismatching rows between gff_df and merged_df when renaming scaffolds from {rename_from} to {rename_to}")
    return merged_df

# Create an empty Name column if doesnt exist
def create_name_column(gff_df):
    if "Name" not in gff_df.columns:
        gff_df["Name"] = None
    return gff_df


# Copies features ID attributes into Name attributes
def infer_and_overwrite_name_attribute_from_id(gff_df):
    gff_df = create_name_column(gff_df)
    gff_df["Name"] = gff_df["ID"]
    return gff_df


# As infer_and_overwrite_name_attribute_from_id, but should leave non-empty
# Name fields unchanged (keeps original Name found in .gff3.)
def infer_name_attribute_from_id(gff_df):

    gff_df = create_name_column(gff_df)
        
    no_name_mask = gff_df.loc[:, "Name"].isna()

    gff_df.loc[no_name_mask, "Name"] = gff_df.loc[no_name_mask, "ID"]

    gff_df["Name"] = gff_df["ID"]
    return gff_df


# Updates CDS features on a scaffold to point to a transcript as parent
# Assumes the CDSs originally pointed to a gene as a parent and that 
# transcripts were extrapolated from said genes.
def _update_scaffold_cds_parent_attribute(gff_df, scaffold):
    scaffold_mask = gff_df.loc[:,"scaffold"] == scaffold
    scaffold_df = gff_df.loc[scaffold_mask]

    scaffold_transcripts = get_all_transcript_features(scaffold_df)
    scaffold_cds = get_all_cds_features(scaffold_df)

    gene_transcript_dict = dict(zip(scaffold_transcripts["Parent"], scaffold_transcripts["ID"]))

    scaffold_cds["Parent"] = scaffold_cds["Parent"].map(gene_transcript_dict)

    scaffold_cds_mask = (gff_df.loc[:,"scaffold"] == scaffold) & (gff_df.loc[:,"type"] == "CDS")

    # Replace the original features for given scaffold in the .gff
    gff_df.drop(gff_df.loc[scaffold_cds_mask].index, inplace = True)

    gff_df = pd.concat([gff_df, scaffold_cds])

    return gff_df

# Creates transcript features for given scaffold based on extant genes.
# Updates extant CDSs to point to transcripts as parents, returns df.
def extrapolate_scaffold_transcripts_from_genes(gff_df, scaffold, suffix = "_mRNA"):
    scaffold_mask = gff_df.loc[:,"scaffold"] == scaffold
    scaffold_df = gff_df.loc[scaffold_mask]

    scaffold_genes = get_all_gene_features(scaffold_df)

    # Extrapolate transcript features from gene features
    scaffold_transcripts = scaffold_genes.copy(deep=True).reset_index(drop=True)

    scaffold_transcripts["type"] = "mRNA"
    # Set transcript's parent to gene it is based on.
    scaffold_transcripts["Parent"] = scaffold_transcripts["ID"]

    # Modify transcript ID to clarify it is mRNA.
    scaffold_transcripts["ID"] += suffix
    scaffold_transcripts["Name"] = scaffold_transcripts["ID"]

    # Add extrapolated transcripts to original gff dataframe
    gff_df = pd.concat([gff_df, scaffold_transcripts]).reset_index()

    # Update scaffold CDS to point to extrapolated transcripts as parents
    gff_df = _update_scaffold_cds_parent_attribute(gff_df, scaffold)

    return gff_df


# Creates gene features for given scaffold based on extant transcripts.
# Updates extant transcripts to point to genes as parents, returns df.
def extrapolate_scaffold_genes_from_transcripts(gff_df, scaffold, suffix = "_gene"):
    scaffold_mask = gff_df.loc[:,"scaffold"] == scaffold
    scaffold_df = gff_df.loc[scaffold_mask]

    scaffold_transcripts = get_all_transcript_features(scaffold_df)

    # Extrapolate gene features from transcript fetures.
    scaffold_genes = scaffold_transcripts.copy(deep = True).reset_index(drop=True)

    scaffold_genes["type"] = "gene"

    # Modify gene ID to clarify it is a gene and make it distinct.
    scaffold_genes["ID"] += suffix

    # Update transcripts' parent attributes to point to extrapolated genes.
    scaffold_transcripts["Parent"] = scaffold_transcripts["ID"] + suffix

    scaffold_genes["Name"] = scaffold_genes["ID"]

    scaffold_transcript_mask = (gff_df.loc[:,"scaffold"] == scaffold) & (gff_df.loc[:,"type"].isin(transcript_types))

    # Replace the original transcripts for scaffold in the .gff
    gff_df.drop(gff_df.loc[scaffold_transcript_mask].index, inplace = True)

    # Add extrapolated genes to .gff
    gff_df = pd.concat([gff_df, scaffold_genes]).reset_index()

    return gff_df


# Creates exon features for given scaffold based on extant CDSs.
# Exons get the same parent as the CDS. Returns df.
def extrapolate_scaffold_exons_from_cds(gff_df, scaffold):
    # Grab all features on scaffold
    scaffold_mask = gff_df.loc[:,"scaffold"] == scaffold
    scaffold_df = gff_df.loc[scaffold_mask]

    scaffold_cds = get_all_cds_features(scaffold_df)

    # Extrapolate exon features from CDS features
    scaffold_exons = scaffold_cds.copy(deep=True).reset_index(drop=True)

    scaffold_exons["type"] = "exon"

    # Name exon names after their parent transcript, appending a number.
    scaffold_exons["ID"] = scaffold_exons["Parent"]
    scaffold_exons["ID"] += "-" + (scaffold_exons.groupby("ID").cumcount()+1).astype(str)
    scaffold_exons["Name"] = scaffold_exons["ID"]
    # Exons don't have a phase, unlike CDS features.
    scaffold_exons["phase"] = "."

    gff_df = pd.concat([gff_df, scaffold_exons]).reset_index()

    return gff_df


# Creates CDS features for given scaffold based on extant exons.
# CDSs get the same parent as the exon. Returns df.
def extrapolate_scaffold_cds_from_exons(gff_df, scaffold):
    # Grab all features on scaffold
    scaffold_mask = gff_df.loc[:,"scaffold"] == scaffold
    scaffold_df = gff_df.loc[scaffold_mask]

    scaffold_exons = get_all_exon_features(scaffold_df)

    # Extrapolate CDS features from exon features
    scaffold_cds = scaffold_exons.copy(deep=True).reset_index(drop=True)

    scaffold_cds["type"] = "CDS"

    scaffold_cds["ID"] = scaffold_exons["Parent"]
    scaffold_cds["ID"] += "-" + (scaffold_cds.groupby("ID").cumcount()+1).astype(str)
    scaffold_cds["Name"] = scaffold_cds["ID"]

    gff_df = pd.concat([gff_df, scaffold_cds]).reset_index()

    return gff_df


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

def make_coding_transcripts_mRNA(gff_df):
    gff_df.loc[gff_df['type'].isin(coding_transcript_types), ['type']] = 'mRNA'
    return(gff_df)


# Replaces names in source column with provided argument.
def rename_sources(gff_df, new_source):
    gff_df["source"] = new_source
    return gff_df


# Returns dataframe consisting only of gene features.
def get_all_gene_features(gff_df):
    gene_mask = gff_df.loc[:, "type"].isin(gene_types)
    gene_df = gff_df.loc[gene_mask]
    return gene_df


# Returns dataframe consisting only of transcript features.
def get_all_transcript_features(gff_df):
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

# Returns a dataframe containing only exon and CDS features.
def get_all_exon_and_cds_features(gff_df, transcripts=None):
    if transcripts is None:
        exon_cds_df = gff_df[(gff_df["type"] == "CDS") | (gff_df["type"] == "exon")]
    else:
        exon_cds_df = gff_df[((gff_df["type"] == "CDS") | (gff_df["type"] == "exon")) & (gff_df["Parent"].isin(transcripts))]
    return exon_cds_df

# This function takes a dataframe grouped by the "Parent" column and
# returns a pandas Series with the exon/cds metrics below.
# This function is being used by the transcripts_have_exons_and_valid_cds datacheck.
def _count_cds_exons_per_transcript(gff_df):
    exons_cds_stats = {
        'no_of_exons': gff_df[gff_df['type']=='exon']['Parent'].count(),
        'no_of_cds': gff_df[gff_df['type']=='CDS']['Parent'].count(),
        'exons_max_coord': max([gff_df[gff_df['type']=='exon']["start"].max(),gff_df[gff_df['type']=='exon']["end"].max()]),
        'exons_min_coord': min([gff_df[gff_df['type']=='exon']["start"].min(),gff_df[gff_df['type']=='exon']["end"].min()]),
        'cds_max_coord': max([gff_df[gff_df['type']=='CDS']["start"].max(),gff_df[gff_df['type']=='CDS']["end"].max()]),
        'cds_min_coord': min([gff_df[gff_df['type']=='CDS']["start"].min(),gff_df[gff_df['type']=='CDS']["end"].min()]),
    }
    return pd.Series(exons_cds_stats)

def get_parent_gene_for_all(gff_df):
    gff_df.set_index("ID", drop = False, inplace = True)

    gene_df = get_all_gene_features(gff_df)
    transcript_df = get_all_transcript_features(gff_df)
    exon_df = get_all_exon_features(gff_df)
    cds_df = get_all_cds_features(gff_df)

    # Create a "Parent Gene" column for each feature to allow sorting.
    gene_df = _get_parent_gene_for_genes(gene_df)

    transcript_df  = _get_parent_gene_for_transcripts(transcript_df)

    exon_df = _get_parent_gene_for_exons_or_cds(exon_df, transcript_df)
    # exon_df.to_csv("exond_df.tsv", sep = '\t')

    cds_df = _get_parent_gene_for_exons_or_cds(cds_df, transcript_df)
    # cds_df.to_csv("cds_df.tsv", sep = '\t')

    # Stick all features together and sort them by parent gene.
    reordered_gff = pd.concat([gene_df, transcript_df, exon_df, cds_df])

    return(reordered_gff)


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
    exon_or_cds_df.reset_index(drop = True, inplace = True)
    transcript_df.reset_index(drop = True, inplace = True)
    exon_or_cds_tr_merged = exon_or_cds_df.merge(transcript_df[["ID","Parent"]], left_on="Parent", right_on="ID", suffixes=(None, "_transcript"))
    exon_or_cds_tr_merged.rename(columns={'Parent_transcript': 'Parent_Gene'}, inplace=True)
    exon_or_cds_tr_merged.drop(columns=['ID_transcript'])
    return exon_or_cds_tr_merged

# Gets a condition (or a semicolon-separated list of them) in the format of <attribute_field>=<something> like
# gene_status=other and outputs the input gff without the features having the condition specified in their attributes
# as well as their children features. It also outputs a separate gff with these filtered features.
def extract_genes_and_features_with_gene_attribute_value(gff_df, attribute_condition_input):

    attribute_conditions = attribute_condition_input.split(';')
    attribute_conditions_dict = {attribute_field.split('=')[0].strip():attribute_field.split('=')[1].strip() for attribute_field in attribute_conditions}

    query = ' & '.join(['{}=="{}"'.format(k, v) for k, v in attribute_conditions_dict.items()]) + ' & type == @gene_types'

    extracted_genes = list(set(gff_df.query(query)['ID']))
    mask_var = (gff_df['ID'].isin(extracted_genes)) | (gff_df['Parent_Gene'].isin(extracted_genes))
    extracted_gff_df = gff_df[mask_var]
    original_gff_df = gff_df[~mask_var]

    return original_gff_df, extracted_gff_df

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
        if unique_type not in ([gene_types + ["exon", "CDS"] + transcript_types]):
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
    # Make sure that ID/Name/Parent attributes will be omitted if NA
    gff_df["outID"] = np.where(gff_df["ID"].isna(), '', 'ID=' + gff_df["ID"] + ';')
    gff_df["outName"] = np.where(gff_df["Name"].isna(), '', 'Name=' + gff_df["Name"] + ';')
    gff_df["outParent"] = np.where(gff_df["Parent"].isna(), '', 'Parent=' + gff_df["Parent"] + ';')

    # Synthesise the final attributes column
    gff_df["attributes"] = gff_df['outID'] + gff_df['outName'] + gff_df['outParent']

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

# Checks if pseudogenes have CDSs and returns a list with the pseudogenes that have CDSs
def pseudogenes_with_CDS(gff_df):
    pseudogenes = list(set(gff_df[gff_df["type"].isin(pseudogene_types)]["ID"].to_list()))

    pseudogenic_transcripts = list(set(gff_df[gff_df["type"].isin(transcript_types) & gff_df["Parent"].isin(pseudogenes)]["ID"].to_list()))

    pseudogenic_CDS_mask = gff_df["type"].isin(cds_types) & gff_df["Parent"].isin(pseudogenic_transcripts)

    # check if there are any pseudogenic_CDS_mask in the gff_df
    if pseudogenic_CDS_mask.any():
        return True
    else:
        return False

# Removes CDSs of pseudogenes from the gff_df.
def remove_pseudogenic_CDS(gff_df):
    pseudogenes = list(set(gff_df[gff_df["type"].isin(pseudogene_types)]["ID"].to_list()))

    pseudogenic_transcripts = list(set(gff_df[gff_df["type"].isin(transcript_types) & gff_df["Parent"].isin(pseudogenes)]["ID"].to_list()))

    pseudogenic_CDS_mask = gff_df["type"].isin(cds_types) & gff_df["Parent"].isin(pseudogenic_transcripts)

    # Remove the gff_df rows from pseudogenic_CDS_mask
    gff_df = gff_df[~pseudogenic_CDS_mask]

    return(gff_df) 

# Switches the type of transcripts of pseudogenes to "pseudogenic_transcript"
def make_pseudogenic_transcripts(gff_df):
    pseudogenes = list(set(gff_df[gff_df["type"].isin(pseudogene_types)]["ID"].to_list()))

    pseudogenic_transcript_mask = gff_df["type"].isin(transcript_types) & gff_df["Parent"].isin(pseudogenes)

    # Rename all the gff_df["type"] values to "pseudogenic_transcript" in for pseudogenic_transcript_mask 
    gff_df.loc[pseudogenic_transcript_mask, "type"] = "pseudogenic_transcript"

    return(gff_df)

# Switches the type of pseudogenes to "genes":
def pseudogenes_to_genes(gff_df):
    pseudogenes = list(set(gff_df[gff_df["type"].isin(pseudogene_types)]["ID"].to_list()))

    pseudogenes_mask = gff_df["type"].isin(gene_types) & gff_df["ID"].isin(pseudogenes)

    # Rename all the gff_df["type"] values to "gene" in for pseudogenes 
    gff_df.loc[pseudogenes_mask, "type"] = "gene"

    return(gff_df)


# Outputs a finalised .gff dataframe that contains only gene, transcrsipt,
# exon and CDS features. Orders the dataframe so that features exist in a
# hierarchy: gene, gene's transcripts, transcript's exons/CDSs
def reorder_gff_features(gff_df):
    gff_df.set_index("ID", drop=False, inplace=True)

    reordered_gff = get_parent_gene_for_all(gff_df)

    reordered_gff["type"] = pd.Categorical(reordered_gff["type"], list((gene_types + transcript_types + ["exon", "CDS"])))

    # TODO: Sort is purely alphabetic (e.g. g10 comes after g1, should be g2)
    reordered_gff = reordered_gff.sort_values(by=["Parent_Gene", "type"])
    return reordered_gff

# Class to handle genome's FASTA file
class FASTA:
    def __init__(self, fasta_path):
        self.path = fasta_path
    def scaffold_names(self):
        """Returns a list with all the unique scaffold names of the gff dataframe"""
        ids = []
        for record in SeqIO.parse(self.path, "fasta"):
            ids.append(record.id)
        if len(ids)!=len(set(ids)):
            sys.exit("There are non unique scaffold names in FASTA file: " + self.path)
        return (ids)

# Class tha performs datachecks on a gff dataframe
class DATACHECK:
    def __init__(self, gff_df, fasta, args, outprefix=""):
        self.gff_df = gff_df
        self.fasta = fasta
        self.args = args
        self.outprefix = outprefix

    def perform_datachecks(self):
        if self.args.dc_fasta_gff_scaffold or self.args.dc_cds_but_no_exons or self.args.dc_cds_within_exons:
            print_info("Performing DCs: " + ("for " + self.outprefix if self.outprefix!="" else ""))
        self.genes_have_names()
        if self.args.dc_fasta_gff_scaffold:
            self.gff_and_fasta_scaffold_names_match()
        if self.args.dc_cds_but_no_exons or self.args.dc_cds_within_exons:
            self.transcripts_have_exons_and_valid_cds(dc_cds_but_no_exons=self.args.dc_cds_but_no_exons,
                                                           dc_cds_but_no_exons_fix=self.args.dc_cds_but_no_exons_fix,
                                                           dc_cds_within_exons=self.args.dc_cds_within_exons,
                                                           dc_cds_within_exons_fix=self.args.dc_cds_within_exons_fix,
                                                           dc_coding_transcripts_with_cds=self.args.dc_coding_transcripts_with_cds,
                                                           dc_coding_transcripts_with_cds_fix=self.args.dc_coding_transcripts_with_cds_fix,
                                                           prefix=self.outprefix)

    def genes_have_names(self):
        # Filtering the DataFrame
        gene_df = get_all_gene_features(self.gff_df)
        # If gene_df doesn't have a "Name" column, exit with error:
        if not "Name" in gene_df.columns or gene_df["Name"].isna().any():
            exit_with_error("Some/All genes do not have a Name column. Please add a Name column to your GFF file or use the --infer_gene_name/-n option.")

    
    def gff_and_fasta_scaffold_names_match(self):
        """Checks if the GFF's and FASTA file's scaffold names match"""
        gff_scaffolds = list(self.gff_df["scaffold"].unique())
        fasta_scaffolds = self.fasta.scaffold_names()

        gff_to_fasta_discrepancies = [str(x) for x in gff_scaffolds if x not in fasta_scaffolds]
        fasta_to_gff_discrepancies = [str(x) for x in fasta_scaffolds if x not in gff_scaffolds]

        if len(gff_to_fasta_discrepancies) > 0:
            with open(self.outprefix+'in_GFF_not_in_FASTA_scaffolds.txt', mode='wt') as myfile:
                myfile.write('\n'.join(gff_to_fasta_discrepancies))
            exit_with_error("These scaffolds exist in your final gff but not in the final fasta file:\n" + "\n".join(gff_to_fasta_discrepancies))

        if len(fasta_to_gff_discrepancies) > 0:
            print_warning("These scaffolds were found in the FASTA file but not in your final GFF file. Is this ok?\n" + "\n".join(fasta_to_gff_discrepancies))
            with open(self.outprefix+'in_FASTA_not_in_GFF_scaffolds.txt', mode='wt') as myfile:
                myfile.write('\n'.join(fasta_to_gff_discrepancies))

    def transcripts_have_exons_and_valid_cds(self, dc_cds_but_no_exons, dc_cds_but_no_exons_fix, dc_cds_within_exons, dc_cds_within_exons_fix,
                                             dc_coding_transcripts_with_cds, dc_coding_transcripts_with_cds_fix, prefix):
        """Performs 2 datachecks: 1) Ensures that all transcripts with CDSs have exons,
         2) Ensures that all CDSs are within their exons boundaries. In both cases the script
         outputs a list of the offending transcripts and also write these to files. If the fix
         options are selected, the script returns a dataframe with the offending transcripts' type
         set to 'nontranslating_transcript'."""

        # Get a DF with all exons and CDS, group it by transcript and calculate useful metrics
        coding_transcripts = list(self.gff_df[self.gff_df["type"].isin(coding_transcript_types)]["ID"].unique())
        exons_cds_df = get_all_exon_and_cds_features(self.gff_df, transcripts=coding_transcripts)

        # Check if exons_cds_df is empty
        if exons_cds_df.empty:
            exit_with_error("No exons and CDS features were found. Please check your GFF file, \
                            maybe exons/CDSs do not have the right ID field connecting them to their parent transcripts.")
        
        counts_per_transcript_df = exons_cds_df.groupby('Parent').apply(_count_cds_exons_per_transcript)


        # Datacheck
        if dc_cds_but_no_exons:
            dc1_offending_transcripts=counts_per_transcript_df[(counts_per_transcript_df['no_of_exons']==0) & (counts_per_transcript_df['no_of_cds']>0)].index.values.tolist()
            if dc1_offending_transcripts:
                print_warning("The following transcripts have CDS but they don't have exons. These transcripts "
                              "will be written in cds_not_exons_transcripts.txt")
                print("\n".join(dc1_offending_transcripts))
                with open(prefix+'cds_not_exons_transcripts.txt', mode='wt') as myfile:
                    myfile.write('\n'.join(dc1_offending_transcripts))
                if dc_cds_but_no_exons_fix:
                    print_info("Changing the type of the offending transcripts to 'nontranslating_transcript'")
                    self.gff_df.loc[(self.gff_df['type'].isin(coding_transcript_types)) &
                                    (self.gff_df['ID'].isin(dc1_offending_transcripts)) , 'type'] = "nontranslating_transcript"
                    self.gff_df.drop(self.gff_df.loc[(self.gff_df['type'] == 'CDS') &
                                                     (self.gff_df['Parent'].isin(dc1_offending_transcripts))].index, inplace=True)

        # Datacheck
        if dc_cds_within_exons:
            dc2_offending_transcripts=counts_per_transcript_df[(counts_per_transcript_df['exons_max_coord'] < counts_per_transcript_df['cds_max_coord']) |
                                                 (counts_per_transcript_df['exons_min_coord'] > counts_per_transcript_df['cds_min_coord'])].index.values.tolist()
            if dc2_offending_transcripts:
                print_warning("The following transcripts have CDS which are not within exon boundaries. These transcripts "
                              "will be written in cds_not_within_exons_transcripts.txt")
                print("\n".join(dc2_offending_transcripts))
                with open(prefix+'cds_not_within_exons_transcripts.txt', mode='wt') as myfile:
                    myfile.write('\n'.join(dc2_offending_transcripts))
                if dc_cds_within_exons_fix:
                    print_info("Changing the type of the offending transcripts to 'nontranslating_transcript'")
                    self.gff_df.loc[(self.gff_df['type'].isin(coding_transcript_types)) &
                                    (self.gff_df['ID'].isin(dc2_offending_transcripts)), 'type'] = "nontranslating_transcript"
                    self.gff_df.drop(self.gff_df.loc[(self.gff_df['type'] == 'CDS') &
                                                     (self.gff_df['Parent'].isin(dc2_offending_transcripts))].index, inplace=True)

        # Datacheck
        if dc_coding_transcripts_with_cds:
            dc3_offending_transcripts=counts_per_transcript_df[counts_per_transcript_df['no_of_cds']==0].index.values.tolist()
            if dc3_offending_transcripts:
                print_warning("The following coding transcripts do not have CDS features. These transcripts "
                              "will be written in without_CDS_transcripts.txt")
                print("\n".join(dc3_offending_transcripts))
                with open(prefix+'without_CDS_transcripts.txt', mode='wt') as myfile:
                    myfile.write('\n'.join(dc3_offending_transcripts))
                if dc_coding_transcripts_with_cds_fix:
                    print_info("Changing the type of the offending transcripts to 'nontranslating_transcript'")
                    self.gff_df.loc[(self.gff_df['type'].isin(coding_transcript_types)) &
                                    (self.gff_df['ID'].isin(dc3_offending_transcripts)), 'type'] = "nontranslating_transcript"

# Updates attributes field and drops unecessary columns, then writes
# processed .gff to a specified output file.
def write_output_gff(gff_df, output_file):
    gff_df = finalise_attributes_column(gff_df)
    gff_df = make_coding_transcripts_mRNA(gff_df)
    gff_df = _drop_useless_columns(gff_df)
    gff_df = reorder_gff_features(gff_df)
    gff_df = _drop_superfluous_columns(gff_df)

    with open(output_file, 'w') as output_gff:
         output_gff.write("##gff-version 3\n")

    gff_df.to_csv(output_file, mode = 'a', header = False, sep = '\t', index = False)

# Unused functions that might inspire us in the future
# def link_exons_and_cds_features(gff_df):
#     exon_df = get_all_exon_features(gff_df)
#     exon_df["max_coord"] = exon_df[["start", "end"]].max(axis=1)
#     exon_df["min_coord"] = exon_df[["start", "end"]].min(axis=1)
#
#     cds_df = get_all_cds_features(gff_df)
#     cds_df["max_coord"] = cds_df[["start", "end"]].max(axis=1)
#     cds_df["min_coord"] = cds_df[["start", "end"]].min(axis=1)
#
#     outer_df = exon_df.merge(cds_df, how='inner', on='Parent', suffixes=("_exon","_cds"))
#     linking_df = outer_df[(outer_df["min_coord_cds"]>=outer_df["min_coord_exon"]) & (outer_df["max_coord_cds"]<=outer_df["max_coord_exon"])]
#
#     exon_df_with_cds = exon_df.merge(linking_df[["ID_exon","ID_cds"]], how='left', left_on="ID", right_on="ID_exon").drop("ID_exon", axis=1)
#     cds_df_with_exon = cds_df.merge(linking_df[["ID_exon", "ID_cds"]], how='left', left_on="ID", right_on="ID_cds").drop("ID_cds", axis=1)
#
#     return exon_df_with_cds, cds_df_with_exon
