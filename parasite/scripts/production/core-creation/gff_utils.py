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
        try:
            attribute_fields = row[9].split(';')
        except AttributeError:
            print(row)
            exit_with_error("Couldn't split the 9th column of the row above.")

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


# Reads in a .gff file as a Pandas dataframe.
def _load_gff(input_gff, gff_column_names):
    col_names = gff_column_names
    gff_df = pd.read_csv(input_gff, sep = '\t', comment = '#', names = col_names, index_col = False)
    return gff_df


# Parsing the loaded GFF pandas dataframe. Rows and values are not
# modified but the attributes column is parsed into separate columns,
# one for each identified attribute (e.g. ID, Name, Parent).
# Attributes not present in the original .gff file are NOT added.
def parse_gff(input_gff, gff_column_names): #, intron_types, exon_types, args):
    gff_df = _load_gff(input_gff, gff_column_names)

    print_info("Checking for redundant GFF types")
    redundant_types = dc_redundant_gff_types(gff_df)
    if redundant_types:
        if (intron_types in redundant_types) & (exon_types not in list(set(gff_df["type"]))):
            print_warning("Only introns and not exons exists in your GFF file.")
        print_info("Removing redundant GFF types: "+", ".join(redundant_types))
        gff_df = remove_redundant_gff_types(gff_df, redundant_types)
    
    gff_df = _parse_attributes_field(gff_df)
    gff_df = _make_ids_unique(gff_df)

    return gff_df

def parse_synonyms(synonyms_file):
    colnames = ["toplevel", "CommunityID", "INSDCID", "INSDC"]
    scaffold_df = pd.read_csv (synonyms_file, sep = '\t', names = colnames, usecols = ["CommunityID", "INSDCID"])
    # scaffold_df.set_index("INSDCID", inplace = True)
    return scaffold_df

def count_features(gff_df):
    """This function calculates all features of the GFF3 dataframe"""
    no_of_genes = len(get_all_gene_features(gff_df))
    no_of_tr = len(get_all_transcript_features(gff_df))
    no_of_exons = len(get_all_exon_features(gff_df))
    no_of_cds = len(get_all_cds_features(gff_df))
    counts = {"genes": no_of_genes, "transcripts": no_of_tr, "exons": no_of_exons, "cds": no_of_cds}
    return(counts)

def report_counts(counts_dict):
    report = "Count of features in the initial GFF:\t"
    for key, value in counts_dict.items():
        report += f"{key.capitalize()}: {value}\t"
    return(report)

def scaffolds_needs_renaming(gff_df, fasta):
    # If gff_df and fasta have the same scaffold names then no renaming is needed
    # Extract the unique scaffold values from the gff_df DataFrame
    unique_scaffolds_in_df = gff_df["scaffold"].unique()
    # Get the scaffold names from the fasta object
    scaffold_names_in_fasta = fasta.scaffold_names()
    # Check if all values from gff_df["scaffold"] exist in scaffold_names_in_fasta
    mismatching_scaffolds = [scaffold for scaffold in unique_scaffolds_in_df if scaffold not in scaffold_names_in_fasta]
    if len(mismatching_scaffolds) == 0:
        return False
    else:
        return mismatching_scaffolds

def gff_rename_scaffolds_to(gff_df, synonyms_df):
    # if gff_df scaffolds match scaffold_df["CommunityID"] then rename them to scaffold_df["INSDC"] if
    # they match scaffold_df["INSDC"] then rename them to scaffold_df["CommunityID"]
    if gff_df["scaffold"].isin(synonyms_df["CommunityID"]).all():
        return "INSDCID"
    elif gff_df["scaffold"].isin(synonyms_df["INSDCID"]).all():
        return("CommunityID")
    else:
        exit_with_error("GFF scaffolds do not match CommunityID or INSDCID columns in synonyms file. "\
                        "Cannot rename scaffolds.")

def rename_scaffolds_in_gff(gff_df, synonyms_df, rename_to):
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
        return(merged_df)

def rename_scaffolds_in_fasta(fasta, synonyms_df, rename_to):
    # Rename the scaffold names in the fasta object
    rename_from = [x for x in synonyms_df.columns.to_list() if x !=rename_to][0]
    print_info(f"Renaming fasta scaffolds from {rename_from} to {rename_to}")
    fasta.rename_scaffolds(synonyms_df, rename_from, rename_to)


# Replaces INSDCID scaffold names with community scaffold names using
# synonyms mapping provided in a synonyms.tsv file.
def rename_scaffolds(gff_df, synonyms_df, rename_gff_scaffolds_if_community):
    rename_to = gff_rename_scaffolds_to(gff_df, synonyms_df)
    if rename_to=="CommunityID":
        gff_df = rename_scaffolds_in_gff(gff_df, synonyms_df, rename_to)
    elif rename_to=="INSDCID":
        print_info("It looks like that the GFF3 file has correct "
                   "Community IDs but your FASTA file has the INSDC ids.")
        if rename_gff_scaffolds_if_community:
            gff_df = rename_scaffolds_in_gff(gff_df, synonyms_df, rename_to)
        else:
            print_info("Your GFF3 scaffold names are INSDC IDs but they are not matching with your fasta file. "
                            "If you would still want to rename the scaffold names of your "
                            "GFF3 file to INSDC IDs, please set the --rename_gff_scaffolds_even_if_community flag.")
    return(gff_df)



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
    parents_not_ids = list(gff_df[(~gff_df["Parent"].isin(gff_df["ID"].unique())) & (~gff_df["Parent"].isnull())]["Parent"].unique())
    if len(parents_not_ids) > 0:
        print_warning("These Parents in your GFF file are not IDs: " + ",".join(parents_not_ids))
    else:
        gff_df["ID"] = prefix + gff_df["ID"]
        gff_df["Parent"] = prefix + gff_df["Parent"]
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
    gene_df = gff_df.loc[gene_mask].copy(deep=True)
    return gene_df


# Returns dataframe consisting only of transcript features.
def get_all_transcript_features(gff_df):
    transcript_mask = gff_df.loc[:, "type"].isin(transcript_types)
    transcript_df = gff_df.loc[transcript_mask].copy(deep=True)
    return transcript_df


# Returns dataframe consisting only of exon features.
def get_all_exon_features(gff_df):
    exon_mask = gff_df.loc[:, "type"].isin(exon_types)
    exon_df = gff_df.loc[exon_mask].copy(deep=True)
    return exon_df

# Returns dataframe consisting only of exon features.
def get_all_intron_features(gff_df):
    intron_mask = gff_df.loc[:, "type"].isin(intron_types)
    intron_df = gff_df.loc[intron_mask].copy(deep=True)
    return intron_df

# Returns dataframe consisting only of CDS features.
def get_all_cds_features(gff_df):
    cds_mask = gff_df.loc[:, "type"].isin(cds_types)
    cds_df = gff_df.loc[cds_mask].copy(deep=True)
    return cds_df

# Returns a dataframe containing only exon and CDS features.
def get_all_exon_and_cds_features(gff_df, transcripts=None):
    if transcripts is None:
        exon_cds_df = gff_df[(gff_df["type"] == "CDS") | (gff_df["type"] == "exon")].copy(deep=True)
    else:
        exon_cds_df = gff_df[((gff_df["type"] == "CDS") | (gff_df["type"] == "exon")) & (gff_df["Parent"].isin(transcripts))].copy(deep=True)
    return exon_cds_df

# Returns a dataframe containing transcript, exon + CDS + transcript features.
def get_all_exon_cds_and_transcript_features(gff_df, transcripts=None):
    if transcripts is None:
        nongenes_df = gff_df[(gff_df["type"].isin(cds_types) | (gff_df["type"].isin(exon_types)) | (gff_df["type"].isin(transcript_types)))].copy(deep=True)
    else:
        nongenes_df = gff_df[(gff_df["type"].isin(cds_types) | (gff_df["type"].isin(exon_types)) | (gff_df["type"].isin(transcript_types))) & (gff_df["Parent"].isin(transcripts))].copy(deep=True)
    return nongenes_df

# This function takes a dataframe grouped by the "Parent" column and
# returns a pandas Series with the exon/cds metrics below.
# This function is being used by the transcripts_have_exons_and_valid_cds datacheck.
def _count_cds_exons_per_transcript(gff_df):
    exons_df = gff_df[gff_df['type'].isin(exon_types)]
    cds_df = gff_df[gff_df['type'].isin(cds_types)]
    exons_cds_stats = {
        'no_of_exons': exons_df['Parent'].count(),
        'no_of_cds': cds_df['Parent'].count(),
        'exons_max_coord': max([exons_df["start"].max(),exons_df["end"].max()]),
        'exons_min_coord': min([exons_df["start"].min(),exons_df["end"].min()]),
        'cds_max_coord': max([cds_df["start"].max(),cds_df["end"].max()]),
        'cds_min_coord': min([cds_df["start"].min(),cds_df["end"].min()]),
    }
    return pd.Series(exons_cds_stats)



def get_parent_gene_for_all(gff_df):
    """
    This function creates a Parent_Gene column for all features.
    This column contains the ID of the feature's parent gene. 
    """
    count_features_before = len(gff_df[gff_df["type"].isin(transcript_types+exon_types+cds_types)])
    all_genes = set(gff_df[(gff_df["type"].isin(gene_types))]["ID"].to_list())
    
    id_parents_df = gff_df[(~gff_df["type"].isin(gene_types)) & (gff_df["type"].isin(allowed_types))][["ID","Parent"]]
    gene_parents_df = gff_df[(gff_df["type"].isin(gene_types))][["ID"]]
    
    # Check if there are any features without a parent. Exit if there are
    if id_parents_df["Parent"].isnull().any():
        parentless_features = list(set(id_parents_df[~id_parents_df["Parent"].isnull()]["ID"]))
        exit_with_error("There are some transcripts/exons/CDSs without a Parent. These will be written at "
                        "parentless_features.txt. Please check the GFF file, correct it and re-run.")
        with open(os.path.join('parentless_features.txt'), mode='wt') as myfile:
            myfile.write('\n'.join(parentless_features))
    
    
    id_parents_dict = dict(zip(id_parents_df['ID'], id_parents_df['Parent']))
    gene_parents_dict = dict(zip(gene_parents_df["ID"], gene_parents_df["ID"]))
    id_parents_dict.update(gene_parents_dict)

    # Add the Parent_Gene ID
    gff_df['Parent_Gene'] = gff_df['Parent'].map(id_parents_dict)

    # Do it for the genes as well
    gene_rows = gff_df[gff_df['type'].isin(gene_types)]
    gff_df.loc[gene_rows.index, 'Parent_Gene'] = gene_rows['ID']

    # Check if there are any 'Parent_Gene' IDs that do not belong to the all_genes. Exit if there are
    if not gff_df[(~gff_df["Parent_Gene"].isnull())]["Parent_Gene"].isin(all_genes).all()==True:
        # print failing rows
        exit_with_error("Couldn't sucessfully assign Parent_Gene column. Please check the GFF file, correct it and re-run.")
    
    # Check if there are any features without a Parent_Gene. Exit if there are.
    if gff_df[gff_df["type"].isin(gene_types+transcript_types+exon_types+cds_types)]["Parent_Gene"].isnull().any():
        exit_with_error("Couldn't sucessfully assign Parent_Gene column. Please check the GFF file, correct it and re-run.")

    # Make sure no features are deleted:
    count_features_after = len(gff_df[gff_df["type"].isin(transcript_types+exon_types+cds_types)])
    if count_features_after != count_features_before:
        exit_with_error("Couldn't sucessfully assign Parent_Gene column. Please check the GFF file, correct it and re-run.")

    return(gff_df)

def all_exon_starts_less_than_ends(gff_df):
    exon_df = get_all_exon_features(gff_df)
    # Check if there are any rows where the start is greater than the end. Exit if there are
    if exon_df[(exon_df["start"] > exon_df["end"])]["start"].any():
        return(False)
    else:
        return(True)

def merge_adjacent_exons(gff_df):
    exon_df = get_all_exon_features(gff_df)
    df = exon_df
    newexon_df = df.sort_values(["Parent","start"]).groupby([df.Parent,((df.end.shift()+1)-df.start).lt(0).cumsum()]).agg(
        {k: "first" if k not in ["start","end"] else "min" if k == "start" else "max" for k in df.columns}).reset_index(drop=True)
    new_gff_df = pd.concat([gff_df[~gff_df["type"].isin(exon_types)],newexon_df])
    return(new_gff_df)

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
    
    # Remove any single/double quotation from the final attributes column
    gff_df["attributes"] = gff_df["attributes"].str.replace("'", "")
    gff_df["attributes"] = gff_df["attributes"].str.replace('"', '')
    
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
def switch_transcripts_to_pseudogenic_transcripts(gff_df):
    pseudogenes = list(set(gff_df[gff_df["type"].isin(pseudogene_types)]["ID"].to_list()))

    pseudogenic_transcript_mask = gff_df["type"].isin(transcript_types) & gff_df["Parent"].isin(pseudogenes)

    # Rename all the gff_df["type"] values to "pseudogenic_transcript" in for pseudogenic_transcript_mask 
    gff_df.loc[pseudogenic_transcript_mask, "type"] = "pseudogenic_transcript"

    return(gff_df)

# Creates new transcript features for pseudogenes
def create_pseudogenic_transcripts(gff_df, gene_ids):
    pseudogenes_df = gff_df.loc[gff_df["type"].isin(pseudogene_types) & gff_df["ID"].isin(gene_ids)].copy(deep=True)
    pseudogenes_df["type"] = "pseudogenic_transcript"
    pseudogenes_df["Parent"] = pseudogenes_df["ID"]
    pseudogenes_df["ID"] = pseudogenes_df["ID"] + "_transcript"
    gff_df = pd.concat([gff_df, pseudogenes_df])
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
    def rename_scaffold_names(self, synonyms_df):
        """Renames the FASTA scaffold names based on the synonyms dataframe"""
        synonyms_dict = {
            old_scaffold_name: new_scaffold_name
            for old_scaffold_name, new_scaffold_name in synonyms_df.itertuples()
        }
        for record in SeqIO.parse(self.path, "fasta"):
            for record in SeqIO.parse(self.path, "fasta"):
                record.id = synonyms_dict[record.id]
   
def dc_genes_have_names(gff_df):
    # Filtering the DataFrame
    gene_df = get_all_gene_features(gff_df)
    # If gene_df doesn't have a "Name" column, exit with error:
    if not "Name" in gene_df.columns or gene_df["Name"].isna().any():
        exit_with_error("Some/All genes do not have a Name column. Please add a Name column to your GFF file or use the --infer_gene_name/-n option.")

def deprecated_counts_per_transcript(gff_df):
    """Calculates the number of exons and CDS per transcript and returns it as a dataframe"""
    # Get a DF with all exons and CDS, group it by transcript and calculate useful metrics
    coding_transcripts = set(list(gff_df.loc[gff_df["type"].isin(coding_transcript_types), "ID"].unique()))
    exons_cds_df = gff_df.loc[(gff_df["type"].isin(exon_types + cds_types)) & (gff_df["Parent"].isin(coding_transcripts))]

    # Check if exons_cds_df is empty
    if exons_cds_df.empty:
        exit_with_error("No exons and CDS features were found. Please check your GFF file, \
                        maybe exons/CDSs do not have the right ID field connecting them to their parent transcripts.")
    counts_per_transcript_df = exons_cds_df.groupby('Parent').apply(_count_cds_exons_per_transcript)
    return(counts_per_transcript_df)

def counts_per_transcript(gff_df):
    """Calculates the number of exons and CDS, and their min/max coordinates per transcript and returns it as a dataframe"""
    exons_cds_df = get_all_exon_and_cds_features(gff_df)
    count_df=exons_cds_df.groupby(['Parent','type']).aggregate(
        {"type":"count","start":["min","max"], "end":["min","max"]}).unstack(fill_value=0).stack().reset_index()
    count_df=count_df[count_df["Parent"].isin(get_all_transcript_features(gff_df)["ID"])]
    count_df.columns = [y.strip("_") for y in ["_".join(x) for x in list(count_df.columns)]]
    pivot_counts_df = count_df.pivot_table(index="Parent",columns="type")
    pivot_counts_df.columns = [y.strip("_") for y in ["_".join(x) for x in list(pivot_counts_df.columns)]]
    pivot_counts_df["cds_max_coord"] = pivot_counts_df[["end_max_CDS","start_max_CDS"]].max(axis=1)
    pivot_counts_df["cds_min_coord"] = pivot_counts_df[["end_min_CDS","start_min_CDS"]].min(axis=1)
    pivot_counts_df["exons_max_coord"] = pivot_counts_df[["end_max_exon","start_max_exon"]].max(axis=1)
    pivot_counts_df["exons_min_coord"] = pivot_counts_df[["end_min_exon","start_min_exon"]].min(axis=1)
    return(pivot_counts_df)

def counts_per_gene(gff_df):
    not_genes_df = get_all_exon_cds_and_transcript_features(gff_df)
    count_df=not_genes_df.groupby(['Parent_Gene','Parent','type']).aggregate(
        {"type":"count","start":["min","max"], "end":["min","max"]}).unstack(fill_value=0).stack().reset_index()


def dc_transcripts_parentage(gff_df, outprefix):
    """Transcripts should have genes as Parents.
    Anything else is wrong"""
    all_genes = set(list(get_all_gene_features(gff_df)["ID"].unique()))
    tr_df = get_all_transcript_features(gff_df)
    wrong_parent_df = tr_df.loc[~tr_df["Parent"].isin(all_genes)]
    offending_transcripts = list(set(wrong_parent_df["ID"]))
    if offending_transcripts:
        print_warning("There are transcripts without genes as Parents. These transcripts will "
                      "be written in the transcripts_with_wrong_parents.txt")
        with open(os.path.join(outprefix,'transcripts_with_wrong_parents.txt'), mode='wt') as myfile:
            myfile.write('\n'.join(offending_transcripts))
        return offending_transcripts
    else:
        return False

def dc_exons_parentage(gff_df, outprefix):
    """CDS/exon should have transcripts as Parents.
    Anything else is wrong"""
    all_transcripts = set(list(get_all_transcript_features(gff_df)["ID"].unique()))
    exon_cds_df = get_all_exon_features(gff_df)
    wrong_parent_df = exon_cds_df.loc[~exon_cds_df["Parent"].isin(all_transcripts)]
    offending_features = list(set(wrong_parent_df["ID"])) 
    if offending_features:
        print_warning("There are exons/CDSs without transcripts as Parents. These features will "
                      "be written in the exons_or_CDS_with_wrong_parents.txt")
        with open(os.path.join(outprefix,'exons_or_CDS_with_wrong_parents.txt'), mode='wt') as myfile:
            myfile.write('\n'.join(offending_features))
        return offending_features
    else:
        return False

def dc_cds_parentage(gff_df, outprefix):
    """CDS/exon should have transcripts as Parents.
    Anything else is wrong"""
    all_transcripts = set(list(get_all_transcript_features(gff_df)["ID"].unique()))
    exon_cds_df = get_all_cds_features(gff_df)
    wrong_parent_df = exon_cds_df.loc[~exon_cds_df["Parent"].isin(all_transcripts)]
    offending_features = list(set(wrong_parent_df["ID"]))
    if offending_features:
        print_warning("There are CDSs without transcripts as Parents. These features will "
                      "be written in the exons_or_CDS_with_wrong_parents.txt")
        with open(os.path.join(outprefix,'exons_or_CDS_with_wrong_parents.txt'), mode='wt') as myfile:
            myfile.write('\n'.join(offending_features))
        return offending_features
    else:
        return False

def dc_cds_but_no_exons(gff_df, outprefix):
    """Checks if there are CDSs but no exons in your gff_df. It returns a list with the offending transcripts"""

    counts_per_transcript_df = counts_per_transcript(gff_df)

    dc1_offending_transcripts=set(list(counts_per_transcript_df[(counts_per_transcript_df['type_count_exon']==0) & (counts_per_transcript_df['type_count_CDS']>0)].index.values.tolist()))
    if dc1_offending_transcripts:
        print_warning("There are transcripts with CDS but without exons. These transcripts "
                            "will be written in cds_not_exons_transcripts.txt")
        with open(outprefix+'cds_not_exons_transcripts.txt', mode='wt') as myfile:
            myfile.write('\n'.join(dc1_offending_transcripts))
        return(dc1_offending_transcripts)
    else:
        return False

def dc_exons_but_no_cds(gff_df, outprefix):
    """Checks if there are exons but no CDS in your gff_df. It returns a list with the offending transcripts"""

    counts_per_transcript_df = counts_per_transcript(gff_df)

    dc1_offending_transcripts=set(list(counts_per_transcript_df[(counts_per_transcript_df['type_count_exon']>0) & (counts_per_transcript_df['type_count_CDS']==0)].index.values.tolist()))
    if dc1_offending_transcripts:
        print_warning("There are transcripts with exons but without CDS. These transcripts "
                            "will be written in cds_not_exons_transcripts.txt")
        with open(outprefix+'cds_not_exons_transcripts.txt', mode='wt') as myfile:
            myfile.write('\n'.join(dc1_offending_transcripts))
        return(dc1_offending_transcripts)
    else:
        return False

def dc_coding_transcripts_with_exons_cds(gff_df, outprefix):
    """Checks if there are CDSs and exons for coding transcripts in your gff_df. It returns a list with the offending transcripts"""
    counts_per_transcript_df = counts_per_transcript(gff_df)
    non_coding_transcripts=set(list(gff_df.loc[~gff_df["type"].isin(pseudogene_types),"ID"]))
    transcripts_df=get_all_transcript_features(gff_df)
    dc1_offending_transcripts=list(set(counts_per_transcript_df[(counts_per_transcript_df['type_count_exon']==0) & (counts_per_transcript_df['type_count_CDS']==0)].index.values.tolist()))
    dc1_offending_transcripts+=list(set(transcripts_df[~transcripts_df["ID"].isin(counts_per_transcript_df.index)]["ID"]))
    dc1_offending_transcripts = [x for x in dc1_offending_transcripts if x not in non_coding_transcripts]
    if dc1_offending_transcripts:
        print_warning("There are transcripts with CDS but without exons. These transcripts "
                            "will be written in cds_not_exons_transcripts.txt")
        with open(outprefix+'cds_not_exons_transcripts.txt', mode='wt') as myfile:
            myfile.write('\n'.join(dc1_offending_transcripts))
        return(dc1_offending_transcripts)
    else:
        return False


def transcripts_to_nontranslating_transcripts(gff_df, offending_transcripts):
    """This function changes the type of the offending transcripts to 'nontranslating_transcript'"""
    gff_df.loc[(gff_df['type'].isin(coding_transcript_types)) &
                                (gff_df['ID'].isin(offending_transcripts)) , 'type'] = "nontranslating_transcript"
    gff_df.drop(gff_df.loc[(gff_df['type'] == 'CDS') &
                            (gff_df['Parent'].isin(offending_transcripts))].index, inplace=True)
    return(gff_df)
    
def dc_cds_within_exons(gff_df, outprefix):
    """Checks if there are transcripts with CDSs exceeding the exon boundaries. It returns a list with the offending transcripts"""
    counts_per_transcript_df = counts_per_transcript(gff_df)
    dc2_offending_transcripts=list(set(counts_per_transcript_df[ ((counts_per_transcript_df['type_count_CDS'] > 0) & (counts_per_transcript_df['type_count_exon'] > 0)) &
                                                       ((counts_per_transcript_df['exons_max_coord'] < counts_per_transcript_df['cds_max_coord']) |
                                                        (counts_per_transcript_df['exons_min_coord'] > counts_per_transcript_df['cds_min_coord']))].index.values.tolist()))
    if dc2_offending_transcripts:
        print_warning("There are transcripts with CDSs within their exons. These transcripts "
                            "will be written in cds_within_exons_transcripts.txt")
        with open(outprefix+'cds_within_exons_transcripts.txt', mode='wt') as myfile:
            myfile.write('\n'.join(dc2_offending_transcripts))
        return(dc2_offending_transcripts)
    else:
        return False

def dc_coding_transcripts_with_cds(gff_df, outprefix):
    cds_parents = set(list(gff_df[gff_df["type"].isin(cds_types)]["Parent"].unique()))
    coding_transcripts = set(list(gff_df[gff_df["type"].isin(coding_transcript_types)]["ID"].unique()))
    dc3_offending_transcripts = [x for x in coding_transcripts if x not in cds_parents]
    if dc3_offending_transcripts:
        print_warning("There are coding transcripts without any CDSs. "\
                      "These will be written in cds_not_exons_transcripts.txt. "
                      "If these are pseudogenes please mark them as such and re-run the parser.")
        with open(outprefix+'cds_not_exons_transcripts.txt', mode='wt') as myfile:
            myfile.write('\n'.join(dc3_offending_transcripts))
        return(dc3_offending_transcripts)
    else:
        return False

def dc_genes_without_transcripts(gff_df, outprefix):
    transcript_parents = set(list(gff_df[gff_df["type"].isin(transcript_types)]["Parent"].unique()))
    all_genes = list(gff_df[gff_df["type"].isin(gene_types)]["ID"].unique())
    genes_without_transcripts = [gene for gene in all_genes if gene not in transcript_parents]
    if genes_without_transcripts:
        print_warning("Some genes in your GFF file do not have any associated transcripts."
                    "These genes will be written in genes_without_transcripts.txt")
        with open(outprefix+'genes_without_transcripts.txt', mode='wt') as myfile:
                myfile.write('\n'.join(genes_without_transcripts))
        return(genes_without_transcripts)
    else:
        return False

def remove_redundant_gff_types(gff_df, redundant_types):
    gff_df = gff_df[~gff_df["type"].isin(redundant_types)]
    return(gff_df)

def dc_redundant_gff_types(gff_df):
    all_types = list(set(gff_df["type"]))
    redundant_types = [x for x  in all_types if x not in allowed_types]
    if redundant_types:
        print_warning("These types in your GFF files are not supported by this parser:")
        print_w_indent("\n\t".join(redundant_types))
        print_w_indent("These types will be ignored while processing it. If you do not want to "                      "and they will be ignored while processing it. If you do not want to "
                      "ignore these types please add them to their relevant list in the "
                      "reformat_gff_dependencies.py.")
        return(redundant_types)
    else:
        return False

def extrapolate_transcripts_from_gene_coordinates(gff_df, gene_ids, transcript_type, transcript_id_suffix):
    """
    This function extrapolates transcript features from gene features.
    It is used when the genes_without_transcripts datacheck fails.
    """
    # Extrapolate transcript features from gene features
    off_genes_df = gff_df[gff_df["type"].isin(gene_types) & gff_df["ID"].isin(gene_ids)]
    extrapolated_transcripts = off_genes_df.copy(deep=True).reset_index(drop=True)
    extrapolated_transcripts["type"] = transcript_type
    # Set transcript's parent to gene it is based on.
    extrapolated_transcripts["Parent"] = extrapolated_transcripts["ID"]
    # Modify transcript ID to clarify it is mRNA.
    extrapolated_transcripts["ID"] += transcript_id_suffix
    extrapolated_transcripts["Name"] = extrapolated_transcripts["ID"]
    # Add extrapolated transcripts to original gff dataframe
    cor_df = pd.concat([gff_df, extrapolated_transcripts]).reset_index()
    # gff_df row IDs with type in exon_types or cds_types that their Parent isin genes_without_transcripts 
    # should be renamed to Parent += suffix
    cor_df.loc[cor_df["type"].isin(exon_types + cds_types) & cor_df["Parent"].isin(gene_ids), "Parent"] += transcript_id_suffix
    if not cor_df.loc[cor_df["type"].isin(exon_types + cds_types) & cor_df["Parent_Gene"].isin(gene_ids)]["Parent"].isin(extrapolated_transcripts["ID"]).all():
        print_info("Couldn't extrapolate transcipts from genes correctly. Exiting...")
    return cor_df

def extrapolate_features_from_transcript_coordinates(gff_df, transcript_ids, feat_type, feat_id_suffix):
    """
    This function extrapolates transcript features from gene features.
    It is used when the genes_without_transcripts datacheck fails.
    """
    # Extrapolate transcript features from gene features
    off_df = gff_df[gff_df["type"].isin(transcript_types) & gff_df["ID"].isin(transcript_ids)]
    extrapolated_features = off_df.copy(deep=True).reset_index(drop=True)
    extrapolated_features["type"] = feat_type
    # Set transcript's parent to gene it is based on.
    extrapolated_features["Parent"] = extrapolated_features["ID"]
    # Modify transcript ID to clarify it is mRNA.
    extrapolated_features["ID"] += feat_id_suffix
    # Add extrapolated transcripts to original gff dataframe
    cor_df = pd.concat([gff_df, extrapolated_features]).reset_index()
    return cor_df

def extrapolate_CDS_from_exons_for_transcripts(gff_df, transcript_ids):
    """
    This function extrapolates transcript features from gene features.
    It is used when the genes_without_transcripts datacheck fails.
    """
    # Extrapolate transcript features from gene features
    off_df = gff_df[gff_df["type"].isin(exon_types) & gff_df["Parent"].isin(transcript_ids)]
    extrapolated_features = off_df.copy(deep=True).reset_index(drop=True)
    extrapolated_features["type"] = "CDS"
    # Modify transcript ID to clarify it is mRNA.
    extrapolated_features["ID"] += "_CDS"
    # Add extrapolated transcripts to original gff dataframe
    cor_df = pd.concat([gff_df, extrapolated_features]).reset_index(drop=True)
    return cor_df

def extrapolate_exons_from_CDS_for_transcripts(gff_df, transcript_ids):
    """
    This function extrapolates transcript features from gene features.
    It is used when the genes_without_transcripts datacheck fails.
    """
    # Extrapolate transcript features from gene features
    off_df = gff_df[gff_df["type"].isin(cds_types) & gff_df["Parent"].isin(transcript_ids)]
    extrapolated_features = off_df.copy(deep=True).reset_index(drop=True)
    extrapolated_features["type"] = "exon"
    # Modify transcript ID to clarify it is mRNA.
    extrapolated_features["ID"] += "_exon"
    # Add extrapolated transcripts to original gff dataframe
    cor_df = pd.concat([gff_df, extrapolated_features]).reset_index()
    return cor_df

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

def perform_datachecks(gff_df, args, outprefix):

    print_info("Datacheck: Checking if all genes have names")
    dc_genes_have_names(gff_df)

    print_info("Datacheck: Checking if all transcripts have genes as their parents")
    transcripts_with_wrong_parents = dc_transcripts_parentage(gff_df, outprefix)
    if transcripts_with_wrong_parents:
        exit_with_error("Review your GFF and re-run.")

    print_info("Datacheck: Checking if all genes have transcripts")
    genes_without_transcripts = dc_genes_without_transcripts(gff_df, outprefix)
    if genes_without_transcripts:
        if args.dc_remove_genes_without_transcripts:
            print_info("Removing genes without transcripts and all their features.")
            gff_df = gff_df[~(gff_df["Parent_Gene"].isin(genes_without_transcripts))]
        elif args.dc_turn_genes_without_transcripts_into_pseudogenes:
            print_info("Turning genes without transcripts into pseudo-genes.")
            gff_df.loc[(gff_df["type"].isin(gene_types) & gff_df["ID"].isin(genes_without_transcripts)),"type"] = "pseudogene"
            gff_df = create_pseudogenic_transcripts(gff_df, genes_without_transcripts)
            gff_df = pseudogenes_to_genes(gff_df)
        else:
            if args.dc_genes_without_transcripts_fix:
                print_info("Extrapolating transcripts from gene coordinates for genes without transcripts.")
                gff_df = extrapolate_transcripts_from_gene_coordinates(gff_df, genes_without_transcripts, "mRNA", "_mRNA")
                if args.dc_genes_without_transcripts_into_nontranslating:
                    print_info("Changing the type of extrapolated transcripts to 'nontranslating_transcript'")
                    gff_df.loc[gff_df["type"].isin(transcript_types) & gff_df["Parent"].isin(genes_without_transcripts),"type"] = "nontranslating_transcript"
                if dc_genes_without_transcripts(gff_df, outprefix):
                    exit_with_error("Couldn't extrapolate transcipts for all failing genes correctly. Exiting...")
            else:
                print_warning("Detected genes without transcripts but the --no_fix_genes_without_transcripts_dc flag is set."\
                    "The script will continue but the missing transcripts will not be extrapolated.")
    
    if gff_df["type"].isin(pseudogene_types).any():
        if pseudogenes_with_CDS(gff_df) and  not args.keep_pseudogenic_CDSs:
            print_info("Found pseudogenes with CDSs. Removing pseudogenes with CDSs. If you would like to keep them " + 
                       "use the keep_pseudogenic_CDSs flag")
            gff_df = remove_pseudogenic_CDS(gff_df)
        print_info("Switching the type of transcripts of pseudogenes to pseudogenic_transcripts")
        gff_df = switch_transcripts_to_pseudogenic_transcripts(gff_df)
        print_info("Switching the type of \"pseudogene\" genes to \"gene\"")
        gff_df = pseudogenes_to_genes(gff_df)


    print_info("Datacheck: Checking for transcripts without any exons and CDSs.")
    transcripts_without_exons_and_CDSs = dc_coding_transcripts_with_exons_cds(gff_df, outprefix)
    if transcripts_without_exons_and_CDSs:
        if args.turn_coding_transcripts_without_exons_and_cds_to_nontranslating:
            print_info("Changing the type of the transcripts without exons and CDSs to 'nontranslating_transcript'")
            gff_df = transcripts_to_nontranslating_transcripts(gff_df, transcripts_without_exons_and_CDSs)
        elif args.extrapolate_missing_exons_and_CDS:
            print_info("Extrapolating missing exons and CDSs from transcript coordinates")
            gff_df = extrapolate_features_from_transcript_coordinates(gff_df, transcripts_without_exons_and_CDSs, "exon", "_exon")
            gff_df = extrapolate_CDS_from_exons_for_transcripts(gff_df, transcripts_without_exons_and_CDSs) 
        else:
            if args.remove_coding_transcripts_without_exons_and_cds:
                print_info("Removing transcripts without exons and CDSs")
                gff_df = gff_df[~(gff_df["type"].isin(coding_transcript_types) & gff_df["ID"].isin(transcripts_without_exons_and_CDSs))]
            else: print_warning("The --do_not_fix_coding_transcripts_without_exons_and_cds flag is set."\
                            "The script will continue but the offending transcripts will not be removed.")

    print_info("Datacheck: Checking for transcript with CDS but no exons")
    transcripts_with_CDS_but_no_exons = dc_cds_but_no_exons(gff_df, outprefix)
    if transcripts_with_CDS_but_no_exons:
        if args.turn_coding_transcripts_without_exons_but_with_cds_to_nontranslating:
            print_info("Changing the type of the transcripts with CDSs but without exons to 'nontranslating_transcript'")
            gff_df = transcripts_to_nontranslating_transcripts(gff_df, transcripts_with_CDS_but_no_exons)
        elif args.extrapolate_missing_exons_from_CDS:
            print_info("Extrapolate missing exons from CDS coordinates")
            gff_df = extrapolate_exons_from_CDS_for_transcripts(gff_df, transcripts_with_CDS_but_no_exons)
        else: print_warning("The --do_not_fix_coding_transcripts_without_exons_but_with_cds flag is set."\
                            "The script will continue but the offending transcripts will not be changed to 'nontranslating_transcript'")


    print_info("Datacheck: Checking for transcript with exons but without CDS")
    transcripts_with_exons_but_no_CDS = dc_exons_but_no_cds(gff_df, outprefix)
    if transcripts_with_exons_but_no_CDS:
        if args.turn_coding_transcripts_without_CDS_but_with_exons_to_nontranslating:
            print_info("Changing the type of the transcripts with exons but without CDSs to 'nontranslating_transcript'")
            gff_df = transcripts_to_nontranslating_transcripts(gff_df, transcripts_with_exons_but_no_CDS)
        elif args.extrapolate_missing_CDS_from_exons:
            print_info("Extrapolate missing CDS from exons coordinates")
            gff_df = extrapolate_CDS_from_exons_for_transcripts(gff_df, transcripts_with_exons_but_no_CDS)
        else: print_warning("The --do_not_fix_coding_transcripts_with_exons_but_without_CDS flag is set."\
                            "The script will continue but the offending transcripts will not be changed to 'nontranslating_transcript'")
            

    print_info("Datacheck: Checking if all exons have transcript as their parents")
    exons_with_wrong_parents = dc_exons_parentage(gff_df, outprefix)
    if exons_with_wrong_parents:
        if args.dc_exons_parentage_fix:
            print_info("Removing exons with wrong parentage")
            # Remove rows that have "type" in exon_types and "ID" in exons_with_wrong_parents
            gff_df = gff_df[~(gff_df["type"].isin(exon_types) & gff_df["ID"].isin(exons_with_wrong_parents))]
            if dc_exons_parentage(gff_df, outprefix):
                exit_with_error("Couldn't remove exons with wrong parentage. Exiting...")
        else: print_warning("Detected exons with wrong parentage but the --no_fix_exons_parentage flag is set."\
                            "The script will continue but the exons will not be removed.")
            
    print_info("Datacheck: Checking if all CDSs have transcript as their parents")
    cds_with_wrong_parents = dc_cds_parentage(gff_df, outprefix)
    if cds_with_wrong_parents:
        if args.dc_cds_parentage_fix:
            print_info("Removing CDSs with wrong parentage")
            # Remove rows that have "type" in cds_types and "ID" in cds_with_wrong_parents
            gff_df = gff_df[~(gff_df["type"].isin(cds_types) & gff_df["ID"].isin(cds_with_wrong_parents))]
            if dc_cds_parentage(gff_df, outprefix):
                exit_with_error("Couldn't remove CDSs with wrong parentage. Exiting...")
        else: print_warning("Detected CDSs with wrong parentage but the --no_fix_cds_parentage flag is set."\
                            "The script will continue but the CDSs will not be removed.")
    

    print_info("Datacheck: Checking for wrong CDS boundaries")
    transcripts_with_wrong_CDS_boundaries = dc_cds_within_exons(gff_df, outprefix)
    if transcripts_with_wrong_CDS_boundaries:
        if args.dc_cds_within_exons_fix:
            print_info("Changing the type of the transcripts with CDSs within exons to 'nontranslating_transcript'")
            gff_df = transcripts_to_nontranslating_transcripts(gff_df, transcripts_with_wrong_CDS_boundaries)
        else: print_warning("The --no_fix_coding_transcripts_with_exons_cds_dc flag is set."\
                            "The script will continue but the offending transcripts will not be removed.")

    print_info("Datacheck: Checking for coding transcripts without CDSs.")
    transcripts_without_CDSs = dc_coding_transcripts_with_cds(gff_df, outprefix)
    if transcripts_without_CDSs:
        if args.dc_coding_transcripts_with_cds_fix:
            print_info("Changing the type of the transcripts without CDSs to 'nontranslating_transcript'")
            gff_df = transcripts_to_nontranslating_transcripts(gff_df, transcripts_without_CDSs)
        else: print_warning("The --no_fix_coding_transcripts_with_cds_dc flag is set."\
                            "The script will continue but the offending transcripts will not be changed to 'nontranslating_transcript'")
    
    print("end of DataChecks")
    return gff_df

