#!/hps/software/users/wormbase/parasite/shared/.pyenv/versions/p395/bin/python
from gff_utils import *
from Bio import SeqIO
from reformat_gff_dependencies import *
from ProductionUtils import *
import re
import sys
import glob
from argparse import ArgumentParser

pd.options.mode.chained_assignment = None  # default='warn'


# Grabs and return command line arguments. Only input.gff and output
# destination are required arguments. -f and -s options allow for 
# specifying a  scaffolds .fasta file and a .tsv mapping NCBI accessions 
# to community scaffold names, but are not required. Other options 
# specify what processing steps should be taken (will vary) 
def get_args():
    parser = ArgumentParser()
    parser.add_argument("input_gff",
        help = "Unprocessed .gff file which requires parsing.")
    parser.add_argument("output_gff",
        help = "Output path to write modified .gff file to.")
    # source_to_WB is set to True by default
    parser.add_argument("--keep_source", action = "store_false", dest = "source_to_WB",
        help = "Leave the source column unchanged. By default, it is changed to WormBase_imported")
    # scaffold_rename is set to True by default.
    parser.add_argument("--keep_scaffold_names", action = "store_false", dest = "scaffold_rename",
        help = "Leave scaffold names unchanged. By default, they are mapped to an alternative name in the synonyms file and renamed.")
    # has_parentage is set to True by default.
    parser.add_argument("--no_parentage", action = "store_false", dest = "has_parentage",
        help = "Indicate input .gff file has no parentage information to link features. By default, assumes parentage information is included.")
    parser.add_argument("-g", "--gene_prefix", default = False,
        help = "A string to append to each gene and ID name, in case original names are lacking. No changes made by default.")

    parser.add_argument("-p", "--prefixes", required = False, default = None, nargs = '+',
        help = "List of prefixes to remove from .gff fields (e.g. parts of gene IDs that are constant across all genes).")
    parser.add_argument("--name_prefixes", required=False, default=None, nargs='+',
                        help="List of prefixes to remove from .gff fields (e.g. parts of gene Name that are constant across all genes).")

    # Set of arguments to extrapolate missing features for a given scaffold.

    # Extrapolate transcripts from genes, or vice versa.
    gene_transcript_extrapolation = parser.add_mutually_exclusive_group()
    gene_transcript_extrapolation.add_argument("-T", "--extrapolate_transcripts_for_scaffold", required = False, default = None,
        help = "Use gene information to create missing transcript features for a specified scaffold (provide community scaffold name).")
    gene_transcript_extrapolation.add_argument("-G", "--extrapolate_genes_for_scaffold", required = False, default = None,
        help = "Use transcript information to create missing gene features for a specified scaffold (provide community scaffold name).")

    # Extrapolate exons from CDSs, or vice versa
    exon_cds_extrapolation = parser.add_mutually_exclusive_group()
    exon_cds_extrapolation.add_argument("-E", "--extrapolate_exons_for_scaffold", required = False, default = None,
        help = "Use CDS information to create missing exon features for a specified scaffold (provide community scaffold name).")
    exon_cds_extrapolation.add_argument("-C", "--extrapolate_CDSs_for_scaffold", required = False, default = None,
        help = "Use exon information to create missing CDS features for a specified scaffold (provide community scaffold name).")

    name_inference_group = parser.add_mutually_exclusive_group()
    name_inference_group.add_argument("-n", "--infer_gene_names", default = False, action = "store_true",
        help = "Use feature ID attribute to create Name attribues for features that lack a name (directy copy)")
    name_inference_group.add_argument("-N", "--overwrite_gene_names", default = False, action = "store_true",
        help = "Use feature ID attribute to create Name attribues for all features (directy copy). Overwrites original names.")

    parser.add_argument("-f", "--fasta", required = False, default = None, type = str,
        help = ".fasta file containing assembly scaffolds. If not specified, .fasta file will be inferred from directory contents.")
    parser.add_argument("-s", "--synonyms_file", required = False, default = None, type = str,
        help = ".tsv file containing mappings between \"community\" scaffold names and NCBI scaffold names. If not specified, .tsv file will be inferred from directory contents.")

    # Set of arguments to split the input gff based on an attribute field
    parser.add_argument("--split_gff_when_gene_attribute", required = False, default = None, type = str,
        help = "Condition to split the input gff on in this format <attribute_field>=<something> for example gene_status=other. For multiple filters use a semicolon separated list.")

    # Datachecks
    parser.add_argument("--no_gff_fasta_scaffold_match_dc", action="store_false", dest="dc_fasta_gff_scaffold",
                        help="Do not perform the DataCheck that compares the scaffolds between the gff and fasta. By default, it is implemented."
                             "The offending scaffold names (if any) are stored in in_GFF_not_in_FASTA_scaffolds.txt and/or in_FASTA_not_in_GFF_scaffolds.txt")
    cds_but_no_exons_dc_group = parser.add_mutually_exclusive_group()
    cds_but_no_exons_dc_group.add_argument("--no_cds_with_exons_dc", action="store_false", dest="dc_cds_but_no_exons",
                        help="Do not perform the DataCheck that ensures that all genes with CDS have exons. By default, it is implemented. Offending "
                             "transcripts (if any) are stored in cds_not_exons_transcripts.txt")
    cds_but_no_exons_dc_group.add_argument("--no_fix_cds_with_exons_dc", action="store_false", dest="dc_cds_but_no_exons_fix",
                        help="Do not perform the fix for offending transcripts of the DataCheck that ensures that all genes with CDS have exons. "
                             "By default, the fix is implemented and the offending transcripts are marked as 'non-translating' in the final gff. "
                             "Offending transcripts (if any) are stored in cds_not_exons_transcripts.txt")

    cds_withing_exon_boundaries_group = parser.add_mutually_exclusive_group()
    cds_withing_exon_boundaries_group.add_argument("--no_cds_within_exons_dc", action="store_false", dest="dc_cds_within_exons",
                        help="Do not perform the DataCheck that ensures that all CDS are within exon limits. By default, it is implemented. Offending "
                             "transcripts (if any) are stored in cds_not_within_exons_transcripts.txt")
    cds_withing_exon_boundaries_group.add_argument("--no_fix_cds_within_exons_dc", action="store_false", dest="dc_cds_within_exons_fix",
                        help="Do not perform the fix for offending transcripts of the DataCheck that ensures that all CDS are within exon limits. "
                             "By default, the fix is implemented and the offending transcripts are marked as 'non-translating' in the final gff. "
                             "Offending transcripts (if any) are stored in cds_not_within_exons_transcripts.txt")

    transcripts_with_cds_group = parser.add_mutually_exclusive_group()
    transcripts_with_cds_group.add_argument("--no_coding_transcripts_with_cds_dc", action="store_false",
                                                   dest="dc_coding_transcripts_with_cds",
                                                   help="Do not perform the DataCheck that ensures that all coding transcripts have CDS features. By default, it is implemented. Offending "
                                                        "transcripts (if any) are stored in without_cds_transcripts.txt")
    transcripts_with_cds_group.add_argument("--no_fix_coding_transcripts_with_cds_dc", action="store_false",
                                                   dest="dc_coding_transcripts_with_cds_fix",
                                                   help="Do not perform the fix for offending transcripts of the DataCheck that ensures that ensures that all coding transcripts have CDS features. "
                                                        "By default, the fix is implemented and the offending transcripts are marked as 'non-translating' in the final gff. "
                                                        "Offending transcripts (if any) are stored in without_cds_transcripts.txt")

    args = parser.parse_args()
    return args


# Infers an assembly .fasta file from working directory
# Exits if none is found.
def infer_fasta_file():
    fastas = glob.glob("./*.fa")
    if len(fastas) < 1:
        exit_with_error("No sequence fasta was found")
    else:
        fasta = fastas[0]
    return fasta

# Infers a scaffold synonym file from working directory
# Exits if none or more than 1 found.
def infer_synonyms_file():
    synonyms_files = glob.glob("./*seq_region_synonyms.tsv")
    if len(synonyms_files) != 1:
        exit_with_error("Expected 1 but got " + str(len(synonyms_files)) + " files ending with seq_regions_synonyms.tsv in the directory")
    else:
        synonyms_file = synonyms_files[0]
    return synonyms_file


def main():
    # Grab command-line arguments.
    args = get_args()

    input_gff = args.input_gff
    output_gff = args.output_gff

    # If no .fasta file has been provided, infer one from directory.
    if args.fasta is None:
        fasta_path = infer_fasta_file()
    else:
        fasta_path = args.fasta
    print_info("FASTA file: "+fasta_path)

    fasta = FASTA(fasta_path)

    gff_df = parse_gff(input_gff)

    # TODO: Find a way to not hardcode WormBase_imported
    if args.source_to_WB is True:
        print_info("Switching gff source to : " + wormbase_source)
        gff_df = rename_sources(gff_df, new_source = wormbase_source)

    if args.overwrite_gene_names is True:
        print_info("Inferring and overwriting gene names (Name=) from the ID field (ID=).")
        gff_df = infer_and_overwrite_name_attribute_from_id(gff_df)
    elif args.infer_gene_names is True:
        print_info("Inferring gene names (Name=) from the ID attribute field (ID=).")
        gff_df = infer_name_attribute_from_id(gff_df)

    if args.gene_prefix:
        print_info("Adding the "+args.gene_prefix+" to gene Name and ID attribute fields.")
        gff_df = add_prefix_to_id(gff_df, prefix = args.gene_prefix)
        gff_df = add_prefix_to_name(gff_df, prefix = args.gene_prefix)

    if args.scaffold_rename:
        # If no synonyms.tsv file has been provided, infer one from directory.
        if args.synonyms_file is None:
            synonyms_file = infer_synonyms_file()
        else:
            synonyms_file = args.synonyms_file
        print_info("SYNONYMS file: " + synonyms_file)
        print_info("Renaming GFF scaffolds.")
        gff_df = rename_scaffolds(gff_df, synonyms_file)

    if args.prefixes is not None:
        print_info("Removing prefixes ("+", ".join(args.prefixes)+") from all the fields of the gff file")
        gff_df = remove_prefixes_from_column_values(gff_df, args.prefixes)

    if args.name_prefixes is not None:
        print_info("Removing prefixes ("+", ".join(args.name_prefixes)+") from the Name= field")
        gff_df = remove_prefixes_from_name_column(gff_df, args.name_prefixes)

    if args.extrapolate_transcripts_for_scaffold is not None:
        print_info("Extrapolating transcripts from genes for scaffold: "+rgs.extrapolate_transcripts_for_scaffold)
        gff_df = extrapolate_scaffold_transcripts_from_genes(gff_df, args.extrapolate_transcripts_for_scaffold)
    elif args.extrapolate_genes_for_scaffold is not None:
        print_info("Extrapolating genes from transcripts for scaffold: " + args.extrapolate_genes_for_scaffold)
        gff_df = extrapolate_scaffold_genes_from_transcripts(gff_df, args.extrapolate_transcripts_for_scaffold)

    if args.extrapolate_exons_for_scaffold is not None:
        print_info("Extrapolating exons from CDSs for scaffold: " + args.extrapolate_exons_for_scaffold)
        gff_df = extrapolate_scaffold_exons_from_cds(gff_df, args.extrapolate_exons_for_scaffold)
    elif args.extrapolate_CDSs_for_scaffold is not None:
        print_info("Extrapolating CDSs from exons for scaffold: " + args.extrapolate_CDSs_for_scaffold)
        gff_df = extrapolate_scaffold_cds_from_exons(gff_df, args.extrapolate_CDSs_for_scaffold)

    print_info("Getting Parent Gene for all features ")
    gff_df = get_parent_gene_for_all(gff_df)

    if args.split_gff_when_gene_attribute:
        print_info("Splitting GFF based on the gene attribute field(s): " + " & ".join(args.split_gff_when_gene_attribute.split(";")))
        gff_df, extracted_gff_df = extract_genes_and_features_with_gene_attribute_value(gff_df, args.split_gff_when_gene_attribute)
        extracted_datacheck = DATACHECK(extracted_gff_df, fasta, args, outprefix="extracted_")
        extracted_datacheck.perform_datachecks()
        extracted_output_gff = "extracted.gff3"
        write_output_gff(extracted_gff_df, extracted_output_gff)

    print_info("Running Datachecks")
    datacheck = DATACHECK(gff_df, fasta, args)
    datacheck.perform_datachecks()

    write_output_gff(gff_df, output_gff)


if __name__ == '__main__':
    main()