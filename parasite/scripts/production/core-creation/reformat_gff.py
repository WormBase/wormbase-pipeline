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
    parser.add_argument("--keep_scaffold_names", action = "store_true", dest = "keep_scaffolds",
        help = "Leave scaffold names unchanged. By default, they are only renamed if the do not match the fasta file's scaffolds.")
    parser.add_argument("--rename_gff_scaffolds_even_if_community", action = "store_true", dest = "rename_gff_scaffolds_if_community",
                        help = "By default, GFF3 scaffold names are only renamed from INSDC IDs to CommunityIDs if they are not "
                        "matching the FASTA file scaffold names. Use this flag to "
                        "enable the renaming of the scaffolds even if they are already Community IDs")
    parser.add_argument("--keep_pseudogenic_CDSs", action = "store_true", dest = "keep_pseudogenic_CDSs",
        help = "If your GFF file contains pseudogenes with CDSs, this option will keep them. By default, they are removed.")
    # has_parentage is set to True by default.
    parser.add_argument("--no_parentage", action = "store_false", dest = "has_parentage",
        help = "Indicate input .gff file has no parentage information to link features. By default, assumes parentage information is included.")
    parser.add_argument("-g", "--gene_prefix", default = False,
        help = "A string to append to each gene and ID name, in case original names are lacking. No changes made by default.")
    parser.add_argument("-p", "--prefixes", required = False, default = None, nargs = '+',
        help = "List of prefixes to remove from .gff fields (e.g. parts of gene IDs that are constant across all genes). \
            Specify multiple prefixes as a space-separated list.")
    parser.add_argument("--name_prefixes", required=False, default=None, nargs='+',
                        help="List of prefixes to remove from .gff fields (e.g. parts of gene Name that are constant across all genes). \
                             Specify multiple prefixes as a space-separated list.")
    parser.add_argument("--merge_adjacent_exons", dest="merge_adjacent_exons", required=False, action="store_true", 
                       help="Merge adjacent exon features that have the same transcript as a parent")
    name_inference_group = parser.add_mutually_exclusive_group()
    name_inference_group.add_argument("-n", "--infer_gene_names", default = False, 
                                      action = "store_true",
                                      help = "Use feature ID attribute to create Name attribues for \
                                      features that lack a name (directy copy)")
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
    parser.add_argument("--no_fix_cds_with_exons_dc", action="store_false", dest="dc_cds_but_no_exons_fix",
                        help="Do not perform the fix for offending transcripts of the DataCheck that ensures that all genes with CDS have exons. "
                             "By default, the fix is implemented and the offending transcripts are marked as 'non-translating' in the final gff. "
                             "Offending transcripts (if any) are stored in cds_not_exons_transcripts.txt")
    parser.add_argument("--no_fix_cds_within_exons_dc", action="store_false", dest="dc_cds_within_exons_fix",
                        help="Do not perform the fix for offending transcripts of the DataCheck that ensures that all CDS are within exon limits. "
                             "By default, the fix is implemented and the offending transcripts are marked as 'non-translating' in the final gff. "
                             "Offending transcripts (if any) are stored in cds_not_within_exons_transcripts.txt")
    parser.add_argument("--no_fix_coding_transcripts_with_cds_dc", action="store_false",
                                                   dest="dc_coding_transcripts_with_cds_fix",
                                                   help="Do not perform the fix for offending transcripts of the DataCheck that ensures that all coding transcripts have CDS features. "
                                                        "By default, the fix is implemented and the offending transcripts are marked as 'non-translating' in the final gff. "
                                                        "Offending transcripts (if any) are stored in without_cds_transcripts.txt")
    parser.add_argument("--no_fix_genes_without_transcripts_dc", action="store_false",
                                                   dest="dc_genes_without_transcripts_fix",
                                                   help="Do not perform the fix for offending genes of the DataCheck that ensures that all all genes have transcripts. "
                                                        "By default, the fix is implemented and transcripts are extrapolated using the coordinates of the offending genes. "
                                                        "The extrapolated transcript IDs is the gene ID + a suffix that can be set with the "
                                                        "--dc_genes_without_transcripts_fix_suffix flag."
                                                        "Offending genes (if any) are stored in genes_without_transcripts.txt")
    parser.add_argument("--fix_genes_without_transcripts_dc_suffix", dest="dc_genes_without_transcripts_fix_suffix", default = "_mRNA", type = str,
        help = "Suffix for fixing the no_genes_without_transcript_datacheck. See --no_fix_genes_without_transcripts_dc documentation for more information."
               "Default value is: _mRNA")
    parser.add_argument("--no_fix_exons_parentage", action="store_false", dest="dc_exons_parentage_fix",
                        help="Do not remove exons that do not have transcript as their parent. By default, they're removed.")
    parser.add_argument("--no_fix_cds_parentage", action="store_false", dest="dc_cds_parentage_fix",
                        help="Do not remove CDS that do not have transcript as their parent. By default, they're removed.")
    parser.add_argument("--no_fix_coding_transcripts_with_exons_and_cds_dc", action="store_false", dest="dc_coding_transcripts_with_exons_and_cds_fix",
                        help="Do not remove coding transcripts that do not have CDS and exons. By default, they're removed.")

    no_exons_cds_group = parser.add_mutually_exclusive_group()
    no_exons_cds_group.add_argument("--turn_coding_transcripts_without_exons_and_cds_to_nontranslating", dest="turn_coding_transcripts_without_exons_and_cds_to_nontranslating", 
                                    action="store_true", help="Turn coding transcripts without exons and CDS to non-translating. By default, they're removed.")
    no_exons_cds_group.add_argument("--extrapolate_missing_exons_and_CDS", dest="extrapolate_missing_exons_and_CDS", 
                                    action="store_true", help="Extrapolate missing exons and CDS. By default, they're removed.")
    no_exons_cds_group.add_argument("--do_not_fix_coding_transcripts_without_exons_and_cds", dest="remove_coding_transcripts_without_exons_and_cds", 
                                    action="store_false", help="Do not fix transcripts with missing exons and CDS. By default, they're removed.")

    cds_without_exons_group = parser.add_mutually_exclusive_group()
    cds_without_exons_group.add_argument("--turn_coding_transcripts_without_exons_but_with_cds_to_nontranslating", dest="turn_coding_transcripts_without_exons_but_with_cds_to_nontranslating",
                                         action="store_true", help="Turn coding transcripts without exons but with CDS to non-translating. By default, they missing exons are extrapolated.")
    cds_without_exons_group.add_argument("--do_not_fix_coding_transcripts_without_exons_but_with_cds", dest="remove_coding_transcripts_without_exons_but_with_cds",
                                         action="store_false", help="Do not fix transcripts with missing exons and CDS. By default, the missing exons are extrapolated.")
    
    exons_without_cds_group = parser.add_mutually_exclusive_group()
    exons_without_cds_group.add_argument("--extrapolate_missing_CDS_from_exons", dest="extrapolate_missing_CDS_from_exons",
                                         action="store_true", help="Extrapolate missing CDS from exons when CDS are missing. By default, these transcripts are turned into non_translating transcripts.")
    exons_without_cds_group.add_argument("--do_not_fix_coding_transcripts_with_exons_but_without_CDS", dest="turn_coding_transcripts_without_CDS_but_with_exons_to_nontranslating",
                                         action="store_false", help="Do not fix transcripts with exons but without CDS. By default, these transcripts are turned into non_translating transcripts.")    
    
    genes_without_tr_group = parser.add_mutually_exclusive_group()
    genes_without_tr_group.add_argument("--remove_genes_without_transcripts", dest="dc_remove_genes_without_transcripts",
                                        action="store_true", help="Remove genes without transcripts and all their features."\
                                        "By default, their transcripts are extrapolated by their gene coordinates.")
    genes_without_tr_group.add_argument("--turn_genes_without_transcripts_to_pseudogenes", dest="dc_turn_genes_without_transcripts_into_pseudogenes",
                                        action="store_true", help="Turn genes without transcripts into pseudogenes."\
                                        "By default, their transcripts are extrapolated by their gene coordinates.")
    genes_without_tr_group.add_argument("--do_not_fix_genes_without_transcripts", dest="dc_genes_without_transcripts_fix",
                                        action="store_false", help="Do not fix genes without transcripts."\
                                        "By default, their transcripts are extrapolated by their gene coordinates.")
    parser.add_argument("--turn_extrap_transcripts_from_genes_into_nontraslating", dest="dc_genes_without_transcripts_into_nontranslating",
                        action="store_true", help="Turn extrapolated transcripts of transcriptless genes into nontranslating_transcripts.")
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

    # If no synonyms.tsv file has been provided, infer one from directory.
    if args.synonyms_file is None:
        synonyms_file = infer_synonyms_file()
    else:
        synonyms_file = args.synonyms_file
    print_info("SYNONYMS file: " + synonyms_file)

    synonyms_df = parse_synonyms(synonyms_file)

    fasta = FASTA(fasta_path)

    gff_df = parse_gff(input_gff, gff_column_names)

    outprefix = os.path.dirname(input_gff)

    initial_feature_count = count_features(gff_df)

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

    if scaffolds_needs_renaming(gff_df, fasta)!=False and not args.keep_scaffolds:
        print_info("Renaming GFF scaffolds.")
        gff_df = rename_scaffolds(gff_df, synonyms_df, args.rename_gff_scaffolds_if_community)

        scaffolds_still_needs_renaming = scaffolds_needs_renaming(gff_df, fasta)
        # if scaffolds still need renaming something went wrong
        if scaffolds_still_needs_renaming:
            with open(os.path.join(outprefix,'scaffolds_could_not_rename.txt'), mode='wt') as myfile:
                    myfile.write('\n'.join(scaffolds_still_needs_renaming))
                    exit_with_error("Exiting as these scaffolds could not be renamed. "
                                    "These scaffolds have been stored in the scaffolds_could_not_rename.txt file. "
                                    "Please rename them manually and rerun the script.")

    
    if args.prefixes is not None:
        print_info("Removing prefixes ("+", ".join(args.prefixes)+") from all the fields of the gff file")
        gff_df = remove_prefixes_from_column_values(gff_df, args.prefixes)

    if args.name_prefixes is not None:
        print_info("Removing prefixes ("+", ".join(args.name_prefixes)+") from the Name= field")
        gff_df = remove_prefixes_from_name_column(gff_df, args.name_prefixes)

    # if args.extrapolate_transcripts_for_scaffold is not None:
    #     print_info("Extrapolating transcripts from genes for scaffold: "+args.extrapolate_transcripts_for_scaffold)
    #     gff_df = extrapolate_scaffold_transcripts_from_genes(gff_df, args.extrapolate_transcripts_for_scaffold)
    # elif args.extrapolate_genes_for_scaffold is not None:
    #     print_info("Extrapolating genes from transcripts for scaffold: " + args.extrapolate_genes_for_scaffold)
    #     gff_df = extrapolate_scaffold_genes_from_transcripts(gff_df, args.extrapolate_transcripts_for_scaffold)

    # if args.extrapolate_exons_for_scaffold is not None:
    #     print_info("Extrapolating exons from CDSs for scaffold: " + args.extrapolate_exons_for_scaffold)
    #     gff_df = extrapolate_scaffold_exons_from_cds(gff_df, args.extrapolate_exons_for_scaffold)
    # elif args.extrapolate_CDSs_for_scaffold is not None:
    #     print_info("Extrapolating CDSs from exons for scaffold: " + args.extrapolate_CDSs_for_scaffold)
    #     gff_df = extrapolate_scaffold_cds_from_exons(gff_df, args.extrapolate_CDSs_for_scaffold)

    print_info("Getting Parent Gene for all features ")
    gff_df = get_parent_gene_for_all(gff_df)

    if gff_df["type"].isin(pseudogene_types).any():
        if pseudogenes_with_CDS(gff_df) and  not args.keep_pseudogenic_CDSs:
            print_info("Found pseudogenes with CDSs. Removing pseudogenes with CDSs. If you would like to keep them " + 
                       "use the keep_pseudogenic_CDSs flag")
            gff_df = remove_pseudogenic_CDS(gff_df)
        print_info("Switching the type of transcripts of pseudogenes to pseudogenic_transcripts")
        gff_df = make_pseudogenic_transcripts(gff_df)
        print_info("Switching the type of \"pseudogene\" genes to \"gene\"")
        gff_df = pseudogenes_to_genes(gff_df)
        

    if args.split_gff_when_gene_attribute:
        print_info("Splitting GFF based on the gene attribute field(s): " + " & ".join(args.split_gff_when_gene_attribute.split(";")))
        gff_df, extracted_gff_df = extract_genes_and_features_with_gene_attribute_value(gff_df, args.split_gff_when_gene_attribute)
        extracted_gff_df = perform_datachecks(extracted_gff_df, args, outprefix=outprefix+".extracted_")
        extracted_output_gff = "extracted.gff3"
        write_output_gff(extracted_gff_df, extracted_output_gff)

    if args.merge_adjacent_exons:
        if all_exon_starts_less_than_ends(gff_df):
            print_info("Merging adjacent exons")
            gff_df = merge_adjacent_exons(gff_df)
        else:
            print_warning("Could not merge adjacent exons as there are exons with start > end. " + 
                          "Please check the GFF file and try again.")
    print_info("Running Datachecks")
    gff_df = perform_datachecks(gff_df, args, outprefix)

    write_output_gff(gff_df, output_gff)

    final_feature_count = count_features(gff_df)
    print_info(report_counts(initial_feature_count))
    print_info(report_counts(final_feature_count))



if __name__ == '__main__':
    main()