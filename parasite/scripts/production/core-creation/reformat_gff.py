#!/hps/software/users/wormbase/parasite/shared/.pyenv/versions/p395/bin/python
from gff_utils import *
from Bio import SeqIO
import re
import sys
import glob
from argparse import ArgumentParser


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

	parser.add_argument("-F", "--fabricate_transcripts_for_scaffold", required = False, default = None,
		help = "Use gene information to create missing transcript features for a specified scaffold (provide scaffold name).")


	name_inference_group = parser.add_mutually_exclusive_group()
	name_inference_group.add_argument("-n", "--infer_gene_names", default = False, action = "store_true",
		help = "Use feature ID attribute to create Name attribues for features that lack a name (directy copy)")
	name_inference_group.add_argument("-N", "--overwrite_gene_names", default = False, action = "store_true",
		help = "Use feature ID attribute to create Name attribues for all features (directy copy). Overwrites original names.")


	parser.add_argument("-f", "--fasta", required = False, default = None, type = str,
		help = ".fasta file containing assembly scaffolds. If not specified, .fasta file will be inferred from directory contents.")
	parser.add_argument("-s", "--synonyms_file", required = False, default = None, type = str,
		help = ".tsv file containing mappings between \"community\" scaffold names and NCBI scaffold names. If not specified, .tsv file will be inferred from directory contents.")

	args = parser.parse_args()
	return args 


# Infers an assembly .fasta file from working directory
# Exits if none is found.
def infer_fasta_file():
	fastas = glob.glob("./*.fa")
	if len(fastas) < 1:
		sys.exit("No sequence fasta was found")
	else:
		fasta = fastas[0]
	return fasta

# Infers a scaffold synonym file from working directory
# Exits if none or more than 1 found.
def infer_synonyms_file():
    synonyms_files = glob.glob("./*seq_region_synonyms.tsv")
    if len(synonyms_files) != 1:
       	sys.exit("Expected 1 but got " + str(len(synonyms_files)) + " files ending with seq_regions_synonyms.tsv in the directory")
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
		fasta = infer_fasta_file()
	else:
		fasta = args.fasta

	# If no synonyms.tsv file has been provided, infer one from directory.
	if args.synonyms_file is None:
		synonyms_file = infer_synonyms_file()
	else:
		synonyms_file = input.synonyms_file

	gff_df = parse_gff(input_gff)

	# TODO: Find a way to not hardcode WormBase_imported
	if args.source_to_WB is True:
		gff_df = rename_sources(gff_df, new_source = "WormBase_imported")

	if args.overwrite_gene_names is True:
		gff_df = infer_and_overwrite_name_attribute_from_id(gff_df)
	elif args.infer_gene_names is True:
		gff_df = infer_name_attribute_from_id(gff_df)

	if args.gene_prefix:
		gff_df = add_prefix_to_id(gff_df, prefix = args.gene_prefix)
		gff_df = add_prefix_to_name(gff_df, prefix = args.gene_prefix)

	if args.scaffold_rename:
		gff_df = rename_scaffolds(gff_df, synonyms_file)

	if args.prefixes is not None:
		print(args.prefixes)
		gff_df = remove_prefixes_from_column_values(gff_df, args.prefixes)

	if args.fabricate_transcripts_for_scaffold is not None:
		gff_df = extrapolate_missing_transcripts_for_scaffold(gff_df, args.fabricate_transcripts_for_scaffold)

	write_output_gff(gff_df, output_gff)
	

if __name__ == '__main__':
	main()