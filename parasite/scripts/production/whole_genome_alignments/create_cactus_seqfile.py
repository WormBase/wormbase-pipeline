import os
import re
import sys
import shutil
import glob
import argparse
from ProductionUtils import *
from ProductionMysql import *
from handlers import *
from LSF import *
from Bio import Phylo
import validators


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-r","--release_run",required = False, default = False, action=argparse.BooleanOptionalAction,
        help = "If you want to run cactus a part of the release process switch this on and you don't need to worry about"
               " the other flags.")
    parser.add_argument("-p","--project_name", required = False, default = "WBPS" + os.environ.get("PARASITE_VERSION"),
        help = "Project name that will be used for directory naming etc. Ignore if you used --release_run.")
    parser.add_argument("-c", "--core_db", type = cactus_argparse_core_db, required = False, action = 'append', default = None, nargs = '+',
        help = "Core database name from the current or previous staging you would like to involve in the alignment."
               "Use this flag multiple times for multiple entries. Ignore if you used --release_run."
               "Example: --core_db haemonchus_contortus_prjeb506_core_17_105_1 --core_db schistosoma_mansoni_prjea36577_core_18_108_1")
    parser.add_argument("-f", "--fasta_file", type = cactus_argparse_fasta_file, required = False, action = 'append', default = None, nargs = '+',
        help = "Fasta file path or URL of FASTA file (.gz is supported) you would like to include in the whole genome cactus alignment. "
               "You should also "
               "provide the genome name you would like to be used in the seqfile. Format: <genome_name>:<fasta_file_path>. "
               "Use this flag multiple times for multiple entries."
               "Example: --fasta_file homo_sapiens_38"+FASTA_ARG_DELIMITER+"/nfs/ftp/ensembl/homo_sapiens_38.softmask.fasta "
               "--fasta_file schistosoma_mansoni_prjea36577"+FASTA_ARG_DELIMITER+"https://mansoni/url/genome.fasta.gz. Ignore if you used --release_run.")
    parser.add_argument("-d", "--parent_directory", required = False, default=os.path.join(os.environ.get("PARASITE_SCRATCH"),"cactus_alignments"),
        help = "Parent directory where the analysis will run and the output files will be created. Ignore if you used --release_run.")
    parser.add_argument("-t", "--tree", required = False, default=os.path.join(os.environ.get("PARASITE_CONF"),'cactus.wbparasite.tre'),
        help = "Tree that will be used to infer tree relationships from. Ignore if you used --release_run.")
    args = parser.parse_args()
    return args

args = get_args()

FASTA_ARG_DELIMITER=";"
PARASITE_FTP_URL_BACKBONE="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases"
PARASITE_FTP_URL_PREVIOUS_RELEASE_SPECIES=os.path.join(PARASITE_FTP_URL_BACKBONE,"WBPS{0}".format(os.environ.get("PREVIOUS_PARASITE_VERSION")),"species")
FASTA_PATH_SUFFIX=".genomic_softmasked.fa.gz"
CHECKSUM_FILE_URL="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS{0}/CHECKSUMS".format(os.environ.get("PREVIOUS_PARASITE_VERSION"))
PARASITE_DATA_DIR=os.environ.get("PARASITE_DATA")
PARASITE_VERSION=os.environ.get("PARASITE_VERSION")
PREVIOUS_PARASITE_VERSION=os.environ.get("PREVIOUS_PARASITE_VERSION")
staging_cores = staging.core_databases
previous_staging_cores = previous_staging.core_databases
all_staging_cores = staging_cores + previous_staging_cores
download_psftp_fasta_bash_script = "{0}/parasite/scripts/production/whole_genome_alignments/download_psftp_fasta.sh".format(os.environ.get("WORM_CODE"))
project_dir=os.path.join(args.parent_directory, args.project_name)
print_info("PROJECT DIRECTORY:" + project_dir)
input_dir=os.path.join(project_dir, "input", "run")
fasta_dir = os.path.join(input_dir, "fasta")
dump_genomes_log_dir=os.path.join(project_dir, "dump_genomes_log")
print_info("INPUT DIRECTORY:" + input_dir)
seqfile=os.path.join(input_dir, args.project_name + ".seqfile")

if not os.path.exists(fasta_dir):
    os.makedirs(fasta_dir)
if not os.path.exists(dump_genomes_log_dir):
    os.makedirs(dump_genomes_log_dir)

# Run on release mode for production
if args.release_run:

    print_info("Running in release mode")
    # Total core dbs
    core_dbs_to_align = staging.release_core_databases

    # Split between core dbs which have been updated (we need to dump their fasta files) and
    # core dbs which have not been updated and we can hence download their fasta files for the
    # ftp of the live release.
    core_dbs_to_align_notupdated, core_dbs_to_align_updated = updated_or_not_updated_core_dbs(core_dbs_to_align, PARASITE_DATA_DIR, staging)
    release_fasta_args = [core_db_to_align_notupdated_to_fasta_args(x, PARASITE_VERSION, PREVIOUS_PARASITE_VERSION,
                                                                    previous_staging, PARASITE_FTP_URL_PREVIOUS_RELEASE_SPECIES,
                                                                    FASTA_ARG_DELIMITER, FASTA_PATH_SUFFIX) for x in core_dbs_to_align_notupdated]

    print_w_indent("{0} updated and {1} not updated core dbs have been identified".format(str(len(core_dbs_to_align_updated)),str(len(core_dbs_to_align_notupdated))))

# Merge with other cores or fasta files that might have been specified through the command line options
dumpable_core_dbs = (["{0}:{1}".format(x,staging.host) for x in core_dbs_to_align_updated] if args.release_run else []) + \
                    (flatten(args.core_db) if args.core_db else [])
fasta_args = (release_fasta_args if args.release_run else []) + \
             (args.fasta_file if args.fasta_file else [])

# Create bash commands for either dumping or downloading fasta files and store them in dictionaries
print_info("Generating bash commands")
fasta_commands_dict = {x.split(FASTA_ARG_DELIMITER)[0]:create_fasta_processing_command(x, fasta_dir, FASTA_ARG_DELIMITER, PARASITE_FTP_URL_BACKBONE, download_psftp_fasta_bash_script) for x in fasta_args}
core_dump_commands_dict = {parasite_core2genome(x.split(":")[0]):create_core_dumping_command(x, fasta_dir) for x in dumpable_core_dbs}

# Unique genomes QC
all_genomes = list(fasta_commands_dict.keys()) + list(core_dump_commands_dict.keys())
if len(all_genomes) != len(set(all_genomes)):
    print_error("Duplicate genomes found:")
    seen = set()
    uniq = [x for x in all_genomes if x not in seen and not seen.add(x)]
    print_w_indent(" ".join(uniq))
    exit_with_error("")

# Parse Tree
print_info("Parsing tree")
tree =  read_newick_tree_file(args.tree)
clades = tree_clades(tree)

# Tree genomes match FASTA genomes QC
# if set(all_genomes) != set(clades):
#     exit_with_error("The production names for the genomes you selected do not match with the ones in the tree file: "+args.tree+".")

# Submitting jobs
print_info("Submitting FASTA download/dump jobs to LSF")
jobs_dict = {}
for genome in fasta_commands_dict:
    job_id = lsf_submit(command=fasta_commands_dict[genome],
                        jobprefix=genome + "_dump",
                        cpu=1, mem="1gb",
                        cwd=dump_genomes_log_dir,
                        queue="short",
                        only_write=False)
    jobs_dict[genome]=job_id

for genome in core_dump_commands_dict:
    job_id = lsf_submit(command=core_dump_commands_dict[genome],
                        jobprefix=genome + "_dump",
                        cpu=1, mem="8gb",
                        cwd=dump_genomes_log_dir,
                        queue="production",
                        only_write=False)
    jobs_dict[genome] = job_id
print_info("Jobs submitted.")

# LSF scheduler monitoring
batch_jobs_check_status_periodically(jobs_dict)

# Write final output
print_info("Writing seqfile {0}".format(seqfile))
write_seqfile(tree, all_genomes, seqfile, fasta_dir)

print_info("Done")