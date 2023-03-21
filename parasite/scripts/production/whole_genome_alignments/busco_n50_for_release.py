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

release = sys.argv[1]

FASTA_ARG_DELIMITER=";"
PARASITE_FTP_URL_BACKBONE="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases"
PARASITE_FTP_URL_PREVIOUS_RELEASE_SPECIES=os.path.join(PARASITE_FTP_URL_BACKBONE,"WBPS{0}".format(os.environ.get("PREVIOUS_PARASITE_VERSION")),"species")
FASTA_GENOMIC_SUFFIX=".genomic.fa.gz"
FASTA_PROTEIN_SUFFIX=".protein.fa.gz"
CHECKSUM_FILE_URL="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/WBPS{0}/CHECKSUMS".format(os.environ.get("PREVIOUS_PARASITE_VERSION"))
PARASITE_DATA_DIR=os.environ.get("PARASITE_DATA")
PARASITE_VERSION=os.environ.get("PARASITE_VERSION")
PARASITE_SOFTWARE=os.environ.get("PARASITE_SOFTWARE")
PREVIOUS_PARASITE_VERSION=os.environ.get("PREVIOUS_PARASITE_VERSION")
BUSCO3_CONTAINER=os.environ.get("BUSCO3_CONTAINER")
AUGUSTUS_PATH=os.path.join(PARASITE_SOFTWARE, "Augustus", "config")
BUSCO_LINEAGES=os.environ.get("BUSCO_LINEAGES")
NEMATODA_DB="nematoda_odb9"
METAZOA_DB="metazoa_odb9"
BUSCO_EXEC_COMMAND_WITH_AUGUSTUS = f"SINGULARITYENV_AUGUSTUS_CONFIG_PATH={AUGUSTUS_PATH} singularity exec {BUSCO3_CONTAINER} run_BUSCO.py"
BUSCO_EXEC_COMMAND = f"singularity exec {BUSCO3_CONTAINER} run_BUSCO.py"
staging_cores = staging.core_databases
previous_staging_cores = previous_staging.core_databases
all_staging_cores = staging_cores + previous_staging_cores
download_psftp_fasta_bash_script = "{0}/parasite/scripts/production/whole_genome_alignments/download_psftp_fasta.sh".format(os.environ.get("WORM_CODE"))
project_dir="/hps/nobackup/flicek/wormbase/parasite/paper_stats/"

fasta_dir=os.path.join(project_dir,"WBPS"+str(release),"fasta")
busco_dir=os.path.join(project_dir,"WBPS"+str(release),"busco")
logs_dir=os.path.join(project_dir,"WBPS"+str(release),"logs")
if not os.path.exists(project_dir): os.makedirs(project_dir)
if not os.path.exists(fasta_dir): os.makedirs(fasta_dir)
if not os.path.exists(fasta_dir): os.makedirs(busco_dir)
if not os.path.exists(logs_dir): os.makedirs(logs_dir)

def parameter_for_augustus(phylum):
    if phylum=="Nematoda":
        return "caenorhabditis"
    elif phylum=="Platyhelminthes":
        return "schistosoma"
    else:
        exit_with_error("Don't know what to do with this phylum: "+phylum+". Exiting.")

def busco_library(phylum, busco_lineages_dir, nematoda_odb_dir, platyhelminthes_odb_dir):
    if phylum=="Nematoda":
        return os.path.join(busco_lineages_dir,nematoda_odb_dir)
    elif phylum=="Platyhelminthes":
        return os.path.join(busco_lineages_dir,platyhelminthes_odb_dir)
    else:
        exit_with_error("Don't know what to do with this phylum: "+phylum+". Exiting.")


def run_busco_command(busco_path, fasta_file, output_dir, outname, run_mode, database, num_threads=4, lineage_dataset=None, augustus_species='auto', extra_options="", tee_log_path=""):
    """
    returns a bash command to execute the run of BUSCO on a fasta file.
    :param fasta_file: The path to the input fasta file.
    :param output_dir: The path to the output directory.
    :param database: The name of the BUSCO database to use.
    :param num_threads: The number of threads to use for the BUSCO run.
    :param lineage_dataset: The path to a directory containing lineage-specific information for the BUSCO run.
    :param augustus_species: The species for which to configure Augustus gene finding. Default is 'auto'.
    :param extra_options: Any additional options to include in the final command.
    :return: A string containing the bash command.
    """
    str_num_threads=str(num_threads)
    command = f"BUSCO_TMP={output_dir}; mkdir -p {output_dir}; cd {output_dir}\n{busco_path} -i {fasta_file} -o {outname} -l {database} -m {run_mode} -c {str_num_threads}"
    if lineage_dataset is not None:
        command += f" -f {lineage_dataset}" + (" -r" if run_mode!="proteins" else "")
    if augustus_species != 'auto':
        command += f" -sp {augustus_species}"
    command += f" {extra_options}"
    if tee_log_path!="":
        command += f" | tee {tee_log_path}"
    command += ";"
    return command

def busco_cmd_wrap(genome, phylum, analysis_suffix, fasta_suffix, run_mode, busco_dir=busco_dir, BUSCO_EXEC_COMMAND_WITH_AUGUSTUS=BUSCO_EXEC_COMMAND_WITH_AUGUSTUS, BUSCO_LINEAGES=BUSCO_LINEAGES,
                                NEMATODA_DB=NEMATODA_DB, METAZOA_DB=METAZOA_DB):
    genome_busco_command = run_busco_command(busco_path=BUSCO_EXEC_COMMAND_WITH_AUGUSTUS,
                                             fasta_file=os.path.join(fasta_dir,genome+fasta_suffix),
                                             output_dir=os.path.join(busco_dir, genome),
                                             run_mode=run_mode,
                                             outname=genome + "_" + analysis_suffix,
                                             database=busco_library(phylum, BUSCO_LINEAGES,
                                                                    NEMATODA_DB, METAZOA_DB),
                                             num_threads=8, lineage_dataset="",
                                             augustus_species=parameter_for_augustus(phylum),
                                             tee_log_path=os.path.join(busco_dir, genome,
                                                      analysis_suffix,
                                                                       "run_busco.out"))
    return genome_busco_command

fasta_commands_dict={}
core_dump_commands_dict={}

if str(release) == str(PARASITE_VERSION):
    core_dbs = staging.release_core_databases
    core_dbs_notupdated, core_dbs_updated = updated_or_not_updated_core_dbs(core_dbs, PARASITE_DATA_DIR, staging)
    release_fasta_args = [[core_db_to_align_notupdated_to_fasta_args(x, PARASITE_VERSION, PREVIOUS_PARASITE_VERSION,
                                                                    previous_staging, PARASITE_FTP_URL_PREVIOUS_RELEASE_SPECIES,
                                                                    FASTA_ARG_DELIMITER, FASTA_GENOMIC_SUFFIX),
                           core_db_to_align_notupdated_to_fasta_args(x, PARASITE_VERSION, PREVIOUS_PARASITE_VERSION,
                                                                     previous_staging,
                                                                     PARASITE_FTP_URL_PREVIOUS_RELEASE_SPECIES,
                                                                     FASTA_ARG_DELIMITER, FASTA_PROTEIN_SUFFIX)] for x in core_dbs_notupdated]

    dumpable_core_dbs=["{0}:{1}".format(x,staging.host) for x in core_dbs_updated]

    fasta_commands_dict = {x[0].split(FASTA_ARG_DELIMITER)[0]:[create_fasta_processing_command(x[0], fasta_dir, FASTA_ARG_DELIMITER, PARASITE_FTP_URL_BACKBONE, download_psftp_fasta_bash_script, file_suffix="genomic"),
                                                               create_fasta_processing_command(x[1], fasta_dir, FASTA_ARG_DELIMITER, PARASITE_FTP_URL_BACKBONE, download_psftp_fasta_bash_script, file_suffix="protein")] for x in release_fasta_args}
    core_dump_commands_dict = {parasite_core2genome(x.split(":")[0]):[create_core_dumping_command(x, fasta_dir, file_suffix="genomic", genome_dump=True, protein_dump=False, canonical_only=False),
                                                                      create_core_dumping_command(x, fasta_dir, file_suffix="protein", genome_dump=False, protein_dump=True, canonical_only=True)
                                                                      ] for x in dumpable_core_dbs}

    for genome in fasta_commands_dict:
        fasta_commands_dict[genome].append(busco_cmd_wrap(genome, staging.species_phylum(genome), analysis_suffix="assembly",
                                                                       fasta_suffix=".genomic.fa", run_mode="genome"))
        fasta_commands_dict[genome].append(busco_cmd_wrap(genome, staging.species_phylum(genome), analysis_suffix="annotation",
                                        fasta_suffix=".protein.fa", run_mode="proteins"))
        fasta_commands_dict[genome].append(calculate_n50_command(fasta=os.path.join(fasta_dir, genome + ".genomic.fa"),
                                                                 outfile=os.path.join(busco_dir, genome, "n50.txt")))

    for genome in core_dump_commands_dict:
        core_dump_commands_dict[genome].append(
            busco_cmd_wrap(genome, staging.species_phylum(genome), analysis_suffix="assembly",
                                        fasta_suffix=".genomic.fa", run_mode="genome"))
        core_dump_commands_dict[genome].append(
            busco_cmd_wrap(genome, staging.species_phylum(genome), analysis_suffix="annotation",
                                        fasta_suffix=".protein.fa", run_mode="proteins"))
        core_dump_commands_dict[genome].append(calculate_n50_command(fasta=os.path.join(fasta_dir, genome + ".genomic.fa"),
                                                                 outfile=os.path.join(busco_dir, genome, "n50.txt")))


if int(release) < int(PARASITE_VERSION):
    release_genomes = get_previous_release_genomes_list(release)

    release_fasta_args = [[x+FASTA_ARG_DELIMITER+get_genome_ftp_file_url_from_ftp_release(x, version=release, filetype="gfasta"),
                           x+FASTA_ARG_DELIMITER+get_genome_ftp_file_url_from_ftp_release(x, version=release, filetype="pfasta")] for x in release_genomes]

    fasta_commands_dict = {x[0].split(FASTA_ARG_DELIMITER)[0]: [
        create_fasta_processing_command(x[0], fasta_dir, FASTA_ARG_DELIMITER, PARASITE_FTP_URL_BACKBONE,
                                        download_psftp_fasta_bash_script, file_suffix="genomic"),
        create_fasta_processing_command(x[1], fasta_dir, FASTA_ARG_DELIMITER, PARASITE_FTP_URL_BACKBONE,
                                        download_psftp_fasta_bash_script, file_suffix="protein")] for x in
                           release_fasta_args}

    for genome in fasta_commands_dict:
        species = "_".join(genome.split("_")[0:2])
        fasta_commands_dict[genome].append(
            busco_cmd_wrap(genome, staging.species_phylum(species), analysis_suffix="assembly",
                           fasta_suffix=".genomic.fa", run_mode="genome"))
        fasta_commands_dict[genome].append(
            busco_cmd_wrap(genome, staging.species_phylum(species), analysis_suffix="annotation",
                           fasta_suffix=".protein.fa", run_mode="proteins"))
        fasta_commands_dict[genome].append(calculate_n50_command(fasta=os.path.join(fasta_dir,genome+".genomic.fa"),outfile=os.path.join(busco_dir,genome,"n50.txt")))


jobs_dict = {}
for genome in fasta_commands_dict:
    job_id = lsf_submit(command="\n\n".join(fasta_commands_dict[genome]),
                        jobprefix=genome,
                        cpu=8, mem="32gb",
                        cwd=logs_dir,
                        queue="production",
                        only_write=False)
    jobs_dict[job_id]=genome

for genome in core_dump_commands_dict:
    job_id = lsf_submit(command="\n\n".join(core_dump_commands_dict[genome]),
                        jobprefix=genome,
                        cpu=8, mem="32gb",
                        cwd=logs_dir,
                        queue="production",
                        only_write=False)
    jobs_dict[genome]=job_id

all_jobs = list(jobs_dict.keys())
while len(all_jobs)>0:
    for job_id in jobs_dict:
        ref = jobs_dict[job_id]
        if lsf_job_status(job_id)!="RUN" and lsf_job_status(job_id)!="PEND":
            if job_id in all_jobs: all_jobs.remove(job_id)
            # print_info(job_id +" is DONE")
            if lsf_job_status(job_id) == "EXIT":
                print_info(ref+" - "+job_id + " has FAILED.")
    print_info("Number of jobs still running: " + str(len(all_jobs)))
    time.sleep(60)