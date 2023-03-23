from ProductionMysql import *
from ProductionUtils import *
import os
import subprocess

repeat_file_suffix = "repeat-families.fa"
PARASITE_VERSION = os.environ["PARASITE_VERSION"]
PREVIOUS_PARASITE_VERSION = os.environ["PREVIOUS_PARASITE_VERSION"]
ENSEMBL_VERSION = os.environ["ENSEMBL_VERSION"]
PARASITE_FTP_DIR = os.environ["PARASITE_FTP"]
PREVIOUS_PARASITE_FTP_DIR = os.path.join(PARASITE_FTP_DIR,"releases",f"WBPS{PREVIOUS_PARASITE_VERSION}","species")

dump_dir = f"/hps/nobackup/flicek/wormbase/parasite/dumps/WBPS{PARASITE_VERSION}/repeats/"
new_repeats_dir = f"/hps/nobackup/flicek/wormbase/parasite/repeats/WBPS{PARASITE_VERSION}/work"
release_cores = staging.release_core_databases
release_genomes = staging.release_genomes
new_repeats_genomes = os.listdir(new_repeats_dir)

if not os.path.exists(dump_dir):
    os.makedirs(dump_dir)

for nrc in new_repeats_genomes:
    if nrc not in release_genomes:
        print_error(f"The {nrc} path from {new_repeats_dir} doesn't correspond to a core db for this release. Exiting.")

for rg in release_genomes:

    print_info(f"Processing {rg}")

    if rg in new_repeats_genomes:

        print_w_indent(f"{rg} is a new genome in this release. Trying to find the families.fa.gz file")
        repeat_dir = os.path.join(new_repeats_dir,rg)
        r_file = os.path.join(repeat_dir, f"{rg}-families.fa")

    elif rg not in new_repeats_genomes:
        print_w_indent(f"{rg} is not a new genome. Copying from previous release")
        if rg=="panagrolaimus_ju765_prjeb32708":
            rg_ftp_fix = Core(previous_staging.host, "propanagrolaimus_ju765_prjeb32708").ftp_filename_n_filename_without_version()
        elif rg=="panagrolaimus_sp1159_prjeb32708":
            rg_ftp_fix = Core(previous_staging.host,"panagrolaimus_ps1159_prjeb32708").ftp_filename_n_filename_without_version()
        else:
            rg_ftp_fix = Core(staging.host,rg).ftp_filename_n_filename_without_version()
        # rg_ftp_fix = Core(staging.host, rg).ftp_filename_n_filename_without_version()
        r_file = os.path.join(PREVIOUS_PARASITE_FTP_DIR,f"{rg_ftp_fix}WBPS{PREVIOUS_PARASITE_VERSION}.{repeat_file_suffix}.gz")

    if not os.path.exists(r_file):
        exit_with_error(f"{r_file} does not exist")

    target_file = os.path.join(dump_dir,f"{rg}-families.fa") + (".gz" if r_file.endswith(".gz") else "")

    print_w_indent(f"Copying {r_file} to {target_file}")
    command = f"cp {r_file} {target_file};"
    if not r_file.endswith(".gz"):
        command += f" gzip {target_file};"
    command += "\n"
    run_code = subprocess.call(command, shell=True)
    if run_code!=0:
        exit_with_error(f"Error when executing {command}")
    else:
        print_w_indent("Done\n")





