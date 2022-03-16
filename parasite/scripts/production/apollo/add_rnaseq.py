import sys
import os
import csv
import re
import pysam
import sysrsync
import gzip
import glob
from optparse import OptionParser
from pathlib import Path
import dependencies
import diotools
from diotools import *
from tools import *
from dependencies import *

#Parse user input
parser = OptionParser(usage='usage: %prog [options] arguments')
parser.add_option("-s", "--species", dest="SPECIES",
                  help="Required: Species in the 'genus_species' format. Example strongyloides_stercoralis")
parser.add_option("-e", "--explore", dest="explore", action="store_true", default=False,
                  help="Optional: Do not perform any processing. Only show the available FTP locations for the defined samples. Default: Not enabled.")
parser.add_option("-c", "--copy_externally", dest="cpext", action="store_true", default=True,
                  help="Optional: Copy the result bam/bigwig files externally (rsync). Default: Enabled.")
parser.add_option("-d", "--external_destination", dest="extdest", default=ssh_host,
                  help="Optional: External destination for copying evidence track files in the form of 'digri@codon-login:/path/to/files'."
                       "Deafult: Sanger NGS Server apollo directory (see dependencies.py for path).")
parser.add_option("-t", "--select_studies", dest="selected_studies", default=[],
                  help="Optional: Comma separated list of studies/samples to be processed. You can specify studies (e.g. DRP003063) and/or samples (e.g. ERR233394). Default: All available.")
user_input = ui(parser)
species = Species(user_input.species)


if user_input.explore:
    explore_print(species)
    exit()

# Check that dependent directories exist
check_dependent_dirs()

# Create and print out needed directories if not exist
species.create_dirs()
species.report_input()

#Find which studies/samples are needed based on the users select_studies input
studies_to_use = studies_needed(user_input.selected_studies, species)
samples_to_use = samples_needed(user_input.selected_studies, species)

print("studies to use:")
print(studies_to_use)
print("samples to use:")
print(samples_to_use)

# Iterrate over studies:
for spec_study in species.studies():
    sst = Study(species.species, spec_study)
    sst.create_study_dirs()
    print_info("Copying study's required reference fasta file from FTP.")
    sst.copy_ref_fasta()
    for dek in sst.design_keys():
        ssa = Sample(sst.species, spec_study, dek)
        print(ssa.bigwig)
        ssa.create_sample_dirs()
        print_info("Submitting copy jobs. Copying BAM/CRAM and BigWig files from FTP.")
        copy_job = ssa.copy_command_submit()
        print_info("Submitting individual BAM files processing jobs.")
        process_job = ssa.bam_processing_command_submit(copy_job)
    sst.merge_and_cap_bams_command_submit()
print_info("Done. Goodbye.")




