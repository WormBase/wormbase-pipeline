import sys
import os
import csv
import re
import pysam
import sysrsync
import gzip
import glob
from pathlib import Path
import dependencies
import diotools
from diotools import *
from tools import *
from dependencies import *

#Parse user input
SPECIES = sys.argv[1]
if len(sys.argv) != 2:
    print(dtnow() + ": ERROR - Usage: python create_bigwig_tracks.py SPECIES_NAME\n"
                    "Example: python create_bigwig_tracks.py strongyloides_stercoralis")
    raise ValueError

# SPECIES="strongyloides_stercoralis"
species = Species(SPECIES)

#Check that dependent directories exist
check_dependent_dirs()

#Create and print out needed directories if not exist
species.create_dirs()
species.report_input()

#Iterrate over studies:
for spec_study in species.studies():
    sst = Study(species.species, spec_study)
    sst.create_study_dirs()
    #print(sst.merge_and_cap_bams_command())
    if sst.study_id == "ERP001556":
        print(sst.study_id)
        sst.merge_and_cap_bams_command_submit()
    else:
        pass
    #break
    #sst.create_study_dirs()
    #print_info("Copying study's required reference fasta file from FTP.")
    #sst.copy_ref_fasta()
    # for dek in sst.design_keys():
    #     ssa = Sample(sst.species, spec_study, dek)
    #     print(ssa.bigwig)
        #ssa.create_sample_dirs()
        #print_info("Submitting copy jobs. Copying BAM/CRAM and BigWig files from FTP.")
        #copy_job = ssa.copy_command_submit()
        #print_info("Submitting individual BAM files processing jobs.")
        #process_job = ssa.bam_processing_command_submit(copy_job)
print_info("Done. Goodbye.")




