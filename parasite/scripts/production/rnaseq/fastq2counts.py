#!/hps/software/users/wormbase/parasite/shared/.pyenv/versions/p395/bin/python

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
from tools import *
from dependencies import *

# Parse user input
parser = OptionParser(usage='usage: %prog [options] arguments')
parser.add_option("-s", "--core_db", dest="COREDB",
                  help="Required: Species in the core_db format. Example schistosoma_mansoni_prjea36577_core_17_105_1")
parser.add_option("-e", "--explore", dest="explore", action="store_true", default=False,
                  help="Optional: Do not perform any processing. Only show the available FTP locations for the defined samples. Default: Not enabled.")
parser.add_option("-g", "--genome_reference", dest="genome_ref_dir",
                  help="Optional: Directory with STAR genome reference for the required species. If not provided, we'll try to access it at " + reference_genomes_dir)
parser.add_option("-m", "--gtf", dest="genome_gtf",
                  help="Optional: Directory with gene models gtf for the genome reference used. If not provided, we'll try to access it at " + reference_genomes_dir)
parser.add_option("-c", "--copy_externally", dest="cpext", action="store_true", default=False,
                  help="Optional: Copy the result bam/bigwig files externally (rsync). Default: Enabled.")
parser.add_option("-p", "--deploy_to_apollo", dest="apdep", action="store_true", default=False,
                  help="Optional: Deploy result files to apollo. Used with -a option. Default: Not enabled.")
parser.add_option("-a", "--apollo_instance_name", dest="apname", help="Required when -p is used. Name of the apollo instance that will be deployed.")
parser.add_option("-t", "--select_studies_samples", dest="selects", default=[],
                  help="Optional: Comma separated list of studies/samples to be processed. Needs to be in this format: SECONDARY_STUDY_ID1:RUN_ID1;RUN_ID2,SECONDARY_STUDY_ID2,SECONDARY_STUDY_ID3:RUN_ID1;RUN_ID2,"
                       "Example: ERP001675:ERR278825;ERR278827,ERP113120,ERP002073:ERR667078")

check_dependent_dirs()

user_input = ui(parser)
user_input.print_user_input()

species = Species(user_input)
species.create_dirs()
species_all_studies_dict = all_studies_dict_for_species(species.get_rna_seq_studies_by_taxon_api_url)
species_selected_studies_dict = selected_studies_dict_for_species(user_input, species_all_studies_dict)
selected_species = Species(user_input, species_selected_studies_dict)

for study_id in selected_species.studies:
    study = Study(user_input, species_selected_studies_dict, study_id)
    study.create_dirs()
    print_info("Processing Study " + study.study_id)
    print_w_indent("With samples:\n" + "\n".join(study.sample_ids))
    for sample_id in study.sample_ids:
        print_info("Processing Sample: " + sample_id)
        sample = Sample(user_input, species_selected_studies_dict, study_id, sample_id)
        sample.create_dirs()
        if sample.no_of_fastqs != 2:
            print_warning("Sample "+sample_id+" doesn't have 2 fastq files. Skipping...")
            continue
    #     # sample_processing_command = sample.download_fastqs_command() + \
    #     #                             sample.trim_fastqs_command(delete_previous=False) + \
    #     #                             sample.star_alignment_command(delete_previous=False) + \
    #     #                             sample.sortrefname_bam_command() + \
    #     #                             sample.cap_bam_command(delete_previous=False) + \
    #     #                             sample.bam2bigwig_commnd(delete_previous=False) + \
    #     #                             study.move_all_finals_to_final_dir_command() + \
    #     #                             study.deploy_to_embassy_command()
    #
        sample_processing_command = sample.download_fastqs_command() + \
                                    sample.trim_fastqs_command(delete_previous=True) + \
                                    sample.star_alignment_command(delete_previous=True) + \
                                    sample.index_aligned_bam_command() + \
                                    sample.sortrefname_bam_command() + \
                                    sample.cap_bam_command(delete_previous=True) + \
                                    sample.sort_capped_bam_command(delete_previous=True) + \
                                    sample.index_final_bam_command() + \
                                    sample.bam2bigwig_command(delete_previous=False) + \
                                    sample.move_final_bam_to_apollo_dir_command() + \
                                    sample.move_final_bigwig_to_apollo_dir_command()

        sample_processing_command = sample.download_fastqs_command() + \
                                    sample.trim_fastqs_command(delete_previous=True) + \
                                    sample.star_alignment_command(delete_previous=True) + \
                                    sample.index_aligned_bam_command() + \
                                    sample.bam2bigwig_command(delete_previous=True) + \
                                    sample.move_final_bigwig_to_apollo_dir_command()

        sample_processing_submit = lsf_submit(sample_processing_command,
                                              jobprefix = sample.sample_id + "_newbigwigs",
                                              to_wait_id="",
                                              cpu=6, mem="20gb",
                                              cwd=sample.sample_dir,
                                              queue="production",
                                              only_write=False)

        print_info("Submitted job: " +  sample_processing_submit + " for sample: " + sample.sample_id)
        # study_processing_command = study.move_all_finals_to_final_dir_command() + \
        #                            study.deploy_to_embassy_command()
        #
        # study_processing_submit = lsf_submit(study_processing_command,
        #                                       jobprefix = study.study_id,
        #                                       to_wait_id="",
        #                                       cpu=1, mem="2gb",
        #                                       cwd=study.log_dir,
        #                                       queue="production",
        #                                       only_write=False)

# print(sample.download_fastqs_command())
# print(sample.trim_fastqs_command())
# print(sample.star_alignment_command())
# print(sample.sortrefname_bam_command())
# print(sample.cap_bam_command())
# print(sample.bam2bigwig_commnd())
# print(study.move_all_finals_to_final_dir_command())
# print(study.deploy_to_embassy_command())




# species = Species(user_input)
#
# print(species.species)
# print(species.taxon_id)
# print([x for x in species.all_studies_dict()["ERP001675"] if x["secondary_sample_accession"]=="ERS323732"])
# if user_input.explore:
#     explore_print(user_input)
#     exit()
#
# # Check that dependent directories exist
# check_dependent_dirs()
#
# # Create and print out needed directories if not exist
# species.create_dirs()
# species.report_input()
#
# #Find which studies/samples are needed based on the users select_studies input
# studies_to_use = studies_needed(user_input.selected_studies, species)
# samples_to_use = samples_needed(user_input.selected_studies, species)
