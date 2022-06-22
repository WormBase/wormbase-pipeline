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
import dependencies
import ProductionUtils
from ProductionUtils import *
from tools import *
from dependencies import *

#Parse user input
parser = OptionParser(usage='usage: %prog [options] arguments')
parser.add_option("-s", "--species", dest="SPECIES",
                  help="Required: Species in the core_db format. Example schistosoma_mansoni_prjea36577_core_17_105_1")
parser.add_option("-e", "--explore", dest="explore", action="store_true", default=False,
                  help="Optional: Do not perform any processing. Only show the available FTP locations for the defined samples. Default: Not enabled.")
parser.add_option("-c", "--copy_externally", dest="cpext", action="store_true", default=False ,
                  help="Optional: Copy the result bam/bigwig files externally (rsync). Default: Enabled.")
parser.add_option("-t", "--select_studies", dest="selected_studies", default=[],
                  help="Optional: Comma separated list of studies/samples to be processed. You can specify studies (e.g. DRP003063) and/or samples (e.g. ERR233394). Default: All available.")
user_input = ui(parser)
species = Species(user_input)

print(species.species)


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