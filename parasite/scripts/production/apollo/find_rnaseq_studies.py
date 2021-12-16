from optparse import OptionParser
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
parser = OptionParser(usage='usage: %prog [-s SPECIES options]')
parser.add_option("-s", "--species", dest="SPECIES",
                  help="Required: Species in the 'genus_species' format. Example strongyloides_stercoralis")
(options, args) = parser.parse_args()

print(options)
if options.SPECIES:
    SPECIES = options.SPECIES
else:
    parser.error('species not given')

species = Species(SPECIES)

explore_print(species)
exit()