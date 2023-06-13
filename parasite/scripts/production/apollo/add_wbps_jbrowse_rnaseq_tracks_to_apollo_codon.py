import sys
import os
import json
import tempfile
import argparse
from dependencies import *
from ProductionUtils import *
from tools import *

# Parse user input
parser = argparse.ArgumentParser(usage='usage: %(prog)s [options] arguments')
parser.add_argument("-c", "--species", dest="species", required=True,
                    help="Required: Species in the PS format. Example schistosoma_mansoni_prjea36577")
parser.add_argument("-a", "--apollo_instance_name", dest="apname",
                    help="Name of the apollo instance that will be deployed.")

# Parse the command-line arguments
args = parser.parse_args()

# Access the values of the options
species = args.species
apname = args.apname

core_jbrowse_path = jbrowse_ftp_dir + f"/{species}/data"
core_jbrowse_json = f"{core_jbrowse_path}/trackList.json"
apollo_remote_path = f"{aws_apollo_address}:{aws_apollo_filepath}/data/{apname}"
apollo_remote_path_tl = f"{aws_apollo_address}:{aws_apollo_filepath}/data/{apname}/trackList.json"

jb_tlopen = open(core_jbrowse_json)
jb_tldata = json.load(jb_tlopen)

ap_tlopen = open_remote_file(apollo_remote_path_tl)
ap_tldata = json.load(ap_tlopen)

jb_bigwig_tracks = [x for x in jb_tldata["tracks"] if 'type' in x.keys() and x['type']=='JBrowse/View/Track/Wiggle/XYPlot' and x['urlTemplate'].endswith(".bw")]
ap_tldata["tracks"] += jb_bigwig_tracks

write_to_remote_apollo_trackList_json(ap_tldata, species, apollo_remote_path_tl, write=False)



