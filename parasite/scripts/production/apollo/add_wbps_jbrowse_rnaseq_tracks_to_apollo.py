#!/home/ubuntu/python/v3.9.5/bin/python3 python3.9.5

# This script automates the updating process of a JSON file used in an Apollo instance, 
# which is a genome annotation editor. It fetches additional track information from a 
# WormBase ParaSite JBrowse instance for the specified species and adds it to the existing 
# track list in the Apollo JSON file.

# Please make sure that the genome version in the WBPS release version you specify here 
# matches the release version used in your apollo instance.

# WBPS JBROWSE FTP URL
jbrowse_ftp_url = "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/web_data/jbrowse/releases/release-{0}/{1}/data/"

import os
import json
import requests
import argparse

# Parse user input
parser = argparse.ArgumentParser(usage='usage: %(prog)s [options] arguments')
parser.add_argument("-c", "--species", dest="species", required=True,
                    help="Required: Species in the PS format. Example schistosoma_mansoni_prjea36577")
parser.add_argument("-a", "--apollo_instance_directory", dest="apdir",
                    help="Full path to the apollo instance directory that contains the trackList.json file.")
parser.add_argument("-p", "--wbps_release", dest="release",
                    help="WormBase ParaSite Release Number (Example: 18)")

# Parse the command-line arguments
args = parser.parse_args()

# Access the values of the options
species = args.species
apdir = args.apdir
release = args.release

# Functions
def load_json_from_url(url):
    """
    Loads JSON data from a URL.

    Args:
        url (str): The URL to load JSON data from.

    Returns:
        dict: The JSON data loaded as a Python dictionary, or None if the request fails.
    """
    response = requests.get(url)
    if response.status_code == 200:
        data = response.json()
        return data
    else:
        print(f"Failed to load JSON from URL: {url}")
        return None

def write_to_apollo_trackList_json(data, output_file, write=True):
    """
    Writes dictionary data to an output JSON file with proper indentation and backups the file if it exists.

    Args:
        data (dict): The dictionary data to be written as JSON.
        output_file (str): The file path to write the JSON data to.

    Returns:
        None
    """
    if write:
        # Backup the output file if it exists
        if os.path.isfile(output_file):
            backup_file = output_file + '.bak'
            os.rename(output_file, backup_file)
            print(f"Backed up existing file to: {backup_file}")

        # Write the dictionary to the output JSON file
        with open(output_file, 'w') as f:
            json.dump(data, f, indent=4)
            print(f"Data written to: {output_file}")
    else:
        print(data)

jbrowse_url = jbrowse_ftp_url.format(release,species)
jbrowse_json = f"{jbrowse_url}/trackList.json"
apollo_path = apdir
apollo_json = f"{apdir}/trackList.json"

# Load JBrowse trackList.json from the provided URL
jb_tldata = load_json_from_url(jbrowse_json)

# Open and load the existing Apollo trackList.json
ap_tlopen = open(apollo_json)
ap_tldata = json.load(ap_tlopen)

# Extract JBrowse bigwig tracks and add them to Apollo trackList.json
jb_bigwig_tracks = [x for x in jb_tldata["tracks"] if 'type' in x.keys() and x['type']=='JBrowse/View/Track/Wiggle/XYPlot' and x['urlTemplate'].endswith(".bw")]
ap_tldata["tracks"] += jb_bigwig_tracks

# Write the updated Apollo trackList.json back to the file
write_to_apollo_trackList_json(ap_tldata, apollo_json, write=True)



