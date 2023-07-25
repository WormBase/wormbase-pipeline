# import necessary libraries
from ProductionMysql import *
from ProductionUtils import *

# function to create part of file path to wb functional annotations files from string of the WB database name
#  s = string of db name
def get_path_from_db_name(s, WORMBASE_VERSION):
    organism_genus = s.split("_")[0][0]
    organism_species = s.split("_")[1]
    bioproject_id = s.split("_")[2].upper()
    return f"{organism_genus}_{organism_species}/{bioproject_id}/annotation/{organism_genus}_{organism_species}.{bioproject_id}.WS{WORMBASE_VERSION}.functional_descriptions.txt.gz"

# function that uses get_path_from_db_name() from utils.py to reconstruct whole file path for each wb functional annotations file
def get_file_paths(WORMBASE_FTP, WORMBASE_VERSION, databases):
    annotation_output = [get_path_from_db_name(x, WORMBASE_VERSION) for x in databases]
    file_paths = [os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', path) for path in annotation_output]
    return file_paths


# nested function get_file_path() will take db string and return organism genus _ organism species _ bioproject id for each database
#def get_file_paths(WORMBASE_FTP, WORMBASE_VERSION, databases):
#    def get_file_path(s):
#        organism_genus = s.split("_")[0][0]
#        organism_species = s.split("_")[1]
#        bioproject_id = s.split("_")[2].upper()
#        return f"{organism_genus}_{organism_species}/{bioproject_id}/annotation/{organism_genus}_{organism_species}.{bioproject_id}.WS{WORMBASE_VERSION}.functional_descriptions.txt.gz"

#    annotation_output = [get_file_path(x) for x in databases]

#    file_paths = [os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', path) for path in annotation_output]
#    return file_paths






