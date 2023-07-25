# import necessary libraries
from ProductionMysql import *
from ProductionUtils import *

# main function get_file_paths() 
# nested function get_file_path() will take db string and return organism genus _ organism species _ bioproject id for each database
def get_file_paths(WORMBASE_FTP, WORMBASE_VERSION, databases):
    def get_file_path(s):
        organism_genus = s.split("_")[0][0]
        organism_species = s.split("_")[1]
        bioproject_id = s.split("_")[2].upper()
        return f"{organism_genus}_{organism_species}/{bioproject_id}/annotation/{organism_genus}_{organism_species}.{bioproject_id}.WS{WORMBASE_VERSION}.functional_descriptions.txt.gz"

    annotation_output = [get_file_path(x) for x in databases]

    file_paths = [os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', path) for path in annotation_output]
    return file_paths