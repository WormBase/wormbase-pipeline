# import necessary libraries
from ProductionMysql import *
from ProductionUtils import *
import gzip

# setting up database server
PARASITE_STAGING = os.environ['PARASITE_STAGING_MYSQL']

# get list stored in the variable databases that comprises the organisms shared between WB and WBPS
databases = (staging.core_dbs("core_18_108_285"))

# creating a file path to each of the shared databases
WORMBASE_FTP = os.environ['WORMBASE_FTP'] # '/nfs/ftp/public/databases/wormbase/'
WORMBASE_VERSION = os.environ['WORMBASE_VERSION'] # '285'

# function that will take db string and return organism genus _ organism species _ bioproject id
def get_file_path(s):
    organism_genus = s.split("_")[0][0]
    organism_species = s.split("_")[1]
    bioproject_id = s.split("_")[2].upper()
    return f"{organism_genus}_{organism_species}/{bioproject_id}/annotation/{organism_genus}_{organism_species}.{bioproject_id}.WS{WORMBASE_VERSION}.functional_descriptions.txt.gz"

#s = "brugia_malayi_prjna10729_core_18_108_285"
#print(get_file_path(s))

# list comprehension to feed each element in the databases list into the get_file_path() function
# outputs a list containing the file path for each element in the initial databases list 
annotation_output = [get_file_path(x) for x in databases] 


# store the output of the get_file_path() function for each database in a separate variable.
b_malayi, c_brenneri, c_briggsae, c_elegans, c_japonica, c_remanei, o_volvulus, p_pacificus, s_ratti, t_muris = annotation_output

#b_malayi = annotation_output[0]
#c_brenneri = annotation_output[1]
#c_briggsae = annotation_output[2]
#c_elegans = annotation_output[3]
#c_japonica = annotation_output[4]
#c_remanei = annotation_output[5]
#o_volvulus = annotation_output[6]
#p_pacificus = annotation_output[7]
#s_ratti = annotation_output[8]
#t_muris = annotation_output[9]

# need to add it to the file path you have created before
elegans_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', c_elegans)

# reading the content from the file path
with gzip.open(elegans_path, 'rb') as f:
    file_content = f.read()

# the variable file_content now contains the contents of the c_elegans.PRJNA13758.WS285.functional_descriptions.txt.gz file