# import necessary libraries
from ProductionMysql import *
from ProductionUtils import *
import os
import gzip
import csv

# setting up environmental variable for database server
PARASITE_STAGING = os.environ['PARASITE_STAGING_MYSQL']

# databases list that comprises the organisms shared between WB and WBPS
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

# list comprehension to feed each element in the databases list into the get_file_path() function
# outputs a list containing the file path for each element in the initial databases list 
annotation_output = [get_file_path(x) for x in databases] 

# store the output of the get_file_path() function for each database in a separate variable.
b_malayi, c_brenneri, c_briggsae, c_elegans, c_japonica, c_remanei, o_volvulus, p_pacificus, s_ratti, t_muris = annotation_output

# construct file paths
malayi_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', b_malayi)
brenneri_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', c_brenneri)
briggsae_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', c_briggsae)
elegans_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', c_elegans)
japonica_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', c_japonica)
remanei_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', c_remanei)
volvulus_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', o_volvulus)
pacificus_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', p_pacificus)
ratti_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', s_ratti)
muris_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', t_muris)

# list of file paths
# to feed into for loop below
file_paths = [malayi_path, brenneri_path, briggsae_path, elegans_path, japonica_path, remanei_path, volvulus_path, pacificus_path, ratti_path, muris_path]

for file_path in file_paths:
    # create a unique output file name based on the input file name
    output_file_name = os.path.basename(file_path) + '.tsv'
    
    # reset the data list to an empty list
    data = []
    
    with gzip.open(file_path, 'rt') as file:
        stable_id = None
        description = None

        for line in file:
            line = line.strip()
            if line.startswith('WBGene'):
                stable_id = line.split()[0]
            elif line.startswith('Automated description:'):
                description = line.replace('Automated description:', '').strip()
                for subline in file:
                    subline = subline.strip()
                    if subline.startswith('Gene class description:'):
                        data.append([stable_id, description])
                        break
                    else:
                        description += ' ' + subline
                stable_id = None
                description = None

    # print the output file name and the data for debugging purposes
    print(output_file_name)
    print(data)

    # write .tsv file with extracted data
    with open(output_file_name, 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(['stable_id', 'description'])
        for row in data:
            writer.writerow(row)
