# import necessary libraries
from ProductionMysql import *
from ProductionUtils import *
import os
import gzip
import csv

# setting up environmental variables
PARASITE_STAGING = os.environ['PARASITE_STAGING_MYSQL']
PARASITE_VERSION = os.environ['PARASITE_VERSION']
ENSEMBL_VERSION = os.environ['ENSEMBL_VERSION']
WORMBASE_VERSION = os.environ['WORMBASE_VERSION']
WORMBASE_FTP = os.environ['WORMBASE_FTP'] # '/nfs/ftp/public/databases/wormbase/'
PARASITE_SCRATCH = os.environ["PARASITE_SCRATCH"] 

# databases list that comprises the organisms shared between WB and WBPS
databases_search_key = f"_core_{PARASITE_VERSION}_{ENSEMBL_VERSION}_{WORMBASE_VERSION}"
databases = (staging.core_dbs(databases_search_key))

# creating a file path to each of the shared databases
# main function get_file_paths() 
# nested function get_file_path() will take db string and return organism genus _ organism species _ bioproject id for each database
# move this to ProductionMysql.py
def get_file_paths(WORMBASE_FTP, WORMBASE_VERSION, databases):
    def get_file_path(s):
        organism_genus = s.split("_")[0][0]
        organism_species = s.split("_")[1]
        bioproject_id = s.split("_")[2].upper()
        return f"{organism_genus}_{organism_species}/{bioproject_id}/annotation/{organism_genus}_{organism_species}.{bioproject_id}.WS{WORMBASE_VERSION}.functional_descriptions.txt.gz"

    annotation_output = [get_file_path(x) for x in databases]

    file_paths = [os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', path) for path in annotation_output]
    return file_paths

# call function to get the list of file paths
file_paths = get_file_paths(WORMBASE_FTP, WORMBASE_VERSION, databases)


# creating a file path to each of the shared databases
# function that will take db string and return organism genus _ organism species _ bioproject id 
# move this to ProductionMysql.py
#def get_file_path(s):
#    organism_genus = s.split("_")[0][0]
#    organism_species = s.split("_")[1]
#    bioproject_id = s.split("_")[2].upper()
#    return f"{organism_genus}_{organism_species}/{bioproject_id}/annotation/{organism_genus}_{organism_species}.{bioproject_id}.WS{WORMBASE_VERSION}.functional_descriptions.txt.gz"

# list comprehension to feed each element in the databases list into the get_file_path() function
# outputs a list containing the file path for each element in the initial databases list 
#annotation_output = [get_file_path(x) for x in databases] 

# store the output of the get_file_path() function for each database in a separate variable.
#b_malayi, c_brenneri, c_briggsae, c_elegans, c_japonica, c_remanei, o_volvulus, p_pacificus, s_ratti, t_muris = annotation_output

# construct file paths
#malayi_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', b_malayi)
#brenneri_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', c_brenneri)
#briggsae_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', c_briggsae)
#elegans_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', c_elegans)
#japonica_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', c_japonica)
#remanei_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', c_remanei)
#volvulus_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', o_volvulus)
#pacificus_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', p_pacificus)
#ratti_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', s_ratti)
#muris_path = os.path.join(WORMBASE_FTP, 'releases', 'WS' + WORMBASE_VERSION, 'species', t_muris)

# list of file paths
# to feed into for loop below
#file_paths = [malayi_path, brenneri_path, briggsae_path, elegans_path, japonica_path, remanei_path, volvulus_path, pacificus_path, ratti_path, muris_path]

def process_files(file_paths):
    # variables to store hard coded text
    WB_gene_prefix = 'WBGene'
    AD_prefix = 'Automated description:'
    GD_prefix = 'Gene class description:'

    # create folder within scratch server to store intermediate files
    int_folder = 'gene_descs'
    # path to new folder
    int_folder_path = os.path.join(PARASITE_SCRATCH, int_folder)

    # Use try-except block to handle the case if the folder already exists
    try:
    # Create the folder
        os.mkdir(int_folder_path)
        print(f"Folder '{int_folder}' created successfully at: {int_folder_path}")
    except FileExistsError:
    # If the folder already exists, remove it and then create a new one
        try:
            shutil.rmtree(int_folder_path)  # Remove the existing folder and its contents
            os.mkdir(int_folder_path)       # Create a new folder with the same name
            print(f"Folder '{int_folder}' replaced successfully at: {int_folder_path}")
        except Exception as e:
            print(f"Error occurred while replacing the folder: {e}")
    except Exception as e:
        print(f"Error occurred while creating the folder: {e}")

    for file_path in file_paths:
        # create output file name based on the input functional annotations file name
        output_file_name = os.path.basename(file_path) + '.tsv'

        # reset the data list to an empty list
        data = []

        with gzip.open(file_path, 'rt') as file:
            stable_id = None
            description = None

            for line in file:
                line = line.strip()
                if line.startswith(WB_gene_prefix):
                    stable_id = line.split()[0]
                elif line.startswith(AD_prefix):
                    description = line.replace(AD_prefix, '').strip()
                    for subline in file:
                        subline = subline.strip()
                        if subline.startswith(GD_prefix):
                            data.append([stable_id, description])
                            break
                        else:
                            description += ' ' + subline
                    stable_id = None
                    description = None
            print(output_file_name)
        
        # Write .tsv file with extracted data inside the "gene_descs" folder
        output_file_path = os.path.join(int_folder_path, output_file_name)
        with open(output_file_path, 'w', newline='') as file:
            writer = csv.writer(file, delimiter='\t')
            writer.writerow(['stable_id', 'description'])
            for row in data:
                writer.writerow(row)

# call function
process_files(file_paths)



#for file_path in file_paths:
    # create a unique output file name based on the input file name
    #output_file_name = os.path.basename(file_path) + '.tsv'
    
    # reset the data list to an empty list
    #data = []
    # variables to store hard coded text
    #WB_gene_prefix = 'WBGene'
    #AD_prefix = 'Automated description:'
    #GD_prefix = 'Gene class description:'

    #with gzip.open(file_path, 'rt') as file:
        #stable_id = None
        #description = None

        #for line in file:
            #line = line.strip()
            #if line.startswith(WB_gene_prefix):
                #stable_id = line.split()[0]
            #elif line.startswith(AD_prefix):
                #description = line.replace(AD_prefix, '').strip()
                #for subline in file:
                    #subline = subline.strip()
                    #if subline.startswith(GD_prefix):
                        #data.append([stable_id, description])
                        #break
                    #else:
                        #description += ' ' + subline
                #stable_id = None
                #description = None

    # print the output file name and the data for debugging purposes
    #print(output_file_name)
    #print(data)

# create new directory within scratch to store the intermediate .tsv files


# create a new directory within PARASITE_SCRATCH = os.environ["PARASITE_SCRATCH"] (/hps/nobackup/flicek/wormbase/parasite)
# this will store the intermediate .tsv files generated by this script

    # write .tsv file with extracted data
    #with open(output_file_name, 'w', newline='') as file:
        #writer = csv.writer(file, delimiter='\t')
        #writer.writerow(['stable_id', 'description'])
        #for row in data:
            #writer.writerow(row)
