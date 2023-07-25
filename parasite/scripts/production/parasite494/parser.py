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

# process each functional annotation file and store each stable_id and corresponding description as an output .tsv in a scratch directory
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