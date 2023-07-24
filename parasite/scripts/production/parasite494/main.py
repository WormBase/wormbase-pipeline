# import necessary libraries
from ProductionMysql import *
from ProductionUtils import *
import os
import gzip
import csv
from sqlalchemy import insert, Table, MetaData
from sqlalchemy.orm import sessionmaker

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

# call function 
process_files(file_paths)

int_folder = 'gene_descs'
int_scratch_directory = os.path.join(PARASITE_SCRATCH, int_folder)
gene_attrib_type = 'description'


for file in os.listdir(int_scratch_directory):
    file_path = os.path.join(int_scratch_directory, file)
    df = pd.read_csv(file, sep='\t', names=["stable_id", "description"])
    #print(df)
    stable_ids = df['stable_id'].tolist()
    # connect to db associated with each .tsv file
    for database in databases:
        core_db = Core(STAGING_HOST, database)

        # write queries
        # query gene table and return gene_id associated with each stable_id as a tuple
        GENES_QUERY = "SELECT stable_id, gene_id FROM gene"
        # query attrib_type table and return description attrib_type_id
        ATTRIB_QUERY = "SELECT attrib_type_id FROM attrib_type WHERE name = '{}'".format(gene_attrib_type)
        
        execution = core_db.connect().execute(GENES_QUERY)
        gene_id_rows = execution.fetchall()
        # create a df from the tuples
        gene_id_df = pd.DataFrame(gene_id_rows, columns=['stable_id', 'gene_id'])
        attrib_execution = core_db.connect().execute(ATTRIB_QUERY)
        attrib_type_id = attrib_execution.fetchone()[0]
        # add attrib type id as a column to df
        gene_id_df['attrib_type_id'] = attrib_type_id
        # merge .tsv df and df created from querying db on the stable_id column 
        merged_df = pd.merge(gene_id_df, df, on='stable_id', how='inner')
        merged_df = merged_df[~merged_df['description'].str.contains('none available')]
        merged_df = merged_df.drop('stable_id', axis=1)
        merged_df = merged_df.rename(columns={'description': 'value'})
        # convert df to dictionary to make easier to insert into db?
        data = merged_df.to_dict(orient='records')
        print(data)

        conn = core_db.engine.connect()