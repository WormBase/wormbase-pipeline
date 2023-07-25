# import necessary libraries
from ProductionMysql import *
from ProductionUtils import *
import os
import gzip
import csv
import parser
import utils
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

# call function from utils.py to get the list of file paths
file_paths = utils.get_file_paths(WORMBASE_FTP, WORMBASE_VERSION, databases)

# call function from parser.py to process functional annotations files and store a .tsv for each organism in scratch directory 
parser.process_files(file_paths)

int_folder = 'gene_descs'
int_scratch_directory = os.path.join(PARASITE_SCRATCH, int_folder)
gene_attrib_type = 'description'

for file in os.listdir(int_scratch_directory):
    file_path = os.path.join(int_scratch_directory, file)
    df = pd.read_csv(file_path, sep='\t', names=["stable_id", "description"])
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