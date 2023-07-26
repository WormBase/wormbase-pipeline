# import necessary libraries
from ProductionMysql import *
from ProductionUtils import *
import os
import parsers
import utils
import sqlalchemy as sa
from sqlalchemy.orm import sessionmaker

# setting up environmental variables
PARASITE_STAGING = os.environ['PARASITE_STAGING_MYSQL']
PARASITE_VERSION = os.environ['PARASITE_VERSION']
ENSEMBL_VERSION = os.environ['ENSEMBL_VERSION']
WORMBASE_VERSION = os.environ['WORMBASE_VERSION']
WORMBASE_FTP = os.environ['WORMBASE_FTP'] 
PARASITE_SCRATCH = os.environ["PARASITE_SCRATCH"] 

# databases list that comprises the organisms shared between WB and WBPS
databases_search_key = f"_core_{PARASITE_VERSION}_{ENSEMBL_VERSION}_{WORMBASE_VERSION}"
databases = (staging.core_dbs(databases_search_key))

# call get_file_paths function function to get full list of file paths
file_paths = utils.get_file_paths(WORMBASE_FTP, WORMBASE_VERSION, databases)

# call function from parser.py to process functional annotations files and store a .tsv for each organism in scratch directory 
parsers.process_files(file_paths)

# setting up path to scratch directory containing intermediate files
int_folder = 'gene_descs'
int_scratch_directory = os.path.join(PARASITE_SCRATCH, int_folder)

# variable to hold term description so it can be configured easily
gene_attrib_type = 'description'


# for each file in scratch directory, reconstruct file path, open as a df
for file in os.listdir(int_scratch_directory):
    file_path = os.path.join(int_scratch_directory, file)
    df = pd.read_csv(file_path, sep='\t', names=["stable_id", "description"])
    print(len(df))

    # connect to db associated with each .tsv file
    for database in databases:
        core_db = Core(STAGING_HOST, database)

        # write query statements
        # query gene table and return gene_id associated with each stable_id as a tuple
        GENES_QUERY = "SELECT stable_id, gene_id FROM gene"
        # query attrib_type table and return description attrib_type_id
        ATTRIB_QUERY = "SELECT attrib_type_id FROM attrib_type WHERE name = '{}'".format(gene_attrib_type)
        # query to check if any of the genes already have a description associated with them 
        DESCRIP_QUERY = "SELECT DISTINCT(gene_id) FROM gene JOIN gene_attrib USING (gene_id) JOIN attrib_type USING (attrib_type_id) WHERE attrib_type.name='description'"
        # query to delete values from gene_attrib table where there are already descriptions present. They will be replaced with the new descriptions
        DELETE_QUERY = "DELETE FROM gene_attrib WHERE attrib_type_id = '49'"

        # execute queries
        genes_query_execution = core_db.connect().execute(GENES_QUERY)
        gene_id_rows = genes_query_execution.fetchall()
        attrib_query_execution = core_db.connect().execute(ATTRIB_QUERY)
        attrib_type_id = attrib_query_execution.fetchone()[0]
        descrip_query_execute = core_db.connect().execute(DESCRIP_QUERY)
        descrip_q = descrip_query_execute.fetchall()
        
        # create a df from the tuples that are returned as output to the gene query
        # add attrib_type_id column
        gene_id_df = pd.DataFrame(gene_id_rows, columns=['stable_id', 'gene_id'])
        gene_id_df['attrib_type_id'] = attrib_type_id
    
        # merge .tsv df and df created from querying db on the stable_id column, drop rows with 'none available' description
        merged_df = pd.merge(gene_id_df, df[df['description'] != 'none available'], on='stable_id').drop('stable_id', axis=1).rename(columns={'description': 'value'}) 

        # convert df to dictionary to make easier to insert into db?
        data = merged_df.to_dict(orient='records')

        # database insertion
        core_db_w = Core(STAGING_HOST, database, writable=True)
        attrib_table = 'gene_attrib'
        delete_query_execution = core_db_w.connect().execute(DELETE_QUERY)

        # Create the table object using the MyTableName class
        # metadata object is a container for the db schema info
        # gene_attrib_table object represents the gene_attrib db table 
        # The autoload_with parameter specifies that the table structure will be loaded from the database using the connection defined by core_db_w.engine.
        metadata = sa.MetaData()
        gene_attrib_table = sa.Table(attrib_table, metadata, autoload_with=core_db_w.engine)

        # create a session, an intermediate object for interacting with the database
        Session = sessionmaker(bind=core_db_w.engine)
        session = Session()

        for item in data:
        # insert method used to insert data into the table row by row from the list of dictionaries
            ins = gene_attrib_table.insert().values(gene_id=item['gene_id'], attrib_type_id=item['attrib_type_id'], value=item['value'])
            session.execute(ins)

        session.commit()
        session.close()