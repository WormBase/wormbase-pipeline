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
attrib_id = '49'



# each output .tsv in int_scratch_directory is paired with its corresponding database from the databases list
# returned as a tuple
# does this based on them being in alphabetical order 
file_database_pairs = zip(os.listdir(int_scratch_directory), databases)

records_added_count = {}

for file, database in file_database_pairs:
    file_path = os.path.join(int_scratch_directory, file)
    df = pd.read_csv(file_path, sep='\t', names=["stable_id", "description"])

    # connect to db
    core_db = Core(STAGING_HOST, database)
    #print(file_path)
    #print(core_db)

    # write query statements
    # query gene table and return gene_id associated with each stable_id as a tuple
    GENES_QUERY = "SELECT stable_id, gene_id FROM gene"
    # query attrib_type table and return description attrib_type_id
    ATTRIB_QUERY = "SELECT attrib_type_id FROM attrib_type WHERE name = '{}'".format(gene_attrib_type)
    # query to check if any of the genes already have a description associated with them 
    DESCRIP_QUERY = "SELECT DISTINCT(gene_id) FROM gene JOIN gene_attrib USING (gene_id) JOIN attrib_type USING (attrib_type_id) WHERE attrib_type.name='description'"
    # query to delete values from gene_attrib table where there are already descriptions present. They will be replaced with the new descriptions
    DELETE_QUERY = "DELETE FROM gene_attrib WHERE attrib_type_id = '{}'".format(attrib_id)
    # data check query
    DATA_CHECK_QUERY = "SELECT COUNT(*) FROM gene_attrib WHERE attrib_type_id = '{}'".format(attrib_id)

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

    # Store the count of dropped 'none available' descriptions
    num_dropped_rows = len(df[df['description'] == 'none available'])

    # Print the number of remaining rows after dropping 'none available' values
    print(f"Database: {database}")
    print(f"Total number of automated descriptions to add: {len(merged_df)}")
    print(f"Number of 'none available' values dropped: {num_dropped_rows}")

    # convert df to dictionary to make easier to insert into db?
    data = merged_df.to_dict(orient='records')

    # database insertion
    core_db_w = Core(STAGING_HOST, database, writable=True)
    attrib_table = 'gene_attrib'
    # execute query to delete existing descriptions from the gene_attrib table before inserting the new descriptions
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

    # Execute the SQL query to count the records added for the current database
    num_records_added = core_db_w.connect().execute(DATA_CHECK_QUERY).fetchone()[0]

    # Check if the number of records added matches the length of the merged_df
    matching_length = num_records_added == len(merged_df)

    # Print the count of records added for the current database and whether it matches the length of merged_df
    print(f"Total number of descriptions added: {num_records_added}. This matches initial number of descriptions to add: {len(merged_df)}")


    # Raise SystemError if the records added do not match the length of the merged_df
    if not matching_length:
        raise SystemError(f"ERROR: Number of descriptions to be added ({num_records_added}) does not match the total number expected: ({len(merged_df)})")

    session.commit()
    session.close()