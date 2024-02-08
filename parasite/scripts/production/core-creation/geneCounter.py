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
databases_search_key = f"_core_{PARASITE_VERSION}_{ENSEMBL_VERSION}_1"

# to search in wb species databases use the below line of code instead of the above 
#databases_search_key = f"_core_{PARASITE_VERSION}_{ENSEMBL_VERSION}_{WORMBASE_VERSION}"

databases = (staging.core_dbs(databases_search_key))
print(databases)

total_stable_ids = {}

for database in databases:
    # connect to db
    core_db = Core(STAGING_HOST, database)

    # create sql query to count the number of distinct genes in the genes table 
    GENES_QUERY = "SELECT COUNT(DISTINCT stable_id) FROM gene"
    DESCRIPTIONS_QUERY = "SELECT COUNT(DISTINCT stable_id) FROM gene WHERE description is NOT NULL"

    # execute the query
    genes_query_execution = core_db.connect().execute(GENES_QUERY)

    # fetch the result
    result = genes_query_execution.fetchone()

    # store the result
    total_stable_ids[database] = result[0]

# Output the total count of stable IDs for each database
for database, total_count in total_stable_ids.items():
    print(f"Total distinct stable IDs in {database}: {total_count}")

# Define a variable to store the total count of stable IDs across all databases
total_count_all_databases = 0

# Iterate over the total counts for each database and sum them up
for total_count in total_stable_ids.values():
    total_count_all_databases += total_count

# Output the final count
print(f"Total distinct stable IDs across all databases: {total_count_all_databases}")
