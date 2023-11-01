import os
import pandas as pd
from ProductionMysql import *
from ProductionUtils import *

# setting up environmental variables
PARASITE_STAGING = os.environ['PARASITE_STAGING_MYSQL']
PARASITE_VERSION = os.environ['PARASITE_VERSION']
ENSEMBL_VERSION = os.environ['ENSEMBL_VERSION']
WORMBASE_VERSION = os.environ['WORMBASE_VERSION'] 
PARASITE_SCRATCH = os.environ["PARASITE_SCRATCH"] 

# create file path to ouputGenes file so it can be read in
outputGenes = "/hps/nobackup/flicek/wormbase/parasite/stringdb/outputGenes.txt"

# add in dictionary function 
data_dict = {}

# read in output genes file
with open(outputGenes, 'r') as file:
    # iterate through each line in the file
    for line in file:
        # create variable columns, and assign column 1 and column 2 in the tsv to it
        # should be two tab separated columns
        columns = line.strip().split('\t')
        if len(columns) == 2:
            column1, column2 = columns
            # if column 1 is already in the dictionary, add the column 2 value to it
            # else create a new key in the dictionary and add column 2 as the values
            if column1 in data_dict:
                data_dict[column1].append(column2)
            else:
                data_dict[column1] = [column2]

# Print key value pairs
for key, values in data_dict.items():
    print(f'{key}\t{", ".join(values)}')


# count how many gene matches there are for each genome project
for key, value in data_dict.items():
    #print key and number of associated values 
    print(key, len([item for item in value if item]))

# make a list of dictionary keys
keyList = list(data_dict.keys())
print(len(keyList))
keyList.remove('Genome project')


# Need to find a way to print total number of genes
# problem at the moment is that some of the species are wb (290) and others are wbps (1) -> different endings

# take keys and reconstruct db names from them - call this variable databases so it can be fed through the for loop below 
database_version = f"_core_{PARASITE_VERSION}_{ENSEMBL_VERSION}_1"

keyDBs = []
for item in keyList:
    item_db = item + database_version
    keyDBs.append(item_db)


# lists all dbs in wbps 19
# dont think this is necessary?

#databases = (staging.core_dbs(database_version))
#print(databases)


# for loop to go through each db with shared genes, and return the total number of genes wbps holds for the genome
# can then compare the number of string db genes to the overall number held in parasite.
#for database in keyDBs:
    # connect to db
#    core_db = Core(STAGING_HOST, database)
    # genes count query
#    GENES_COUNT_QUERY = "SELECT COUNT(DISTINCT stable_id) FROM gene;"
    # execute queries
#    genes_count_query_execution = core_db.connect().execute(GENES_COUNT_QUERY)