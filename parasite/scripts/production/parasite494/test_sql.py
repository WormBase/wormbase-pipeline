# import necessary libraries
from ProductionMysql import *
from ProductionUtils import *
import sqlalchemy
from sqlalchemy import insert
import os
import re
from ensembl.database import DBConnection, UnitTestDB
import csv

# setting up environmental variables 
# database staging server etc
STAGING_HOST = os.environ['PARASITE_STAGING_MYSQL']
PARASITE_VERSION = os.environ['PARASITE_VERSION']
ENSEMBL_VERSION = os.environ['ENSEMBL_VERSION']
WORMBASE_VERSION = os.environ['WORMBASE_VERSION']

# list of databases containing WB and WBPS organisms 
databases = (staging.core_dbs("core_18_108_285"))


# STARTING JUST WITH c_japonica
# for each organism: HOW WOULD I IMPLEMENT A FOR LOOP THAT LOADS EACH ORGANISMS .tsv AS A SEPARATE DF?
# import .tsv file containing stable_id and automated description as a pandas dataframe
df=pd.read_csv('c_japonica.PRJNA12591.WS285.functional_descriptions.txt.gz.tsv', sep='\t',names=["stable_id","description",])
# get list of stable_ids that can be used to query database gene table and return gene_id 
stable_id = df['stable_id'].tolist()

# loop through WB/WBPS databases one by one using the stable_id to query the gene table and return stable_id and their associated gene_id as a tuple
# set database 
japonica_db = Core(staging.host, "caenorhabditis_japonica_prjna12591_core_18_108_285")
execution = japonica_db.connect().execute("SELECT stable_id, gene_id FROM gene WHERE stable_id IN (stable_id)")
# return list of the tuples
gene_id = execution.fetchall()
attrib_type = japonica_db.connect().execute("SELECT attrib_type_id from attrib_type WHERE name = 'description' ")
at = attrib_type.fetchall()

# add descriptions to list of tuples 
# first create a dataframe of the stable_ids and the gene_ids
df2 = pd.DataFrame(gene_id, columns=['stable_id', 'gene_id'])
df2['attrib_type_id'] = str(at[0]).split(',')[0].strip('()')
# merge this dataframe with the initial dataframe containing the stable_id and the description
# this creates a 
merged_df = pd.merge(df2, df, on='stable_id', how='inner')
merged_df = merged_df[~merged_df['description'].str.contains('none available')]

print(merged_df)
# convert dataframe into a dictionary
#data = merged_df.to_dict(orient='records')


insert_query = "INSERT INTO your_table (gene_id, attrib_type_id, value) VALUES (%s, %s, %s)"

japonica_db.connect().execute(insert_query, data) # data will be the rows from the merged_df?
