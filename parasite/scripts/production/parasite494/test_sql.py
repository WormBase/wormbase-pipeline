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
# database staging server
STAGING_HOST = os.environ['PARASITE_STAGING_MYSQL']
PARASITE_VERSION = os.environ['PARASITE_VERSION']
ENSEMBL_VERSION = os.environ['ENSEMBL_VERSION']
WORMBASE_VERSION = os.environ['WORMBASE_VERSION']

# list of databases
databases = (staging.core_dbs("core_18_108_285"))

# import .tsv file containing stable_id and automated description as a pandas dataframe?
df=pd.read_csv('data.tsv', sep='\t',names=["stable_id","description",])
stable_id = df['stable_id'].tolist()
# loop through WB/WBPS databases one by one using the stable_id to query the gene table and return the gene_id
# set database 
mansoni_db = Core(staging.host, "schistosoma_mansoni_prjea36577_core_18_108_1")
execution = mansoni_db.connect().execute("SELECT gene_id FROM gene WHERE stable_id IN (stable_id)");
print([x for x in execution])