import os
import pandas as pd
from ProductionMysql import *
from ProductionUtils import *
from sqlalchemy import insert

# Set up environmental variables
STAGING_HOST = os.environ['PARASITE_STAGING_MYSQL']

# List of databases containing WB and WBPS organisms
databases = staging.core_dbs("core_18_108_285")

# function that will take db string and return organism genus _ organism species _ bioproject id
def get_file_path(s):
    organism_genus = s.split("_")[0][0]
    organism_species = s.split("_")[1]
    bioproject_id = s.split("_")[2].upper()
    return f"{organism_genus}_{organism_species}.{bioproject_id}.WS{WORMBASE_VERSION}.functional_descriptions.txt.gz.tsv"

# list comprehension to feed each element in the databases list into the get_file_path() function
# outputs a list containing the file name for each element in the initial databases list 
organism_files = [get_file_path(x) for x in databases]

for file in organism_files:
    # Import .tsv file containing stable_id and automated description as a pandas dataframe
    df = pd.read_csv(file, sep='\t', names=["stable_id", "description"])

    # Get list of stable_ids that can be used to query the gene table and return gene_id
    stable_ids = df['stable_id'].tolist()

    # Loop through WB/WBPS databases one by one using the stable_id to query the gene table and return stable_id and their associated gene_id as a tuple
    for database in databases:
        # Set database
        core_db = Core(STAGING_HOST, database)
        execution = core_db.connect().execute("SELECT stable_id, gene_id FROM gene WHERE stable_id IN :stable_ids", {"stable_ids": stable_ids})

        # Return list of tuples
        gene_id_rows = execution.fetchall()

        # Add descriptions to list of tuples
        # First create a dataframe of the stable_ids and the gene_ids
        gene_id_df = pd.DataFrame(gene_id_rows, columns=['stable_id', 'gene_id'])
        
        # Query attrib_type table to get attrib_type_id for description
        # add attrib_type as a new column to df
        attrib_type_query = core_db.connect().execute("SELECT attrib_type_id FROM attrib_type WHERE name = 'description'")
        attrib_type_id = attrib_type_query.fetchone()[0]
        gene_id_df['attrib_type_id'] = attrib_type_id

        # Merge this dataframe with the initial dataframe containing the stable_id and the description
        merged_df = pd.merge(gene_id_df, df, on='stable_id', how='inner')

        # Filter out rows with 'none available' description
        merged_df = merged_df[~merged_df['description'].str.contains('none available')]

        # Drop the 'stable_id' column
        # rename the description column 'value' so it matches the gene_attrin table
        merged_df = merged_df.drop('stable_id', axis=1)
        merged_df = merged_df.rename(columns={'description': 'value'})
        
        
        # Convert merged_df to a list of tuples
        data = merged_df.to_records(index=False).tolist()
        insert_statment = "INSERT INTO gene_attrib (gene_id, attrib_type_id, value) VALUES (%s, %s, %s)"
        core_db.connect().execute(insert_statment, data)
"""         # convert merged_df to dictionary to insert as new rows in the gene_attrib table
        data = merged_df.to_dict(orient='records')

        # Create the 'gene_attrib' table object using the 'meta' attribute from the 'core_db' object
        insert_statment = "INSERT INTO gene_attrib (gene_id, attrib_type_id, value) VALUES (%s, %s, %s)"
        core_db.connect().execute(insert_statment, data) """