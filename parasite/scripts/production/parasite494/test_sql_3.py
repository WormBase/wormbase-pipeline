import os
import pandas as pd
from ProductionMysql import *
from ProductionUtils import *
from sqlalchemy import insert, create_engine, Table, MetaData
from sqlalchemy.orm import sessionmaker


STAGING_HOST = os.environ['PARASITE_STAGING_MYSQL']
databases = staging.core_dbs("core_18_108_285")


# function to create the file name for each .tsv output from functional annotations parsing script 
def get_file_path(s):
    organism_genus = s.split("_")[0][0]
    organism_species = s.split("_")[1]
    bioproject_id = s.split("_")[2].upper()
    return f"{organism_genus}_{organism_species}.{bioproject_id}.WS{WORMBASE_VERSION}.functional_descriptions.txt.gz.tsv"

# list of the file names 
organism_files = [get_file_path(x) for x in databases] 

# loop through each .tsv file and create a df with stable_id and description columns
for file in organism_files:
    df = pd.read_csv(file, sep='\t', names=["stable_id", "description"])
    stable_ids = df['stable_id'].tolist()

# connect to db associated with each .tsv file
    for database in databases:
        core_db = Core(STAGING_HOST, database)
        #print(core_db)
        # query gene table and return gene_id associated with each stable_id as a tuple
        execution = core_db.connect().execute("SELECT stable_id, gene_id FROM gene WHERE stable_id IN :stable_ids", {"stable_ids": stable_ids})
        gene_id_rows = execution.fetchall()
        # create a df from the tuples
        gene_id_df = pd.DataFrame(gene_id_rows, columns=['stable_id', 'gene_id'])
        # query attrib_type table and return description attrib_type_id
        attrib_type_query = core_db.connect().execute("SELECT attrib_type_id FROM attrib_type WHERE name = 'description'")
        attrib_type_id = attrib_type_query.fetchone()[0]
        # add attrib type id as a column to df
        gene_id_df['attrib_type_id'] = attrib_type_id
        # merge .tsv df and df created from querying db on the stable_id column 
        merged_df = pd.merge(gene_id_df, df, on='stable_id', how='inner')
        merged_df = merged_df[~merged_df['description'].str.contains('none available')]
        merged_df = merged_df.drop('stable_id', axis=1)
        merged_df = merged_df.rename(columns={'description': 'value'})
        # convert df to dictionary to make easier to insert into db?
        data = merged_df.to_dict(orient='records')
        #print(data)

        engine = create_engine(f"mysql://{STAGING_HOST}/{database}") #core.engine()
        metadata = MetaData(bind=engine)
        gene_attrib = Table('gene_attrib', metadata, autoload=True)

        # Insert data into gene_attrib table
        conn = engine.connect()
        conn.execute(gene_attrib.insert().values(data))
        conn.close()


"""         engine = create_engine(core_db)
        Session = sessionmaker(bind=engine)
        session = Session()

        metadata = MetaData(bind=engine)
        table_name = Table('gene_attrib', metadata, autoload=True)

        for row in data:
            insert_stmt = table_name.insert().values(**row)
            session.execute(insert_stmt)

        session.commit() """
