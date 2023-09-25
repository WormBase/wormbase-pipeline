import pandas as pd
import gzip

# function to input a .txt file with headers and return as a pandas df
# Read the contents of the file
def create_stringdf_from_txt(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    # Extract column headers from the first line
    #columns = lines[0].strip().split('\t')
    columns = ['#string_protein_id', 'alias', 'source']
    print(columns)
    # Process the contents (excluding the first line) to create a list of tuples
    data = [line.strip().split('\t') for line in lines[0:]]
    # Create the DataFrame
    df = pd.DataFrame(data, columns=columns)
    return df

def create_genesdf_from_txt(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    # Extract column headers from the first line
    #columns = lines[0].strip().split('\t')
    columns = ['Genome project',  'Gene stable ID']
    print(columns)
    # Process the contents (excluding the first line) to create a list of tuples
    data = [line.strip().split('\t') for line in lines[0:]]
    # Create the DataFrame
    df = pd.DataFrame(data, columns=columns)
    return df

# read in gzipped txt file as a pandas dataframe
def create_df_from_txt_gz(filename):
    with gzip.open(filename, 'r') as file:  # 'rt' mode reads as text
        lines = file.readlines()
    # Extract column headers from the first line
    columns = lines[0].strip().split('\t')
    # Process the contents (excluding the first line) to create a list of tuples
    data = [line.strip().split('\t') for line in lines[1:]]
    # Create the DataFrame
    df = pd.DataFrame(data, columns=columns)
    return df

