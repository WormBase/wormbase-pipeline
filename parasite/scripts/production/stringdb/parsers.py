import pandas as pd
import gzip

# function to input a .txt file with headers and return as a pandas df
# Read the contents of the file
def create_stringdf_from_txt(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    # Extract column headers from the first line
    #columns = lines[0].strip().split('\t')
    columns = ['#string_protein_id', 'Gene stable ID', 'source']
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

# 
# create a blank dataframe matching format of genesdf to store the reults

def create_df(df, match_column, match_list, output_path):
    # Create an empty DataFrame with the same columns as the input DataFrame
    result_df = pd.DataFrame(columns=df.columns)

    # Iterate through each row in the input DataFrame
    for index, row in df.iterrows():
        # Check if the value in the specified match_column is in the match_list
        if row[match_column] in match_list:
            # Append the matching row to the result DataFrame
            result_df = result_df.append(row)

    # Print the resulting DataFrame
    print(result_df)

    # Save the results to an output .txt file
    df.to_csv(output_path, sep='\t', index=False)
    print(f'DataFrame saved as {output_path}')

