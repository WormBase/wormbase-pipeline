import pandas as pd
import gzip

# function to input a .txt file with headers and return as a pandas df
# Read the contents of the file
def create_stringdf_from_txt(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    columns = ['#string_protein_id', 'Useful ID', 'source']
    print(columns)
    # for each line in the file, create list split by tab separated values
    data = [line.strip().split('\t') for line in lines[0:]]
    # Create the DataFrame
    df = pd.DataFrame(data, columns=columns)
    return df

def create_genesdf_from_txt(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    # Extract column headers from the first line
    #columns = lines[0].strip().split('\t')
    columns = ['Genome project',  'Useful ID']
    print(columns)
    # for each line in the file, create list split by tab separated values
    data = [line.strip().split('\t') for line in lines[0:]]
    # Create the DataFrame
    # drop initial file headers from row 0
    df = pd.DataFrame(data, columns=columns).drop(0)
    return df

def create_wbpsproteindf_from_txt(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    columns = ['Genome name',  'Useful ID']
    print(columns)
    # for each line in the file, create list split by tab separated values
    data = [line.strip().split('\t') for line in lines[0:]]
    # Create the DataFrame
    # Drops original file headers in row 0
    df = pd.DataFrame(data, columns=columns).drop(0)
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


# function to iterate through output matching proteins file and for each genome project, print the total number of matches 
# input is a variable that has the full path to the output proteins file

def count_matches(output_file, number_of_matches_file):
    data_dict = {}

    # read in output genes file
    with open(output_file, 'r') as file:
        # iterate through each line in the file
        for line in file:
            # create variable columns, and assign column 1 and column 2 in the tsv to it
            # should be two tab separated columns
            columns = line.strip().split('\t')
            if len(columns) == 2:
                column1, column2 = columns
                # if column 1 is already in the dictionary, add the column 2 value to it
                # else create a new key in the dictionary and add column 2 as the values
                if column1 in data_dict:
                    data_dict[column1].append(column2)
                else:
                    data_dict[column1] = [column2]

    # Print key value pairs
    for key, values in data_dict.items():
        print(f'{key}\t{", ".join(values)}')

    # count how many gene matches there are for each genome project
    for key, value in data_dict.items():
        # print key and number of associated values
        print(key, len([item for item in value if item]))

    # count how many gene matches there are for each genome project
    with open(number_of_matches_file, 'a') as out_file:
        for key, value in data_dict.items():
            # write key and number of associated values to the output TSV file
            out_file.write(f'{key}\t{len(value)}\n')