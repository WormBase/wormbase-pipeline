# import libraries
import parsers
import os
import utils
import pandas as pd

# load environments
PARASITE_SCRATCH = os.environ["PARASITE_SCRATCH"]

# create path to .txt files in scratch directory
int_folder = 'stringdb'

genetxt = os.path.join(PARASITE_SCRATCH, int_folder, 'random_1000_mart.txt')
stringtxt = os.path.join(PARASITE_SCRATCH, int_folder, 'filtered_protein.txt')

# create path to output file containing matching genes
output_path = os.path.join(PARASITE_SCRATCH, int_folder, 'outputGenes.txt')
print(output_path)


# load .txt files as pandas dfs
stringdf = parsers.create_stringdf_from_txt(stringtxt)
genesdf = parsers.create_genesdf_from_txt(genetxt)

# store second column of each df as a list
wbpsgenes = genesdf['Gene stable ID'].to_list()
stringgenes = stringdf['alias'].to_list()

# another way to find common genes
matching_genes = utils.find_matches(wbpsgenes, stringgenes)
print(matching_genes)

# create a blank dataframe matching format of genesdf to store the reults
result_df = pd.DataFrame(columns=genesdf.columns)

# iterate through each row in the genesdf 
for index, row in genesdf.iterrows():
    if row['Gene stable ID'] in matching_genes:
        result_df = result_df.append(row)

# Print the resulting DataFrame
print(result_df)

# save results to output .txt file 
result_df.to_csv(output_path, sep='\t', index=False)
print(f'DataFrame saved as {output_path}')