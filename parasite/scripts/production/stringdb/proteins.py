# import libraries
import parsers
import os
import utils
import pandas as pd

# load intermediate environment
PARASITE_SCRATCH = os.environ["PARASITE_SCRATCH"]

# create path to input and output .txt files in scratch directory
# set variable for intermediate directory 
int_folder = 'stringdb'
# create paths
stringtxt = os.path.join(PARASITE_SCRATCH, int_folder, 'filtered_protein.txt')
protOutput_path = os.path.join(PARASITE_SCRATCH, int_folder, 'outputtrEMBLProteins.txt')
proteintxt = os.path.join(PARASITE_SCRATCH, int_folder, 'uniprot_trEMBL.txt')
totalMatchesTxt = os.path.join(PARASITE_SCRATCH, int_folder,'totalMatches.txt')

# read in input string file as pd df
stringdf = parsers.create_stringdf_from_txt(stringtxt)

## load input protein file as pandas dfs
proteinsdf = parsers.create_wbpsproteindf_from_txt(proteintxt)

## merge the input protein and string df on the stable ID columns where the values match 
merged_df = pd.merge(stringdf, proteinsdf, on='Useful ID')

# keep only Genome project and stable id columns
# drop rows where there are duplicate genome stable IDs 
df = merged_df[['Genome name', 'Useful ID']]
df2 = df.drop_duplicates(subset=['Useful ID'])

# output as a .csv
df2.to_csv(protOutput_path, sep='\t', index=False)
print(f'DataFrame saved as {protOutput_path}')

### PROCESSING OUTPUT FILES ###

# process the output file by counting the total number of matches for each genome project 
outputmatches = parsers.count_matches(protOutput_path, totalMatchesTxt)
print(outputmatches)
