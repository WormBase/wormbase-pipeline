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
genetxt = os.path.join(PARASITE_SCRATCH, int_folder, 'mart_export.txt')
stringtxt = os.path.join(PARASITE_SCRATCH, int_folder, 'filtered_protein.txt')
output_path = os.path.join(PARASITE_SCRATCH, int_folder, 'outputGenes.txt')

# load input .txt files as pandas dfs
stringdf = parsers.create_stringdf_from_txt(stringtxt)
genesdf = parsers.create_genesdf_from_txt(genetxt)

# merge the dataframes on the stable ID columns where the values match 
merged_df = pd.merge(stringdf, genesdf, on='Gene stable ID')

# keep only Genome project and stable id columns
# drop rows where there are duplicate genome stable IDs 
df = merged_df[['Genome project', 'Gene stable ID']]
df2 = df.drop_duplicates(subset=['Gene stable ID'])

# output as a .csv
df2.to_csv(output_path, sep='\t', index=False)
print(f'DataFrame saved as {output_path}')



# GRAVEYARD CODE
# store second column of each df as a list
#wbpsgenes = genesdf['Gene stable ID'].to_list()
#stringgenes = stringdf['alias'].to_list()

# find common genes - either partial substring match or complete match
#matching_genes = utils.find_matches(wbpsgenes, stringgenes)
#print(matching_genes)

# store results of matching genes as a new df and then output as a .txt file 
#matching_column = 'Gene stable ID'
# call function to create output file 
#parsers.create_df(genesdf, matching_column, matching_genes, output_path