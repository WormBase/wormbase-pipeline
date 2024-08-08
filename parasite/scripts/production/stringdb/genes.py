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
geneOutput_path = os.path.join(PARASITE_SCRATCH, int_folder, 'outputGenes.txt')

# read in string file as pd df
stringdf = parsers.create_stringdf_from_txt(stringtxt)

### FOR GENE MATCHES ###
# load input file as pandas df
genesdf = parsers.create_genesdf_from_txt(genetxt)

# merge the input gene and string dataframes on the stable ID columns where the values match 
merged_df = pd.merge(stringdf, genesdf, on='Useful ID')

# keep only Genome project and stable id columns
# drop rows where there are duplicate genome stable IDs 
df = merged_df[['Genome project', 'Useful ID']]
df2 = df.drop_duplicates(subset=['Useful ID'])

# output as a .csv
df2.to_csv(output_path, sep='\t', index=False)
print(f'DataFrame saved as {geneOutput_path}')