# import libraries
import requests, sys
import pandas as pd
import ssl

# get JSON file containing all available species, their aliases, available adaptor groups and data release in WBPS
server = "https://parasite.wormbase.org"
ext = "/rest/info/species?division=EnsemblParasite" 
 
r = requests.get(server+ext, headers={ "Content-Type" : "application/json", "Accept" : ""})
 
if not r.ok:
  r.raise_for_status()
  sys.exit()
 
decoded = r.json()

# from the .JSON file create a list that appends the species names
names = []
for species in decoded['species']:
    names.append(species['name'])

# remove the _bioproject from these names
for i in range(len(names)):
    names[i] = names[i][:names[i].rfind('_')]

# make all the names lower case and remove any '_' 
# this is to ensure the list of names is in the same format as the StringDB names
names = [x.lower() for x in names]
# replace underscores with white space
for i in range(len(names)):
    names[i] = names[i].replace('_', ' ')

# drop duplicates
names = list(dict.fromkeys(names))

# import .txt file of all the organisms hosted by STRING as a pandas dataframe
# overcome SSL error to allow the file to be imported 
ssl._create_default_https_context = ssl._create_unverified_context

# read in .txt file from url
url = 'https://stringdb-static.org/download/species.v11.5.txt'
df = pd.read_csv(url, sep='\t')

# append the organism names from the dataframe to a list
organisms = df['official_name_NCBI'].to_list()
# list comprehension to make all characters lower case to match their format in WBPS
organisms = [x.lower() for x in organisms]
# drop duplicates
organisms = list(dict.fromkeys(organisms))

# see which organisms on WBPS are also on string by comparing the two lists
# make this into a function
common_organisms = [element for element in names if element in organisms]

# print common organisms with each one separated by a new line
print(*common_organisms,sep='\n')