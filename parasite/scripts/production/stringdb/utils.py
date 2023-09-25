import requests
import pandas as pd
import ssl

# Constants
WBPS_SERVER = "https://parasite.wormbase.org"
WBPS_SPECIES_ENDPOINT = "/rest/info/species?division=EnsemblParasite"
STRING_URL = 'https://stringdb-static.org/download/species.v11.5.txt'

# function to print first n elements of list 
def head(lst, n):  
    return lst[:n]

# Fetch species names and taxon_ids from WBPS API, process them and return as list
def fetch_wbps_species_names_taxonids():
    try:
        response = requests.get(WBPS_SERVER + WBPS_SPECIES_ENDPOINT, headers={"Content-Type": "application/json", "Accept": ""})
        response.raise_for_status()
        decoded = response.json()

        species_data = [{"name": species['name'], "taxon_id": species['taxon_id']} for species in decoded['species']]
        unique_species_data = []

        for data in species_data:
            # replace underscores with spaces to match STRING name format
            name = data['name'][:data['name'].rfind('_')].lower().replace('_', ' ')
            taxon_id = data['taxon_id']
            # Check if the name is already in the list, and only append if it's not
            if not any(item['name'] == name for item in unique_species_data):
                unique_species_data.append({"name": name, "taxon_id": taxon_id})

        return unique_species_data
    except requests.exceptions.RequestException as e:
        print(f"Error fetching data from WBPS: {e}")
        return []

# fetch organisms from STRING
# process them by dropping duplicates and converting names to lower case
# return organisms as list 
def fetch_string_organisms():
    try:
        ssl._create_default_https_context = ssl._create_unverified_context
        df = pd.read_csv(STRING_URL, sep='\t')
        organisms = df['official_name_NCBI'].str.lower().drop_duplicates().tolist()
        return organisms
    except Exception as e:
        print(f"Error fetching data from STRING: {e}")
        return []



# find common organisms between WBPS and STRING
def find_common_organisms(wbps_species_data, string_organisms):
    common_organisms = []

    for wbps_species in wbps_species_data:
        name = wbps_species['name'].lower().replace('_', ' ') # think making lower case and replacing _ has already been done above
        taxon_id = wbps_species['taxon_id']

        if name in string_organisms:
            common_organisms.append({"name": name, "taxon_id": taxon_id})

    return common_organisms


# fix this function on monday!

def find_matches(list_a, list_b):
    # Initialize a list to store the matches
    matches = []

    # Iterate through each element in List A
    for element_a in list_a:
        # Flag to check if a match is found
        match_found = False

        # Iterate through each element in List B
        for element_b in list_b:
            # Check if element_a is a complete match or substring of element_b
            if element_a in element_b:
                match_found = True
                break  # Stop searching further in List B if a match is found

        # If a match was found, add element_a to the matches list
        if match_found:
            matches.append(element_a)

    return matches

