import utils

def main():
    # Fetch and process data
    wbps_species_data = utils.fetch_wbps_species_names_taxonids()
    string_organisms = utils.fetch_string_organisms()

    # Find common organisms
    common_organisms = utils.find_common_organisms(wbps_species_data, string_organisms)
    print(len(common_organisms))

    # Print common organisms with their taxon_id
    for organism in common_organisms:
        print(f"Organism: {organism['name']}, Taxon ID: {organism['taxon_id']}")


    # print the taxon IDs of the common organisms 
    # this can be used to filter the string file containing the genes
    for organism in common_organisms:
        print(organism['taxon_id'])

if __name__ == "__main__":
    main()
