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
        #print(organism['name'])
        print(organism['taxon_id'])

if __name__ == "__main__":
    main()


# grep command to run to pre-filter file from string db
# grep -E '^(6337|51028|6184|27835|6290|31234|1147741|6239|29170|334426|6188|147828|6238|6182|42157|6210|34508|157069|144512|387005|135651|7209|6335|281687|6265|70415|29172|6248|6293|51031|42155|6326|282301|6313|48269|6280|6186|174720|268474|6277|79923|6279|6198|45882|36087|1561998|46835|53326|70667|42156|6183|131310|6211|60517|6216|53468|6334|6185|31246|34506|6282|318479|1503980|268475|990121|54126|451379|103827|6205|2018661|1611254|6336|75913)' protein.aliases.v12.0.txt > filtered_protein.txt