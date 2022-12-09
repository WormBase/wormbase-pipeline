from argparse import ArgumentParser
import json 
import re
import pandas as pd

'''
A script to parse the outputs of the Geneset datacheck and extract a list
of nontranslating transcripts for each species (based on failing the
ProteinTranslation datacheck). The list is output as flat .tsv file.
Input parameters:
	input_json: The results_by_species.json file output by the Geneset
	Datacheck. Contains list of failed tests by species.
	species: The coreDB species name (genus_species_bioproject).
	output_file: .tsv file to write transcript names to.
'''

def get_args():
	parser = ArgumentParser()
	parser.add_argument('-i', "--input_json", required = True,
		help = "Input .json file containing Geneset Datacheck outputs.")
	parser.add_argument('-s', "--species", required = True,
		help = "The name of the species to extract nontranslating transcripts for ('genus_species_bioprojectID')")
	parser.add_argument('-o', "--output_file", required = True,
		help = "Output .tsv file to write summary of datachecks to.")
	args = parser.parse_args()
	return args


# Parse nontranslating transcript names from the json entry for a species
def extract_nontranslating_transcripts(json_entry):
	# Parsing json structure. Tests are divided into categories and we
	# want the ProteinTranslation category.
	for datacheck_category in json_entry:
		if datacheck_category == "ProteinTranslation":
			# Parse transcript stable ID out of failed datacheck entry 
			for test in json_entry[datacheck_category]["tests"]:
				nontranslating_transcripts = [re.split(" ",x)[0] for x in json_entry[datacheck_category]["tests"][test] if "has invalid translation" in x]
				return nontranslating_transcripts
		else:
			continue


def main():
	args = get_args()

	# Load in tests results json
	with open(args.input_json) as f:
		datacheck_json = json.load(f)

	# Extract nontranslating transcripts for species in json.
	for key in datacheck_json:
		# Parse species name
		species=re.split(',',key)[0]
		# Skip non-target species
		if species != args.species:
			continue
		else:
			# Extract stable IDs for nontranslating transcripts.
			nontranslating_transcripts = extract_nontranslating_transcripts(datacheck_json[key])
			# Turn names into a dataframe and write to .tsv
			transcript_df = pd.DataFrame(data=nontranslating_transcripts, columns=["Transcript"])
			transcript_df.to_csv(args.output_file, sep = '\t')
			# Done, exit loop.
			break



if __name__ == "__main__":
	main()