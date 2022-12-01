#!/hps/software/users/wormbase/parasite/shared/.pyenv/versions/p395/bin/python

'''
A simple script to parse the 'results_by_species.json' file generated
by the Geneset Datacheck pre-production pipeline https://tinyurl.com/49x27r3k.
It generates a .tsv where each row represents a core database and each column
represents a test. Values represent the number of times a datacheck failed for
a specific core database (e.g. how many genes failed).
'''

from argparse import ArgumentParser
import json 
import re
import pandas as pd

def get_args():
	parser = ArgumentParser()
	parser.add_argument('-i', "--input_json", required = True,
		help = "Input .json file containing Geneset Datacheck outputs.")
	parser.add_argument('-o', "--output_file", required = True,
		help = "Output .tsv file to write summary of datachecks to.")
	args = parser.parse_args()
	return args 

'''
Takes the JSON object for a single species database and extracts the name
of any datachecks it failed along with the number of times it failed said 
datachecks. Returns a dictionary where the key is the name of the test
and the value is the number of times it was failed.
'''
def parse_failed_tests(json_entry):
	failed_tests = {}

	# Datachecks are lumped together into classes so that's the first 
	# level of iteration (e.g "ProteinTranslation", "Boundaries")
	for datacheck_category in json_entry:
		# Iterate over each test and extract number of failures.
		for test in json_entry[datacheck_category]["tests"]:
			# Remove constant text to obtain name of test only.
			trimmed_test = re.sub('not ok [0-9] - ', "", test)
			num_fails = json_entry[datacheck_category]["tests"][test][2]
			# Remove non-numerical fluff.
			num_fails = re.search(r'(?<=(got:)).*', num_fails)
			num_fails = re.sub("'","",num_fails[0]).strip()
			failed_tests[trimmed_test] = num_fails
	return failed_tests


def main():
	args = get_args()

	overall_df = None

	with open(args.input_json) as f:
		datacheck_json = json.load(f)

	# Iterate over species in the datacheck .json and extract
	# the names and number of failed tests for each.
	for key in datacheck_json:
		# Extract species name to use as index in dataframe.
		species = re.split(',', key)[0].strip()

		failed_tests = parse_failed_tests(datacheck_json[key])

		species_df = pd.DataFrame(failed_tests, index = [species])
		# Create dataframe if it's blank, append to it if not.
		if overall_df is None:
			overall_df = species_df
		else:
			overall_df = overall_df.append(species_df)
	overall_df = overall_df.fillna(0)
	overall_df.to_csv(args.output_file, sep = '\t')



if __name__=="__main__":
	main()