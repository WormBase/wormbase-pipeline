from argparse import ArgumentParser
import json 
import re
import pandas as pd
from ProductionMysql import *
import logging

'''
Takes a WBPS core species and a list of transcripts which failed a
protein translation datacheck (nontranslating transcripts). Updates
the transcript entries in the database to indicate that they are 
nontranslating. This includes removing the transcript's corresponding
entry in the translation table and changing its biotype in the transcript
table to 'nontranslating_CDS'. 
Also checks parent genes for nontranslating transcripts. If all of the
parent gene's transcripts are nontranslating, the gene's biotype in the
gene table is also changed to 'nontranslating_CDS'.
Input paramters:
	species: The coreDB species name (genus_species_bioproject).
	transcripts_file: a flat file list of nontranslating transcripts'
	stable IDs.
	log_file: File to log DB changes in.
'''


def get_args():
	parser = ArgumentParser()
	parser.add_argument('-s', "--species", required = True,
		help = "The name of the species to extract nontranslating transcripts for ('genus_species_bioprojectID')")
	parser.add_argument('-t', '--transcripts_file', required = True,
		help = "A .tsv containing a list of nontranslating transcripts.")
	parser.add_argument('-l', "--log_file", required = True,
		help = "Log file to write record of changes to.")
	args = parser.parse_args()
	return args


'''
An extension of the Core class for making changes to the coreDB
via SQL commands.
TODO: May make more sense to insert these functions straight into the 
Core class in ProductionMysql
'''
class Updater(Core):
	def __init__(self, STAGING_HOST, pattern, writable = False):
		super().__init__(STAGING_HOST, pattern = pattern, writable = writable)

	# Given a transcript stable ID, retrieves the primary transcript ID
	def get_id_for_transcript(self, transcript_stable_id):
		TRANSCRIPT_SQL="SELECT transcript_id FROM transcript WHERE stable_id=\"{0}\";".format(transcript_stable_id)
		transcript_id = [x for x in self.connect().execute(TRANSCRIPT_SQL)]
		try:
			transcript_id = transcript_id[0][0]
			return str(transcript_id)
		except IndexError:
			return None 

	# Given a transcript ID, returns its corresponding translation ID
	def get_translation_id_for_transcript(self, transcript_id):
		TRANSLATION_SQL="SELECT translation_id FROM translation WHERE transcript_id=\"{0}\"".format(transcript_id)
		translation_id = [x for x in self.connect().execute(TRANSLATION_SQL)]
		try:
			translation_id = translation_id[0][0]
			return str(translation_id)
		except IndexError:
			return None 


	# TODO: Compress the below functions into one function that just 
	# returns the first query result.

	# Given a transcript ID, returns its parent gene ID
	def get_gene_id_for_transript(self, transcript_id):
		GENE_SQL="SELECT gene_id FROM transcript WHERE transcript_id=\"{0}\"".format(transcript_id)
		gene_id = [x for x in self.connect().execute(GENE_SQL)]
		try:
			gene_id = gene_id[0][0]
			return str(gene_id)
		except IndexError:
			return None

	# Given a transcript ID, returns its biotype
	def get_biotype_for_transcript(self, transcript_id):
		BIOTYPE_SQL="SELECT biotype FROM transcript WHERE transcript_id=\"{0}\"".format(transcript_id)
		biotype = [x for x in self.connect().execute(BIOTYPE_SQL)]
		try:
			biotype = biotype[0][0]
			return str(biotype)
		except IndexError:
			return None 

	# Given a gene ID, returns a list of unique biotypes for its transcripts
	def get_child_transcripts_biotypes_for_gene(self, gene_id):
		CHILD_TRANSCRIPT_BIOTYPES_SQL="SELECT DISTINCT biotype FROM transcript WHERE gene_id=\"{0}\"".format(gene_id)
		child_transcript_biotypes = [x for x in self.connect().execute(CHILD_TRANSCRIPT_BIOTYPES_SQL)]
		try:
			# Get list of biotypes (first element in tuple)
			child_transcript_biotypes = [x[0] for x in child_transcript_biotypes]
			return child_transcript_biotypes
		except IndexError:
			return None 

	# Given a transcript ID, drops its corresponding translation entry
	def delete_translation_row(self, transcription_id):
		DELETE_SQL="DELETE FROM translation WHERE transcription_id=\"{0}\"".format(transcription_id)
		self.connect().execute(DELETE_SQL)

	# Given a transcript ID, updates its biotype to a provided value
	def update_transcript_biotype(self, biotype, transcript_id):
		UPDATE_SQL="UPDATE transcript SET biotype=\"{0}\" WHERE transcript_id=\"{1}\"".format(biotype, transcript_id)
		self.connect().execute(UPDATE_SQL)

	# Given a gene ID, updates its biotype to a provided value
	def update_gene_biotype(self, biotype, gene_id):
		UPDATE_SQL="UPDATE gene SET biotype=\"{0}\" WHERE gene_id=\"{1}\"".format(biotype, gene_id)
		self.connect().execute(UPDATE_SQL)


def main():
	args = get_args()

	# Set up logging
	logging.basicConfig(level=logging.INFO, filename = args.log_file, format = "%(asctime)s - %(levelname)s: %(message)s")
	logging.info("Updating nontranslating transcripts for {0}.".format(args.species))

	# Load in list of nontranslating transcripts
	nontranslating_transcripts = pd.read_csv(args.transcripts_file, sep = '\t', usecols = [1], header = 0)

	# Create database object to read/write from species core database.
	core_db = Updater(STAGING_HOST, args.species, writable=True)

	# Update core database entries for each transcript
	for transcript in nontranslating_transcripts["Transcript"]:
		transcript_id = core_db.get_id_for_transcript(transcript)
		translation_id = core_db.get_translation_id_for_transcript(transcript_id)

		# Delete translation entry for transcript
		logging.info("Deleting translation ID {0} from translation table, corresponding to transcript ID {1} (stable ID {2})".format(translation_id, transcript_id, transcript))
		core_db.delete_translation_row(translation_id)
		# Update biotype to 'nontranslating_CDs' for transcript
		logging.info("Setting biotype for transcript ID {0} (stable ID {1}) to 'nontranslating_CDS'".format(transcript_id, transcript))
		core_db.update_transcript_biotype("nontranslating_CDS", transcript_id)

		gene_id = core_db.get_gene_id_for_transript(transcript_id)

		gene_child_biotypes = core_db.get_child_transcripts_biotypes_for_gene(gene_id)
		if len(gene_child_biotypes) == 1:
			if gene_child_biotypes[0] == "nontranslating_CDS":
				core_db.update_gene_biotype("nontranslating_CDS", gene_id)
				logging.info("Gene ID {0}'s transcripts are all nontranslating_CDS. Setting biotype for gene ID {0} to 'nontranslating_CDS'".format(gene_id))
		else:
			logging.info("Gene ID {0}'s transcripts have more than one biotype. No changes made its biotype.".format(gene_id))




if __name__ == "__main__":
	main()