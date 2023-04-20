"""
This script corrects mismatched transcript and translation_attr
tables for a given species which is failing the RunXrefCriticalDatacheck
of the XrefUpdate pipeline. Any transcripts with canonical_translation_id
that doesn't show up in the translation table have their
canonical_translation_id set to NULL. Any translation_attr rows that
don't have a matching entry in the translation table are dropped.
However, this is only done if these transcripts are classed as
nontranslating CDS (via the biotype column in transcript table).
"""

from argparse import ArgumentParser
import re
import pandas as pd
from ProductionMysql import *
import logging

def get_args():
	parser = ArgumentParser()
	parser.add_argument('-s', "--species", required = True,
		help = "The name of the species to correct core_db for.")
	parser.add_argument('-l', "--log_file", required = True,
		help = "Log file to write record of changes to.")
	args = parser.parse_args()
	return args

# Extends the Core class for core_db access.
# Functions for retrieving lists of unmatched transcripts
# and updating/dropping them.
class Updater(Core):
	def __init__(self, STAGING_HOST, pattern, writable = False):
		super().__init__(STAGING_HOST, pattern = pattern, writable = writable)

	# Get biotypes and ids of transcripts which don't have a matching 
	# translation entry in the translation table. 
	def get_failing_transcript_ids_and_biotype_from_transcript_table(self):
		TRANSCRIPT_SQL="SELECT t1.transcript_id, biotype FROM transcript t1 LEFT JOIN translation t2 ON t1.canonical_translation_id = t2.translation_id \
			WHERE t1.canonical_translation_id IS NOT NULL AND t2.translation_id IS NULL;"
		offending_transcripts = [x for x in self.connect().execute(TRANSCRIPT_SQL)]
		try:
			# trans_biotype_dict = {transcript_id:biotype for (transcript_id, biotype) in offending_transcripts}
			return offending_transcripts
		except IndexError:
			return None

	# Get biotypes and ids of translation_attr entries which don't have 
	# a matching translation entry in the translation table. 
	def get_failing_translation_ids_and_biotype_from_translation_attrib_table(self):
		TRANSCRIPT_SQL="SELECT t1.translation_id, t3.biotype FROM translation_attrib t1 LEFT JOIN translation t2 ON t1.translation_id = t2.translation_id \
			INNER JOIN transcript t3 ON t1.translation_id = t3.canonical_translation_id WHERE t1.translation_id IS NOT NULL AND t2.translation_id IS NULL;"
		offending_transcripts = [x for x in self.connect().execute(TRANSCRIPT_SQL)]
		try:
			# trans_biotype_dict = {transcript_id:biotype for (transcript_id, biotype) in offending_transcripts}
			return offending_transcripts
		except IndexError:
			return None


	# Set canonical_translation_id to null if transcript has no matching
	# entry in the translation table.
	def correct_canonical_translation_id_for_transcript(self, transcript_id):
		UPDATE_SQL = "UPDATE transcript SET canonical_translation_id = NULL WHERE transcript_id = '{0}'".format(transcript_id)
		self.connect().execute(UPDATE_SQL)


	# Drop entry from translation_attr taable if it has not matching
	# entry in the translation table.
	def drop_translation_id_for_transcript(self, translation_id):
		UPDATE_SQL = "DELETE FROM translation_attrib WHERE translation_id = '{0}'".format(translation_id)
		self.connect().execute(UPDATE_SQL)


# Given a list of transcripts, splits the list into transcripts which
# are nontranslating CDS and transcripts which are not. Returns 2 lists
def classify_transcripts_by_biotype(transcripts):
	nontranslating_transcripts = []
	translating_transcripts = []

	for transcript in transcripts:
		transcript_id = transcript[0]
		biotype = transcript[1]
		if biotype == "nontranslating_CDS":
			nontranslating_transcripts.append(transcript_id)
		elif biotype != "nontranslating_CDS":
			translating_transcripts.append(transcript_id)
	return nontranslating_transcripts, translating_transcripts

def all_transcripts_are_nontranslating(nontranslating_transcripts, translating_transcripts):
	return len(translating_transcripts) == 0	

# Checks if any transcripts appear as both nontranslating and translating
def check_for_ambiguous_transcript_biotypes(nontranslating_transcripts, translating_transcripts):
	transcript_intersection = set(nontranslating_transcripts) & set(translating_transcripts)
	if not transcript_intersection:
		return None
	else:
		return list(transcript_intersection)



def main():
	args = get_args()

	# Create database object to read/write from species core database.
	core_db = Updater(STAGING_HOST, args.species, writable=True)
	db_name = core_db.db()

	# Set up logging
	logging.basicConfig(level=logging.INFO, filename = args.log_file, format = "%(asctime)s - %(levelname)s: %(message)s")
	logging.info("Updating nontranslating translation IDs for {0}.".format(db_name))

	nontranslating_transcripts = []
	translating_transcripts = []

	# Gettranscript entries which lack matching entry in translation table
	offending_transcripts = core_db.get_failing_transcript_ids_and_biotype_from_transcript_table()
	# Bin them in separate lists by biotype (nontranslating CDS or other)
	nontranslating_transcripts, translating_transcripts = classify_transcripts_by_biotype(offending_transcripts)

	# Get translation_attr entries which lack matching entries in the translation table.
	# Note: run this before changing any of the tables as it relies on the transcript 
	# table not being modified.
	offending_translation_ids = core_db.get_failing_translation_ids_and_biotype_from_translation_attrib_table()
	# Bin them in separate lists by biotype (nontranslating CDS or other)
	nontranslating_translation_ids, translating_translation_ids = classify_transcripts_by_biotype(offending_translation_ids)


	# Update the canonical_translation_id to NULL for any nontranslating transcripts
	if len(offending_transcripts) == 0:
		logging.info("No transcripts found to fail the Xref datacheck for this species. Not updating the transcript table.")
	elif all_transcripts_are_nontranslating(nontranslating_transcripts, translating_transcripts):
		# Update the relevant table for each transcript and log changes
		for transcript in nontranslating_transcripts:
			logging.info("Updating canonical_translation_id to NULL for transcript {transcript_id} in table transcript".format(transcript_id=transcript))
			core_db.correct_canonical_translation_id_for_transcript(transcript)
	else:
		# Get list of transcripts that are nontranslating, verify that they don't also show up as translating for some reason
		ambiguous_transcripts = check_for_ambiguous_transcript_biotypes(nontranslating_transcripts, translating_transcripts)
		nonambiguous_nontranslating_transcripts = [transcript for transcript in nontranslating_transcripts if transcript not in ambiguous_transcripts]
		
		# Update the relevant table and log changes
		for transcript in nonambiguous_nontranslating_transcripts:
			logging.info("Updating canonical_translation_id to NULL for transcript {transcript_id} in table transcript".format(transcript_id=transcript))
			core_db.correct_canonical_translation_id_for_transcript(transcript)

		# Log transcript IDs that were not updated since those will still cause the DC to fail
		for transcript in ambiguous_transcripts:
			logging.info("Transcript {transcript_id}'s canonical_translation_id was not updated because it \
				registers as having both a nontranslating and translating biotype.".format(transcript_it=transcript))

		for transcript in translating_transcripts:
			logging.info("Transcript {transcript_id}'s canonical_translation_id was not updated because it \
				registers as having a translating biotype.".format(transcript_it=transcript))
		
	# Drop the translation_attr entry for any nontranslating transcripts.
	if len(offending_translation_ids) == 0:
		logging.info("No translation_ids found to fail the Xref datacheck for this species. Not updating the translation_attr table.")

	if all_transcripts_are_nontranslating(nontranslating_translation_ids, translating_translation_ids):
		nontranslating_translation_ids = set(nontranslating_translation_ids)
		# Update the relevant table for each transcript and log changes
		for translation_id in nontranslating_translation_ids:
			logging.info("Deleting translation_id {translation_id} from table translation_attr".format(translation_id=translation_id))
			core_db.drop_translation_id_for_transcript(translation_id)
	else:
		# Get list of transcripts that are nontranslating, verify that they don't also show up as translating for some reason
		ambiguous_translation_ids = check_for_ambiguous_transcript_biotypes(nontranslating_translation_ids, translating_translation_ids)

		nonambiguous_nontranslating_translation_ids = [translation_id for translation_id in nontranslating_translation_ids if translation_id not in ambiguous_translation_ids]
		nonambiguous_nontranslating_translation_ids = set(nonambiguous_nontranslating_translation_ids)

		# Update the relevant table and log changes
		for translation_id in nonambiguous_nontranslating_translation_ids:
			logging.info("Deleting translation_id {translation_id} from table translation_attr".format(translation_id=translation_id))
			core_db.drop_translation_id_for_transcript(translation_id)

		# Log transcript IDs that were not updated since those will still cause the DC to fail
		for translation_id in ambiguous_translation_ids:
			logging.info("{translation_id}'s translation_id was not deleted because it registers as having \
			both a nontranslating and translating biotype.".format(translation_id=translation_id))

		for translation_id in translating_translation_ids:
			logging.info("{translation_id}'s translation_id was not deleted because it registers as having \
			a translating biotype.".format(translation_id=translation_id))


if __name__ == "__main__":
	main()