from ProductionMysql import *
import logging
from argparse import ArgumentParser


def get_args():
	parser = ArgumentParser()
	parser.add_argument("log_file", help = "Log file to write changelog to.")
	args = parser.parse_args()
	return args


class Taxonomy(Staging):
	def __init__(self, STAGING_HOST):
		super().__init__(STAGING_HOST)
		self.db = "ncbi_taxonomy" 

	def connect(self):
		dbc = DBConnection(self.url+self.db)
		return dbc 

	# Given a species as <Genus species>, retrieves the corresponding
	# TaxId from the NCBI taxonomy database
	def get_taxid_for_species(self, species):
		TAXID_SQL="SELECT taxon_id FROM ncbi_taxa_name WHERE name=\"{0}\";".format(species)
		taxid = [x for x in self.connect().execute(TAXID_SQL)]
		try:
			taxid = taxid[0][0]
			return str(taxid)
		except IndexError:
			return None

	# Given a TaxID, retrieves the corresponding species from NCBI
	# taxonomy database. Species returnes as <Genus species>
	def get_species_for_taxid(self, taxid):
		SPECIES_SQL="SELECT name FROM ncbi_taxa_name WHERE taxon_id=\"{0}\";".format(taxid)
		species = [x for x in self.connect().execute(SPECIES_SQL)]
		print(species)
		try:
			species = species[0][0]
			return str(species)
		except IndexError:
			return None



def main():
	args = get_args()

	# Camel case? In my Python? Eww
	# Start logger
	logging.basicConfig(level=logging.INFO, filename = args.log_file, format = "%(asctime)s - %(levelname)s: %(message)s")

	logging.info("Running TaxID update script.")

	# Handle for ncbi taxonomy database
	tax = Taxonomy(STAGING_HOST)

	# Handle for current ParaSite staging database
	staging = Staging(STAGING_HOST) 

	# Flag to keep track of changes for output message.
	taxids_changed = False

	for db in staging.core_dbs("_core_{0}_{1}".format(PARASITE_VERSION, ENSEMBL_VERSION)):
		core_db = Core(STAGING_HOST+"-w", db)
		wbps_species_name = core_db.meta_value("species.scientific_name")

		wbps_taxid = core_db.meta_value("species.taxonomy_id")
		# Retrieve NCBI version of TaxID using species name
		ncbi_taxid = tax.get_taxid_for_species(wbps_species_name)
	
		# If no TaxID available for WBPS species name, find 
		# corresponding species name from NCBI using WBPS TaxID
		# Log the NCBI species name, if any is found.
		if ncbi_taxid is None:
			ncbi_species_name = tax.get_species_for_taxid(wbps_taxid)
			if ncbi_species_name is None:
				logging.warning("Species not found in NCBI Taxonomy database: {0}. WBPS TaxID {1} also not found.".format(wbps_species_name, wbps_taxid))
			else:
				logging.warning("Species not found in NCBI Taxonomy database: {0}. WBPS TaxID {1} corresponds to following species in NCBI Taxonomy database: {2}".format(wbps_species_name, wbps_taxid, ncbi_species_name))
			continue
		# If WBPS and NCBI TaxIDs for species do not match, update WBPS TaxID
		if wbps_taxid != ncbi_taxid:
			logging.info("Updated TaxID for {0}: Old TaxID: {1}; Updated TaxID: {2}".format(wbps_species_name, wbps_taxid, ncbi_taxid))
			UPDATE_TAXID_SQL = "UPDATE meta SET meta_value = {0} WHERE meta_key = \"{1}\"".format(ncbi_taxid, "species.taxonomy_id")
			core_db.connect().execute(update_taxid_sql)
			taxids_changed = True

	if not taxids_changed:
		logging.info("All TaxIDs were already up to date. No changes made.")

if __name__ == "__main__":
	main()