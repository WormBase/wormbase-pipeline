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
	def get_taxid_for_species(self, species):
		TAXID_SQL="SELECT taxon_id FROM ncbi_taxa_name WHERE name=\"{0}\";".format(species)
		taxid = [x for x in self.connect().execute(TAXID_SQL)]
		try:
			taxid = taxid[0][0]
			return str(taxid)
		except IndexError:
			return None
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
	logging.basicConfig(level=logging.INFO, filename = args.log_file, format = "%(asctime)s - %(levelname)s: %(message)s")

	tax = Taxonomy(STAGING_HOST)

	staging = Staging(STAGING_HOST) 

	missing_species = []

	for db in staging.core_dbs("_core_{0}_{1}".format(PARASITE_VERSION, ENSEMBL_VERSION)):
		core_db = Core(STAGING_HOST+"-w", db)
		wbps_species_name = core_db.meta_value("species.scientific_name")

		wbps_taxid = core_db.meta_value("species.taxonomy_id")
		ncbi_taxid = tax.get_taxid_for_species(wbps_species_name)
	
		if ncbi_taxid is None:
			ncbi_species_name = tax.get_species_for_taxid(wbps_taxid)
			if ncbi_species_name is None:
				logging.warning("Species not found in NCBI Taxonomy database: {0}. WBPS TaxID {1} also not found.".format(wbps_species_name, wbps_taxid))
			else:
				logging.warning("Species not found in NCBI Taxonomy database: {0}. WBPS TaxID {1} corresponds to following species in NCBI Taxonomy database: {2}".format(wbps_species_name, wbps_taxid, ncbi_species_name))
			continue
		if wbps_taxid != ncbi_taxid:
			logging.info("Updated TaxID for {0}: Old TaxID: {1}; Updated TaxID: {2}".format(wbps_species_name, wbps_taxid, ncbi_taxid))
			UPDATE_TAXID_SQL = "UPDATE meta SET meta_value = {0} WHERE meta_key = \"{1}\"".format(ncbi_taxid, "species.taxonomy_id")
			print(update_taxid_sql)
			core_db.connect().execute(update_taxid_sql)


if __name__ == "__main__":
	main()