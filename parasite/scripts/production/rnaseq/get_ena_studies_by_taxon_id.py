import requests
from ProductionUtils import *
from dependencies import *
import time

from optparse import OptionParser

parser = OptionParser(usage='usage: %prog [options] arguments')
parser.add_option("-t", "--taxon_id", dest="TAXID",
                  help="Required: NCBI Taxon ID for the species of interest.")
parser.add_option("-s", "--species", dest="SPECIES",
                  help="Required: Species in the format of genus_species (e.g. schistosoma_mansoni).")

(options, args) = parser.parse_args()

if not options.TAXID:
    sys.exit("No -t/--taxon_id input. please try again.")
taxon_id=str(options.TAXID)

if not options.SPECIES:
    sys.exit("No -s/--species input. please try again.")
species=str(options.SPECIES)

def ena_rnaseq_by_taxon_id_to_dict(url):
    r = requests.get(url)
    r_json = r.json()
    study_dict = {}
    for sample_accession in r_json:
        if sample_accession["secondary_study_accession"] not in study_dict.keys():
            study_dict[sample_accession["secondary_study_accession"]] = []
        study_dict[sample_accession["secondary_study_accession"]].append(sample_accession)
    return(study_dict)

url = ena_rnaseq_by_taxon_url_onlyrnaseq.format(taxon_id)

try:
    studies_dict = ena_rnaseq_by_taxon_id_to_dict(url)
    #print("studies: " + str(len(studies_dict.keys())) + "  runs: " + str(
    #    sum([len(studies_dict[x]) for x in studies_dict])))
    print("\n".join([species+"\t"+x for x in studies_dict]))
except ValueError:
    studies_dict = {}
    print("")
