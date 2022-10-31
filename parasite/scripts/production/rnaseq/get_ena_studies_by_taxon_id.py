import requests
from ProductionUtils import *
from dependencies import *
import time

from optparse import OptionParser

parser = OptionParser(usage='usage: %prog [options] arguments')
parser.add_option("-t", "--taxon_id", dest="TAXID",
                  help="Required: NCBI Taxon ID for the species of interest.")

(options, args) = parser.parse_args()

if not options.TAXID:
    sys.exit("No -t/--taxon_id input. please try again.")
taxon_id=str(options.TAXID)





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
print(url)

studies_dict = ena_rnaseq_by_taxon_id_to_dict(url)

print("studies: "+str(len(studies_dict.keys()))+"  runs: "+str(sum([len(studies_dict[x]) for x in studies_dict])))
print(studies_dict.keys())
