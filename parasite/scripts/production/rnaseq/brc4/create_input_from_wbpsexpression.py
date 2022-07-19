import os
import subprocess
import json
import requests
from optparse import OptionParser
from ProductionUtils import exit_with_error,print_info
from ProductionMysql import *

parser = OptionParser(usage='usage: %prog [options] arguments')
parser.add_option("-s", "--species", dest="SPECIES",
                  help="Required: Species in the <genus>_<species>_<bioproject> format. Example schistosoma_mansoni_prjea36577")
parser.add_option("-t", "--studies", dest="STUDIES",
                  help="Optional: Use this option to specify a specific secondary ENA study accession to be processed (use a comma-separated list for multiple entries")
(options, args) = parser.parse_args()
SPECIES = options.SPECIES.strip()


if (not options.SPECIES) or (len(staging.core_dbs(SPECIES))!=1):
    exit_with_error("--species options was either wrong or not provided. Please try again")
if options.STUDIES:
    STUDIES = options.STUDIES.strip().split(",")
    print_info("--studies options was used. Studies specified: "+", ".join(STUDIES)+".")
else:
    print_info("--studies options was not used. All WbpsExpression studies will be processed (Hopefully).")
    STUDIES = None

GENUS_SPECIES = "_".join(SPECIES.split("_")[0:2])
production_name = Core(staging.host, SPECIES).meta_value("production_name")

previous_wbps_studies_json = os.environ["PARASITE_SCRATCH"]+"/jbrowse/WBPS"+os.environ["PARASITE_VERSION"]+"/WbpsExpression/"+SPECIES+"/"+GENUS_SPECIES+".studies.json"
out_studies_json_dir = os.environ["PARASITE_SCRATCH"]+"/brc4rnaseq/e"+os.environ["RNASEQ_ENSEMBL_VERSION"]+"/"+SPECIES
out_studies_json = out_studies_json_dir + "/datasets.json"
if not os.path.exists(out_studies_json_dir):
    os.makedirs(out_studies_json_dir)
f = open(previous_wbps_studies_json)
previous_wbps_studies = json.load(f)

final_list_tojson = []
for study_dict in previous_wbps_studies:
    if STUDIES:
        if study_dict["study_id"] not in STUDIES:
            continue
    new_study_dict = {}
    new_study_dict["production_name"] = production_name
    new_study_dict["name"] = study_dict["study_id"]
    new_study_dict["species"] = SPECIES
    new_study_dict["runs"] = [{"name":x["run_id"], "accessions":[x["run_id"]]} for x in study_dict["runs"]]
    final_list_tojson.append(new_study_dict)


with open(out_studies_json, "w") as final:
    json.dump(final_list_tojson, final, indent=2)






