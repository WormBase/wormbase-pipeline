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
parser.add_option("-y", "--year_range", dest="YEAR",
                  help="Required: Year range in the form of START-END (e.g. 2019-2022).")

(options, args) = parser.parse_args()

if not options.TAXID:
    sys.exit("No -t/--taxon_id input. please try again.")
taxon_id=str(options.TAXID)

if not options.SPECIES:
    sys.exit("No -s/--species input. please try again.")
species=str(options.SPECIES)

if not options.YEAR or "-" not in options.YEAR:
    sys.exit("No -y/--year input. please try again.")
year_range=str(options.YEAR).split("-")
from_year=min(year_range)
to_year=max(year_range)

def ena_rnaseq_by_taxon_id_to_dict(url):
    r = requests.get(url)
    r_json = r.json()
    study_dict = {}
    for sample_accession in r_json:
        if sample_accession["secondary_study_accession"] not in study_dict.keys():
            study_dict[sample_accession["secondary_study_accession"]] = []
        study_dict[sample_accession["secondary_study_accession"]].append(sample_accession)
    return(study_dict)

def ena_year_range_to_list(url):
    study_list=[]
    r = requests.get(url)
    r_json = r.json()
    study_dict = {}
    for accession in r_json:
        if accession["secondary_study_accession"]=="":
            continue
        else:
            study_list.append(accession["secondary_study_accession"])
    return(study_list)

year_range_url = ena_by_taxon_for_year_range_url.format(taxon_id, from_year, to_year)
rnaseq_url = ena_rnaseq_by_taxon_url_onlyrnaseq.format(taxon_id)

try:
    rnaseq_studies_dict = ena_rnaseq_by_taxon_id_to_dict(rnaseq_url)
    #print("studies: " + str(len(studies_dict.keys())) + "  runs: " + str(
    #    sum([len(studies_dict[x]) for x in studies_dict])))
except ValueError:
    rnaseq_studies_dict = {}

if rnaseq_studies_dict != {}:
    try:
        year_range_studies_list = ena_year_range_to_list(year_range_url)
        print("\n".join([species + "\t" + x for x in rnaseq_studies_dict if x in year_range_studies_list]))
    except ValueError:
        pass

