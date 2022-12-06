import requests
import time

ena_rnaseq_by_taxon_url = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=tax_eq({0})%20AND%20(library_strategy%3D%22RNA-Seq%22%20OR%20library_strategy%3D%22OTHER%22)&fields=study_accession%2Caccession%2Cfastq_ftp%2Csecondary_study_accession%2Csecondary_sample_accession&format=json"
ena_rnaseq_by_taxon_url_onlyrnaseq = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=tax_eq({0})%20AND%20(library_strategy%3D%22RNA-Seq%22)&fields=study_accession%2Caccession%2Cfastq_ftp%2Csecondary_study_accession%2Csecondary_sample_accession&format=json"
ena_isoseq_by_taxon_url = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_study&query=tax_eq({0})%20AND%20(instrument_platform%3D%22OXFORD_NANOPORE%22%20OR%20instrument_platform%3D%22PACBIO_SMRT%22)%20AND%20(%20library_strategy%3D%22RNA-Seq%22%20)&format=json"


from optparse import OptionParser

parser = OptionParser(usage='usage: %prog [options] arguments')
parser.add_option("-t", "--taxon_id", dest="TAXID",
                  help="Required: NCBI Taxon ID for the species of interest.")
parser.add_option("-s", "--species", dest="SPECIES",
                  help="Required: genus_species")
(options, args) = parser.parse_args()

if not options.TAXID:
    sys.exit("No -t/--taxon_id input. please try again.")
taxon_id=str(options.TAXID)

if not options.SPECIES:
    sys.exit("No -s/--species input. please try again.")
species=str(options.SPECIES)



def ena_rnaseq_by_taxon_id(url):
    r = requests.get(url)
    r_json = r.json()
    return(r_json)



url = ena_isoseq_by_taxon_url.format(taxon_id)
try:
    studies_dict_list = ena_rnaseq_by_taxon_id(url)
except ValueError:
    studies_dict_list = []


if studies_dict_list:
    for x in studies_dict_list:
        print(species+"\t"+taxon_id+"\t"+x['study_accession']+"\t"+x["description"])