import requests
import time

ena_scRNASeq_by_taxon_url = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=tax_eq({0})%20AND%20library_source%3D%22TRANSCRIPTOMIC%20SINGLE%20CELL%22&fields=study_accession%2Csecondary_study_accession%2Chost%2Clibrary_source&format=json"


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



url = ena_scRNASeq_by_taxon_url.format(taxon_id)
try:
    studies_dict_list = ena_rnaseq_by_taxon_id(url)
except ValueError:
    studies_dict_list = []


if studies_dict_list:
    for x in studies_dict_list:
        print(species+"\t"+taxon_id+"\t"+x['study_accession'])