import os
import subprocess
import json
import requests
from optparse import OptionParser
from ProductionUtils import exit_with_error, print_info, flatten
from ProductionMysql import *

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from Utils import is_ena_secondary_study, is_ena_run_accession

RELEASE_DIR = os.environ["PARASITE_PRODUCTION"] + "/data/releases/release{parasite_release}"
EXPRESSION_STUDIES_DIR = os.environ["EXPRESSION_CODE"] + "/studies/{species}"
PARASITE_VERSION =  os.environ["PARASITE_VERSION"]
PARASITE_CONF = os.environ["PARASITE_CONF"]
BRC4_DIR = os.environ["PARASITE_PRODUCTION"] + "/data/rnaseq/brc4rnaseq"
DEFAULT_OUTFILE = os.path.join(PARASITE_CONF,"brc4_rnaseq.pre-run.status.json")


parser = OptionParser(usage='usage: %prog [options] arguments')
parser.add_option("-r", "--release", dest="RELEASES",
                  help="Required: Parasite release that you would like to receive updates for. For multiple entries "
                       "use comma-separated list. e.g. 16,17.")
parser.add_option("-o", "--outfile", dest="OUTFILE",
                  help="Required: File where the output will be written to.")

(options, args) = parser.parse_args()


if not options.RELEASES:
    print_info("--release options was not provided. Parasite Release " + PARASITE_VERSION + " will be used instead.")
    RELEASES = [PARASITE_VERSION]
else:
    RELEASES = [x.strip() for x in options.RELEASES.split(",")]

if not options.OUTFILE:
    print_info("No output json file path was provided. The default "+DEFAULT_OUTFILE+" will be used instead. If exists, it will be overwritten.")
    OUTFILE = DEFAULT_OUTFILE
elif options.OUTFILE and not options.OUTFILE.endswith(".json"):
    exit_with_error("The provided output file is not a json file. Exiting.")
else:
    OUTFILE = options.OUTFILE

release_dirs = [RELEASE_DIR.format(parasite_release=x) for x in RELEASES]
release_genomes = list(set([y for y in flatten([os.listdir(x) for x in release_dirs]) if is_parasite_genome(y)]))
latest_release_genomes = staging.release_genomes
release_genomes_in_latest_release = [x for x in release_genomes if x in latest_release_genomes]


def get_run_accessions_from_irap_file(irap_file):
    if not os.path.exists(os.path.dirname(irap_file)):
        return []
    try:
        run_accessions = [x for x in pd.read_csv(irap_file, sep="\t")["Run"].to_list()]
        return run_accessions
    except FileNotFoundError:
        if irap_file.endswith("skipped_runs.tsv"):
            return[]
        else:
            exit_with_error("Could not find "+irap_file+". Exiting.")

def brc4_status(genus_species, assembly, irap_study_ids, brc4_dir=BRC4_DIR):
    genome_brc4_dir = os.path.join(brc4_dir,genus_species,assembly)
    if not os.path.exists(genome_brc4_dir) or irap_study_ids == []:
        return False
    for study_dir in [os.path.join(genome_brc4_dir,x) for x in irap_study_ids]:
        if not os.path.exists(study_dir):
            return False
    return True

genomes_rnaseq_dict = {}
count=0
for genome in release_genomes_in_latest_release:
    print_info(genome)
    genus_species = "_".join(genome.split("_")[0:2])
    assembly = Core(staging.host, genome).meta_value("assembly.default")
    irap_studies_dir = EXPRESSION_STUDIES_DIR.format(species=genus_species)
    if not os.path.exists(irap_studies_dir):
        irap_study_ids=[]
    else:
        irap_study_ids = [os.path.basename(x) for x in os.listdir(irap_studies_dir)]
    if genus_species not in genomes_rnaseq_dict.keys():
        genomes_rnaseq_dict[genus_species] = {}
    if assembly not in genomes_rnaseq_dict[genus_species].keys():
        genomes_rnaseq_dict[genus_species][assembly] = {}
    genomes_rnaseq_dict[genus_species][assembly] = { "genus_species": genus_species,
                            "assembly": assembly,
                            "production_name": genome,
                            "irap": (True if os.path.exists(irap_studies_dir) else False),
                            "brc4": brc4_status(genus_species, assembly, irap_study_ids),
                            "studies": { x: {"run_accessions":  get_run_accessions_from_irap_file(
                                os.path.join(irap_studies_dir, x, x + ".sources.tsv")),
                                "path": os.path.join(irap_studies_dir, x),
                                "skipped_run_accessions": get_run_accessions_from_irap_file(
                                    os.path.join(irap_studies_dir, x, x + ".skipped_runs.tsv"))}
                                for x in irap_study_ids}
                                                     }
    count+=1

with open(OUTFILE, "w") as outfile:
    json.dump(genomes_rnaseq_dict, outfile, indent=4)

