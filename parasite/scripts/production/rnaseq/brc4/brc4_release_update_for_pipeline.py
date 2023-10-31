import os
import json
import glob
from optparse import OptionParser
from ProductionUtils import exit_with_error, print_info, flatten
from ProductionMysql import *

sys.path.insert(1, os.path.join(sys.path[0], '..'))
from Utils import is_ena_secondary_study, is_ena_run_accession

RELEASE_DIR = os.environ["PARASITE_PRODUCTION"] + "/data/releases/release{parasite_release}"
EXPRESSION_STUDIES_DIR = os.environ["EXPRESSION_CODE"] + "/studies/{species}"
PARASITE_VERSION =  os.environ["PARASITE_VERSION"]
PREVIOUS_PARASITE_VERSION = os.environ["PREVIOUS_PARASITE_VERSION"]
PARASITE_CONF = os.environ["PARASITE_CONF"]
BRC4_DIR = os.environ["PARASITE_PRODUCTION"] + "/data/rnaseq/brc4rnaseq"
PREVIOUS_RELEASE_OUTFILE = os.path.join(PARASITE_CONF,"..",f"WBPS{PREVIOUS_PARASITE_VERSION}","brc4_rnaseq.pre-run.status.json")
DEFAULT_OUTFILE = os.path.join(PARASITE_CONF,"brc4_rnaseq.pre-run.status.json")

parser = OptionParser(usage='usage: %prog [options] arguments')
parser.add_option("-r", "--release", dest="RELEASES",
                  help="Required: Parasite release that you would like to receive updates for. For multiple entries "
                       "use comma-separated list. e.g. 16,17.")
parser.add_option("-o", "--outfile", dest="OUTFILE",
                  help="Required: File where the output will be written to.")
parser.add_option("-s", "--studies_tsv", dest="FORCED_STUDIES_TSV", default=os.path.join(PARASITE_CONF,"forced_studies.tsv"),
                  help="Optional: Specify a TSV file with study IDs that should be realigned for this release. " +\
                    "This option is useful when you want to process studies that are not part of the current release, " +\
                    "such as new studies or updates to existing studies. " +\
                    "The TSV file should contain one study ID (e.g., ERP015574) per line. ")

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
release_genus_species = ["_".join(genome.split("_")[0:2]) for genome in release_genomes]
latest_release_genomes = staging.release_genomes
previous_release_genomes = previous_staging.release_genomes
release_genomes_in_latest_release = [x for x in release_genomes if x in latest_release_genomes]
not_updated_previous_release_genomes = [x for x in latest_release_genomes if x in previous_release_genomes and x not in release_genomes]

def extract_studies_from_file(file_path):
    studies = []

    try:
        with open(file_path, 'r') as file:
            for line in file:
                study = line.strip()
                if study:  # Check if the line is not empty
                    studies.append(study)

    except FileNotFoundError:
        print_warning(f"{file_path} does not exist. Will continue without looking for additional studies apart from the ones." +\
            f"that will be analysed for the release genomes: {', '.join(RELEASES)}.")
        return(studies)
    except Exception as e:
        print(f"An error occurred: {e}")
    
    return(studies)

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

def get_irap_study_ids_for_species(species):
    irap_studies_dir = EXPRESSION_STUDIES_DIR.format(species=species)
    if not os.path.exists(irap_studies_dir):
        irap_study_ids=[]
    else:
        irap_study_ids = [os.path.basename(x) for x in os.listdir(irap_studies_dir)]
    return(irap_study_ids)

def generate_rnaseq_dict(genome, irap_study_ids, rnaseq_dict):
    genus_species = "_".join(genome.split("_")[0:2])
    assembly = Core(staging.host, genome).meta_value("assembly.default")
    irap_studies_dir = EXPRESSION_STUDIES_DIR.format(species=genus_species)
    if genus_species not in rnaseq_dict.keys():
        rnaseq_dict[genus_species] = {}
    if assembly not in rnaseq_dict[genus_species].keys():
        rnaseq_dict[genus_species][assembly] = {}
    rnaseq_dict[genus_species][assembly] = { "genus_species": genus_species,
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
    return(rnaseq_dict)

if not options.FORCED_STUDIES_TSV:
    forced_studies = []
else:
    forced_studies = extract_studies_from_file(options.FORCED_STUDIES_TSV)

genomes_rnaseq_dict = {}
for genome in release_genomes_in_latest_release:
    print_info(genome)
    genus_species = "_".join(genome.split("_")[0:2])
    irap_study_ids = get_irap_study_ids_for_species(genus_species)
    genomes_rnaseq_dict = generate_rnaseq_dict(genome, irap_study_ids, genomes_rnaseq_dict)

if len(forced_studies)>0:
    forced_species_studies_dict = {}
    for study_id in forced_studies:
        study_paths = glob.glob(os.path.join(EXPRESSION_STUDIES_DIR.format(species=""),f"**/{study_id}"))
        if len(study_paths)!=1:
            study_paths_str = "\n".join(study_paths)
            print_error(f"Multiple paths found for {study_id}: {study_paths_str}. Please investigate, fix and re-run.")
        study_path = study_paths[0]
        genus_species = os.path.dirname(study_path).split("/")[-1]
        if genus_species not in forced_species_studies_dict.keys():
            forced_species_studies_dict[genus_species] = []
        forced_species_studies_dict[genus_species].append(study_id)

    print_info("Looking for updated studies:")
    for genus_species in forced_species_studies_dict:
        if genus_species in release_genus_species:
            continue
        studies = forced_species_studies_dict[genus_species]
        genomes = ["_".join(x.split("_")[0:3]) for x in staging.release_core_dbs(genus_species)]
        if genus_species not in ["_".join(genome.split("_")[0:2]) for genome in  latest_release_genomes]:
            print_error(f"Doesn't look like {genus_species} (from study {study}) is a species in WBPS{PARASITE_VERSION}. "+\
                "Please check the studies you specified in --studies_tsv.")
        for genome in genomes:
            print_info(genome)
            irap_study_ids = get_irap_study_ids_for_species(genus_species)
            for fstudy in studies:
                if fstudy not in irap_study_ids:
                    print_error("Cannot find the study directory for study {study} for {genus_species}.")
            studies_str = ",".join(studies)
            print_info(f"Found {studies_str}")
            irap_study_ids = studies
            genomes_rnaseq_dict = generate_rnaseq_dict(genome, irap_study_ids, genomes_rnaseq_dict)

with open(OUTFILE, "w") as outfile:
    json.dump(genomes_rnaseq_dict, outfile, indent=4)

