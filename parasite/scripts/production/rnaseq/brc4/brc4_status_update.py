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
BRC4_DATA_DIR = os.environ["PARASITE_PRODUCTION"]

parser = OptionParser(usage='usage: %prog [options] arguments')
parser.add_option("-r", "--release", dest="RELEASES",
                  help="Required: Parasite release that you would like to receive updates for. For multiple enties "
                       "use comma-separated list. e.g. 16,17.")
parser.add_option("-o", "--outfile", dest="OUTFILE",
                  help="Required: File where the output will be written to.")

(options, args) = parser.parse_args()
PARASITE_RELEASE = os.environ["PARASITE_VERSION"]

if not options.RELEASES:
    print_info("--release options was not provided. Parasite Release " + PARASITE_RELEASE + " will be used instead.")
    RELEASES = [PARASITE_RELEASE]
else:
    RELEASES = [x.strip() for x in options.RELEASES.split(",")]

if not options.OUTFILE:
    exit_with_error("No outfile was provided (-o/--outfile). Exiting.")

OUTFILE=options.OUTFILE

release_dirs = [RELEASE_DIR.format(parasite_release=x) for x in RELEASES]
release_genomes = list(set([y for y in flatten([os.listdir(x) for x in release_dirs]) if is_parasite_genome(y)]))


def get_run_accessions_from_irap_file(irap_file):
    try:
        run_accessions = [x for x in pd.read_csv(irap_file, sep="\t")["Run"].to_list() if is_ena_run_accession(x)]
    except FileNotFoundError:
        if irap_file.endswith("skipped_runs.tsv"):
            run_accessions = []
        else:
            raise FileNotFoundError
    return run_accessions


f = open(OUTFILE, "w")
f.write("Genome\tirap_status\tbrc4_status\tNo_of_Studies\tStudies\tNo_of_runs\tNo_of_bigwig-only_runs\n")
f.close()

for genome in release_genomes:
    print_info(genome)
    genus_species = "_".join(genome.split("_")[0:2])
    assembly = Core(staging.host, genome).meta_value("assembly.default")
    irap_studies_dir = EXPRESSION_STUDIES_DIR.format(species=genus_species)
    if not os.path.exists(irap_studies_dir):
        f = open(OUTFILE, "a")
        f.write("\t".join([genome, "no_irap_metadata", "no_brc4", "0", "", "0", "", "\n"]))
        f.close()
        continue
    irap_study_ids = [os.path.basename(x) for x in os.listdir(irap_studies_dir) if
                      is_ena_secondary_study(os.path.basename(x))]
    rnaseq_dict[genome] = { "genus_species": genus_species,
                            "assembly": assembly,
                            "irap_status": "True",
                            "brc4_status": brc4_status(),
                            "studies": { x: {"run_accessions":  get_run_accessions_from_irap_file(
                                                os.path.join(irap_studies_dir, x, x + ".sources.tsv")),
                                             "path": os.path.join(irap_studies_dir, x),
                                             "skipped_run_accessions": get_run_accessions_from_irap_file(
                                                os.path.join(irap_studies_dir, x, x + ".skipped_runs.tsv"))}
                                        for x in irap_study_ids}
                          }


    studies = { x: {"run_accessions":  get_run_accessions_from_irap_file(
                                        os.path.join(irap_studies_dir, x, x + ".sources.tsv")),
                    "skipped_run_accessions": get_run_accessions_from_irap_file(
                                        os.path.join(irap_studies_dir, x, x + ".skipped_runs.tsv"))} for x in irap_study_ids}

                            "studies": { x: {
                                    "run_accessions": get_run_accessions_from_irap_file(
                                        os.path.join(irap_studies_dir, x, x + ".sources.tsv")),
                                    "skipped_run_accessions": get_run_accessions_from_irap_file(
                                        os.path.join(irap_studies_dir, x, x + ".skipped_runs.tsv"))
                                for x in irap_study_ids }}
    }
    all_studies = [x for x in rnaseq_dict]
    all_runs = flatten([rnaseq_dict[x]["run_accessions"] for x in rnaseq_dict])
    all_skipped_runs = flatten([rnaseq_dict[x]["skipped_run_accessions"] for x in rnaseq_dict])
    f = open(OUTFILE, "a")
    f.write("\t".join([genome, "irap_metadata", "no_brc4", str(len(all_studies)), ",".join(all_studies), str(len(all_runs)),
                       str(len(all_skipped_runs)), "\n"]))
    f.close()

