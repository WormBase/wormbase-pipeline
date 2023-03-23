import os
import subprocess
import json
import requests
from itertools import islice
import argparse
from ProductionUtils import exit_with_error,print_info
from ProductionMysql import *

PARASITE_CONF=os.environ["PARASITE_CONF"]
PARASITE_SCRATCH=os.environ["PARASITE_SCRATCH"]
PARASITE_VERSION=os.environ["PARASITE_VERSION"]
DEFAULT_SFILE=os.path.join(PARASITE_CONF,"brc4_rnaseq.pre-run.status.json")
DEFAULT_OUTPUT=os.path.join(PARASITE_SCRATCH,"brc4rnaseq","WBPS"+PARASITE_VERSION,"dataset.json")

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--status_file", dest="SFILE", required=False,
                        default=DEFAULT_SFILE,
                        help="Optional: A json file generated by the brc4_release_update_for_pipeline.py script.")
    parser.add_argument("--species", dest="SPECIES", required=False, type=str,
                        help="Optional: Species in the <genus>_<species>_<bioproject> format. "
                             "Example schistosoma_mansoni_prjea36577")
    parser.add_argument("--studies", dest="STUDIES", required=False, type=str,
                        help="Optional: Use this option to specify a specific secondary ENA study "
                             "accession to be processed (use a comma-separated list for multiple entries")
    parser.add_argument("--max_runs_per_batch", dest="BATCHSIZE", required=False, type=int,
                        help="Optional: How many runs should be included in a pipeline's batch?"
                             "Default: 0 meaning that the pipeline will run in one batch.")
    parser.add_argument("--output_basename", dest="OUTBN", required=False, type=str, default="homes/digri/final_dataset",
                        help="Required: Basename of the output json file. Example: /example/path/example")
    # parser.add_argument("--only_skipped_runs", dest="OSR", required=False, default=False, action=argparse.BooleanOptionalAction,
    #                     help="Required: Basename of the output json file. Example: /example/path/example")
    args = parser.parse_args()
    return(args)

args = get_args()

def slice_dict(data, size=0):
    if size == 0:
        size = len(data)
    it = iter(data)
    for i in range(0, len(data), size):
        yield {k: data[k] for k in islice(it, size)}

def write_output_json(run_dict, outfile):
    dict_of_final_runs = {}
    for run_index in run_dict:
        study = run_dict[run_index]["study"]
        run_name = run_dict[run_index]["run_name"]
        production_name = run_dict[run_index]["production_name"]
        assembly = run_dict[run_index]["assembly"]
        if run_index not in dict_of_final_runs.keys():
            dict_of_final_runs[run_name] = {"name": run_name,
                                                  "species": production_name,
                                                  "production_name": production_name,
                                                  "runs": [{ "name" : run_dict[run_index]["run"],
                                                             "accessions": [run_dict[run_index]["run"]]
                                                            }]
                                           }
        else:
            exit_with_error("The run {0] is not unique in your set.".format(run_index))
        with open(outfile, "w") as final:
            json.dump(list(dict_of_final_runs.values()), final, indent=2)


print_info("PRE-RUN STATUS FILE SELECTED:\n"+args.SFILE)
if args.SPECIES:
    SPECIES=args.SPECIES.split(",")
    print_info("SPECIES SPECIFIED:\n"+"\n\t".join(SPECIES))

if args.STUDIES:
    STUDIES=args.STUDIES.split(",")
    print_info("STUDIES SPECIFIED:\n"+"\n\t".join(args.STUDIES))

if args.BATCHSIZE:
    BATCHSIZE=args.BATCHSIZE
    print_info("MAXIMUM NUMBER OF RUNS PER BATCH:\n"+str(BATCHSIZE))
else: BATCHSIZE=0

if args.OUTBN.endswith(".json"):
    OUTBN=args.OUTBN.split(".json")[0]
else:
    OUTBN=args.OUTBN

OUTDIR = os.path.dirname(OUTBN)
if not os.path.exists(OUTDIR):
   os.makedirs(OUTDIR)

status_file = open(args.SFILE)
status_studies = json.load(status_file)


all_runs = []
all_skipped_runs = []
run_studies = {}
skipped_run_studies = {}
all_studies = []
all_genomes = []
for genus_species in status_studies:
    for assembly in status_studies[genus_species]:
        production_name = status_studies[genus_species][assembly]["production_name"]
        if args.SPECIES and production_name not in SPECIES:
            continue
        all_genomes.append(production_name)
        for study in status_studies[genus_species][assembly]["studies"]:
            if args.STUDIES and study not in STUDIES:
                continue
            all_studies.append(study)
            for run in status_studies[genus_species][assembly]['studies'][study]['run_accessions']:
                run_index = production_name + "_" + study + "_" + run
                run_name = study + "_" + run
                all_runs.append(run_index)
                run_studies[run_index] = {
                                    "run":             run,
                                    "run_name":        run_name,
                                    "study":           study,
                                    "production_name": production_name,
                                    "species":         production_name,
                                    "assembly":        assembly}
            for run in status_studies[genus_species][assembly]['studies'][study]['skipped_run_accessions']:
                if run not in status_studies[genus_species][assembly]['studies'][study]['run_accessions']:
                    run_index = production_name + "_" + study + "_" + run
                    run_name = study + "_" + run
                    all_skipped_runs.append(run_index)
                    skipped_run_studies[run_index] = {
                                        "run":             run,
                                        "run_name":        run_name,
                                        "study":           study,
                                        "production_name": production_name,
                                        "species":         production_name,
                                        "assembly":        assembly}

having_chunks = False
if BATCHSIZE!=0 and BATCHSIZE < len(run_studies):
    having_chunks = True

counter = 1
for chunk in slice_dict(run_studies, size=BATCHSIZE):
    if having_chunks:
        outbn = OUTBN + "_" + str(counter) + ".json"
    else:
        outbn = OUTBN + ".json"
    write_output_json(run_dict=chunk, outfile=outbn)
    counter += 1


path = "/hps/nobackup/flicek/wormbase/parasite/brc4rnaseq/WBPS18/component"
final_list = []
for genome in os.listdir(path):
    gpath = os.path.join(path,genome)
    for study_run in os.listdir(path):
        run = study_run.split("_")[1]
        name = study_run
        species = genome
        production_name = genome
        runs = [{"name": run, "accessions": [run]}]
        final_list.append({"name": name,
                           "species": species,
                           "production_name": production_name,
                           "runs": runs})


