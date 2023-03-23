import os
import json
from ProductionMysql import *
from ProductionUtils import *

# Set input and output directories based on environment variables
PARASITE_SCRATCH = os.getenv("PARASITE_SCRATCH")
PARASITE_VERSION = str(os.getenv("PARASITE_VERSION"))
PARASITE_CONF = os.getenv("PARASITE_CONF")
BRC4_DIR = os.path.join(PARASITE_SCRATCH,"brc4rnaseq",f"WBPS{PARASITE_VERSION}")
GENEXP_DIR = os.path.join(BRC4_DIR,"gene_expression")
COMPONENT_DIR = os.path.join(BRC4_DIR,"component")
ASSEMBLY_JSON_OUTPUT = os.path.join(PARASITE_CONF,"brc4_rnaseq.assemblies-rename.json")
# Create the output directory if it does not exist
if not os.path.exists(GENEXP_DIR):
    os.makedirs(GENEXP_DIR)

# Define an empty dictionary to hold the mapping of original assembly names to new valid assembly names
assembly_rename_dict = {}

# Iterate over the genomes in the input directory
genomes = []
for genome in os.listdir(COMPONENT_DIR):
    # Check if the genome is a staging genome and append it to the list of genomes if it is
    if staging.is_genome(genome):
        genomes.append(genome)
    else:
        exit_with_error(f"{genome} is not a staging genome. Exiting.")
    # Extract species and assembly information from the genome name
    species = "_".join(genome.split("_")[0:2])
    assembly = Core(staging.host, genome).assembly_default()
    assembly_valid_name = assembly
    genome_dir = os.path.join(COMPONENT_DIR, genome)

    if not os.path.isdir(genome_dir):
        print_info(f"{genome_dir} is not a path. Continue.")
        continue

    # If the assembly name is not valid, try to fix it
    if not is_valid_directory_name(assembly):
        assembly_valid_name = fix_directory_name(assembly)

    # Add the original and new assembly names to the dictionary
    if assembly not in assembly_rename_dict.keys():
        assembly_rename_dict[assembly] = assembly_valid_name
    else:
        exit_with_error(f"At least 2 databases have the same assembly name {assembly}. Exiting.")

    # Create the directory structure for the current assembly if it does not exist
    assembly_dir = os.path.join(GENEXP_DIR, species, assembly_valid_name)
    if not os.path.exists(assembly_dir):
        os.makedirs(assembly_dir)

    # Iterate over the runs for the current genome and create soft links for all files in each run directory
    for run_index in os.listdir(genome_dir):
        study = run_index.split("_")[0]
        run = run_index.split("_")[1]
        run_dir = os.path.join(genome_dir, run_index, run)
        dest_run_dir = os.path.join(assembly_dir, study, run)
        if not os.path.exists(dest_run_dir):
            os.makedirs(dest_run_dir)
        create_soft_link_for_all_dir_contents(run_dir, dest_run_dir)

with open(ASSEMBLY_JSON_OUTPUT, "w") as json_file:
    json.dump(assembly_rename_dict, json_file)




