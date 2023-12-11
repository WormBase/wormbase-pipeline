import os.path
from Bio import Phylo
from ProductionUtils import *
from ProductionMysql import *

staging_cores = staging.core_databases
previous_staging_cores = previous_staging.core_databases
all_staging_cores = staging_cores + previous_staging_cores
staging_genomes = staging.release_genomes
previous_staging_genomes = previous_staging.release_genomes
all_staging_genomes = staging_genomes + previous_staging_genomes

def downloadable_core(core_db, all_staging_cores, parasite_data_dir):
    core_in_staging(core_db, error=True)

def ask_for_any_core_db_url(core_db, occasion=""):
    input_url = input("If you know it, please provide a valid {0} url for {1}:".format(occasion, core_db))
    return input_url.strip()

def not_in_any_staging_error(core_db):
    exit_with_error("{0} is not in the current or previous staging server".format(core_db))

def core_in_staging(core_db, all_staging_cores=all_staging_cores, error=False):
    if error:
        if core_db not in all_staging_cores:
            not_in_any_staging_error(core_db)
    return core_db in all_staging_cores

def core_what_staging(core_db, staging_cores=staging_cores, previous_staging_cores=previous_staging_cores):
    if core_in_staging(error=True):
        if cdb_input in staging_cores:
            return "current"
        if cdb_input in previous_staging_cores:
            return "previous"
    else:
        not_in_any_staging_error(core_db)

def cactus_argparse_core_db(cdb_input, staging_server=staging, previous_staging_server=previous_staging):
    return argparse_core_db(cdb_input, staging_server, previous_staging_server)

def cactus_argparse_fasta_file(fasta_input):
    return argparse_fasta_file(fasta_input)

def updated_or_not_updated_core_dbs(core_dbs, parasite_data_dir, staging_server):
    core_dbs_to_align_notupdated = []
    core_dbs_to_align_updated = []
    for core_db in core_dbs:
        if core_db in parasite_release_data_cores(parasite_data_dir, staging_server):
            core_dbs_to_align_updated.append(core_db)
        elif core_db not in parasite_release_data_cores(parasite_data_dir, staging_server):
            core_dbs_to_align_notupdated.append(core_db)
        else:
            exit_with_error("Cannot identify is this {0} has been updated for WBPS{1} or not".format(core_db,PARASITE_VERSION))
    return core_dbs_to_align_notupdated, core_dbs_to_align_updated

def core2download_url(core_db, parasite_version, staging_server, parasite_ftp_url, suffix):
    core = Core(staging_server.host, core_db)
    fasta_url = os.path.join(parasite_ftp_url, "{0}WBPS{1}".format(core.ftp_filename_n_filename_without_version(),parasite_version) + suffix)
    counter=0
    while not url_file_exists(fasta_url):
        print_warning("{0} is not a valid url".format(fasta_url))
        ask_for_any_core_db_url(core_db, suffix)
        counter+=1
        if counter>=3:
            exit_with_error("Sorry. Could not figure out a proper url for your core db. Exiting.")
    return fasta_url
    
def core_db_to_align_notupdated_to_fasta_args(core_db, parasite_version, previous_parasite_version, staging_server, parasite_ftp_url, fasta_arg_delimiter, fasta_suffix):
    genome = parasite_core2genome(core_db)
    if parasite_version=="18":
        if core_db=="panagrolaimus_sp1159_prjeb32708_core_18_108_1":
            core_db="panagrolaimus_ps1159_prjeb32708_core_18_108_1"
        elif core_db=="panagrolaimus_ju765_prjeb32708_core_18_108_1":
            core_db="propanagrolaimus_ju765_prjeb32708_core_18_108_1"
    prev_core_db = parasite_core2previouscore(core_db, staging_server)
    fasta_path = core2download_url(prev_core_db, previous_parasite_version, staging_server, parasite_ftp_url, fasta_suffix)
    return "{0}{1}{2}".format(genome, fasta_arg_delimiter, fasta_path)

def create_fasta_processing_command(fasta_arg, fasta_dir, fasta_arg_delimiter, parasite_ftp_url_backbone, download_bash_script, file_suffix=""):
    genome, fasta_url = fasta_arg.split(fasta_arg_delimiter)
    save_path = os.path.join(fasta_dir, genome + ("."+file_suffix+"." if file_suffix!="." and file_suffix!="" else "") +
    ("." if file_suffix=="" or file_suffix=="." else "") +
    "fa" + (".gz" if fasta_url.endswith(".gz") else ""))
    fasta_name_in_checksums = re.sub(parasite_ftp_url_backbone+"/WBPS\d+/", "",fasta_url)
    checksums_url = re.match(parasite_ftp_url_backbone+"/WBPS\d+/", fasta_url).group(0) + "CHECKSUMS"
    # if not url_file_exists(fasta_url): exit_with_error("{0} is not a valid url.".format(fasta_url))
    if not url_file_exists(checksums_url): exit_with_error("{0} is not a valid url.".format(checksums_url))
    return(" ".join(["sh",download_bash_script,fasta_url,save_path,checksums_url,fasta_name_in_checksums]),save_path)

def create_core_dumping_command(core_db_arg, fasta_dir, file_suffix="",
                                genome_dump=True, softmask_option=True,
                                protein_dump=False, canonical_only=True):
    core_db, host = core_db_arg.split(":")
    genome = parasite_core2genome(core_db)
    save_path = os.path.join(fasta_dir, genome + 
    ("."+file_suffix+"." if file_suffix!="." and file_suffix!="" else "") +
    ("." if file_suffix=="" or file_suffix=="." else "") + "fa")
    core = Core(host, core_db, writable=False)
    if genome_dump == True and protein_dump == True:
        exit_with_error("Error: Either genome_dump or protein_dump should be used.")
    elif genome_dump == True and protein_dump == False:
        return(core.dump_genome_command(save_path, softmask=softmask_option), save_path)
    elif genome_dump == False and protein_dump == True:
        return(core.dump_protein_command(save_path, canonical_only=canonical_only), save_path)

def read_newick_tree_file(tree_path):
    tree = Phylo.read(tree_path, "newick")
    return tree

def tree_clades(tree):
    return [clade.name for clade in tree.find_clades() if clade.name]

def write_seqfile_tree(tree, outfile, tree_format):
    output_buffer = StringIO()
    # Write the tree to the StringIO object
    Phylo.write(tree, output_buffer, tree_format)
    # Get the content as a string
    output_string = output_buffer.getvalue()
    # Close the StringIO object
    output_buffer.close()
    # For example, remove the undesired characters
    output_string = output_string.replace("[&r]", "")
    with open(outfile, "w") as file:
        file.write(output_string)


def write_seqfile(tree, genomes_dict, outfile, fasta_dir, tree_format="newick"):
    write_seqfile_tree(tree, outfile, tree_format)
    with open(outfile, 'a') as fp:
        for genome in genomes_dict:
            fasta_file = genomes_dict[genome]
            uncompressed_fasta_file = fasta_file.rstrip(".gz")
            if os.path.isfile(fasta_file):
                file_to_write = fasta_file
            elif os.path.isfile(uncompressed_fasta_file):
                file_to_write = uncompressed_fasta_file
            else:
                exit_with_error("Fasta file does not exist: {0} or {1}".format(fasta_file,uncompressed_fasta_file))
            
            # write each item on a new line
            fp.write(f"{genome}\t{file_to_write}\n")

