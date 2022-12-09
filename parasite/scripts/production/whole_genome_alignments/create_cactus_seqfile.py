import os
import re
import sys
import shutil
import glob
from argparse import ArgumentParser
from ProductionUtils import *
from ProductionMysql import *
from LSF import *
from Bio import Phylo
import validators

def get_args():
    parser = ArgumentParser()
    parser.add_argument("-p","--project_name", required = True, default = None,
        help = "Project name that will be used for directory naming etc.")
    parser.add_argument("-c", "--core_db", required = False, action = 'append', default = None, nargs = '+',
        help = "Core database you would like to perform a cactus alignments on.")
    parser.add_argument("-f", "--fasta_file", required = False, action = 'append', default = None, nargs = '+',
        help = "Fasta file path or URL you would like to include in the whole genome cactus alignment. You should also"
               "provide the genome name you would like to be used in the seqfile. Format: <genome_name>:<fasta_file_path>"
               "Example: homo_sapiens_38:/nfs/ftp/ensembl/homo_sapiens_38.softmask.fasta")
    parser.add_argument("-d", "--parent_directory", required = False, default=None,
        help = "Parent directory where the analysis will run and the output files will be created.")
    parser.add_argument("-t", "--tree", required = True, default="parasite",
        help = "Tree that will be used to infer tree relationships from.")

    args = parser.parse_args()

    return args

args = get_args()

parent_directory=os.path.join(os.environ.get("PARASITE_SCRATCH"),"cactus_alignments")
if args.parent_directory:
    parent_directory = args.parent_directory


project_dir=os.path.join(parent_directory, args.project_name)
print_info("PROJECT DIRECTORY:" + project_dir)

input_dir=os.path.join(project_dir, "input")
fasta_dir = os.path.join(input_dir, "fasta")
dump_genomes_log_dir=os.path.join(project_dir, "dump_genomes_log")
print_info("INPUT DIRECTORY:" + input_dir)
seqfile=os.path.join(input_dir, args.project_name + ".seqfile")

for dir in [parent_directory, project_dir, input_dir, fasta_dir, dump_genomes_log_dir]:
    if not os.path.exists(dir):
        os.makedirs(dir)

seqfile_dict={}

if args.core_db:
    print_info("Core dbs were defined. Processing...")
    core_dbs = flatten(args.core_db)
    for core_db in core_dbs:
        print_w_indent(core_db)
        server = core_which_staging(core_db)
        core = Core(server.host, core_db, writable=True)
        prod_name = core.meta_value("species.production_name")
        output_fasta_basename = prod_name+".genome.fa"
        output_fasta = os.path.join(fasta_dir, output_fasta_basename)
        dump_fasta_command = core.dump_genome_command(output_fasta, softmask=True)
        if prod_name not in seqfile_dict.keys():
            seqfile_dict[prod_name] = {"output_fasta_basename" : output_fasta_basename,
                                       "output_fasta" : output_fasta,
                                       "dump_fasta_command" : dump_fasta_command}
        else:
            exit_with_error("You want to include "+prod_name+" more than one times in the seqfile. Review your input and try again")

if args.fasta_file:
    print_info("Fasta files were defined. Processing...")
    for fasta_arg in flatten(args.fasta_file):
        print_w_indent(fasta_arg)
        if len(fasta_arg.split(":", 1))!=2:
            exit_with_error("Wrong fasta input for "+fasta_arg+". Please try again.")
        prod_name, fasta_path = [x.strip() for x in fasta_arg.split(":", 1)]
        output_fasta_basename = prod_name+".genome.fa"
        output_fasta = os.path.join(fasta_dir, output_fasta_basename)
        dump_fasta_output = output_fasta + (".gz" if fasta_path.endswith(".gz") else "")
        dump_fasta_command = download_if_url_or_copy_file_command(fasta_path, dump_fasta_output)
        if dump_fasta_output.endswith(".gz"):
            dump_fasta_command += "\n\ngzip -d " + dump_fasta_output + ";"
        if prod_name not in seqfile_dict.keys():
            seqfile_dict[prod_name] = {"output_fasta_basename" : output_fasta_basename,
                                       "output_fasta" : output_fasta,
                                       "dump_fasta_command" : dump_fasta_command}
        else:
            exit_with_error(prod_name+"  will appear more than one times in the seqfile. Review your input and try again")

print_info("Processing tree file: " + args.tree)
tree_file = args.tree
tree = Phylo.read(tree_file, "newick")
tree_clades = [clade.name for clade in tree.find_clades() if clade.name]
if set(seqfile_dict.keys()) != set(tree_clades):
    exit_with_error("The production names for the genomes you selected do not match with the ones in the tree file: "+tree_file+".")

print_info("Submitting jobs: ")
dump_job_ids=[]
for genome in seqfile_dict:
    job_id = lsf_submit(command=seqfile_dict[genome]["dump_fasta_command"],
                                   jobprefix=genome+"_dump",
                                   cpu=1, mem="1gb",
                                   cwd=dump_genomes_log_dir,
                                   queue="production",
                                   only_write=False)
    print_w_indent(job_id)
    dump_job_ids.append(job_id)


while len(dump_job_ids)!=0:
    for job_id in dump_job_ids:
        if lsf_job_status(job_id)!="RUN" and lsf_job_status(job_id)!="PEND":
            dump_job_ids.remove(job_id)
            print_info(job_id +" is DONE")
    print_info("Number of jobs still running: " + str(len(dump_job_ids)))
    print_info("Jobs still running: " + ", ".join(dump_job_ids))
    time.sleep(15)


print_info("Writing output seqfile: " + seqfile)
Phylo.write(tree, seqfile, "newick")
with open(seqfile, 'a') as fp:
    for genome in seqfile_dict:
        # write each item on a new line
        fp.write("%s\t%s\n" % (genome, seqfile_dict[genome]["output_fasta"]))
print('Done')

all_file_exist = [os.path.isfile(seqfile_dict[x]["output_fasta"]) for x in seqfile_dict]
if all(all_file_exist) and os.path.isfile(seqfile):
    shutil.rmtree(dump_genomes_log_dir, ignore_errors=False, onerror=None)
    print_info("Done.")
else:
    print_info("Error. Some files have not been created as expected.")