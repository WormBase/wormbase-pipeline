from __future__ import print_function
import yaml
import sys
import argparse
from os.path import isdir, join
from os import path
from os import listdir
import pprint
import re

def main():
  old_fh, new_fh, data_directory = parse_arguments()
  logic_names = parse_and_print_previous_analysis_descs(old_fh, new_fh)
  new_descs   = create_new_descriptions(data_directory, logic_names, new_fh)

##

def parse_arguments():
  parser = argparse.ArgumentParser(description = "Parse configs in data directory to generate new analysis descriptions")
  parser.add_argument('--old_file', help="File with existing WBPS-specific analysis descriptions")
  parser.add_argument('--new_file', help="New analysis descriptions output file")
  parser.add_argument('--data_directory', help="Path to data directory for new species")

  if not len(sys.argv) > 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

  args=parser.parse_args()
  
  try:
    old_fh = open(args.old_file,"r")
  except:
    sys.exit("Couldn't open " + args.old_file)

  try:
    new_fh = open(args.new_file,"w")
  except:
    sys.exit("Can't open file to write: " + args.new_file)

  try:
    path.exists(args.data_directory)
  except:
    sys.exit("Directory doesn't exist: " + args.data_directory)
   
  return(old_fh, new_fh, args.data_directory)

##

def parse_and_print_previous_analysis_descs(old_fh, new_fh):
  logic_names = []

  for line in old_fh:
    line = line.rstrip()
    fields = line.split("\t")
    if len(fields) != 7:
      sys.exit("Expected old_file to have 7 fields")
    logic_names.append(fields[0])
    print(line, file=new_fh)
  return(logic_names)

##

def create_new_descriptions(data_directory, logic_names, new_fh):
  # locate all species config files
  configs = {}
  for i in listdir(data_directory):
    if isdir(join(data_directory, i)) and re.match("[^_]*_[^_]*_[^_]*$", i) and path.exists(join(data_directory,i,str(i+".conf"))):
      configs[i] = (join(data_directory,i,str(i+".conf")))
    else:
       print("WARNING: couldn't find config file for " + i)
      
  for i in configs:
    logic_names = parse_config_and_print_new(i, configs[i], logic_names, new_fh)

##

def parse_config_and_print_new(species, config, logic_names, new_fh):
  these_lns = {}
  with open(config, "r") as config_stream:
    config_data = yaml.safe_load(config_stream)
    try:
      these_lns["Coding"] = config_data[species]["gff_codinggene_logic"]
    except:
      sys.exit("WARNING: couldn't find a coding gene logic name for "+ species) # we minimally require a coding gene logic name
    try:  # it's ok not to have a logic name for pseudogenes and ncRNAs
      these_lns["Pseudogene"] = config_data[species]["gff_pseudogene_logic"]
      these_lns["Noncoding"] = config_data[species]["gff_noncodinggene_logic"]
    except KeyError:
      pass

  for ln_flavour, ln in these_lns.items():
    if ln in logic_names:
      print("Already have a description for "+ ln + " - skipping")
    else:
      print("Found new name "+ ln + " for " + species)
      logic_names.append(ln)
      description = generate_description(config_data, species)
      print(ln)
      print(description)
      print(new_fh)
      print_to_file(ln, description, new_fh)

  return(logic_names)

##

def generate_description(config_data, species):
  try:
    institute = config_data[species]["meta"]["provider.name"]
    institute_url = config_data[species]["meta"]["provider.url"]
  except:
    sys.exit("Couldn't find provider.name and/or provider.url in config file for "+ species)
  try:
    pm_text = config_data[species]["PM.text"]
    pmid = config_data[species]["PMID"]
  except:  # if there is no paper
    description = "Gene models produced by <a href=\"" + institute_url + "\"> " + institute + "</a>"
  else:
    description = "Gene models produced by <a href=\"" + institute_url + "\"> " + institute + " </a>, as described by <a href=\"https:\/\/pubmed.ncbi.nlm.nih.gov\/" + str(pmid) + "\">" + pm_text + "</a>"
  return(description) 
  
##

def print_to_file(ln, description, new_fh):
  web_data_id = 5115 
  vals = [ln, description, "Genes", 0, 1, web_data_id, 1]
  print(*vals, sep = "\t", file = new_fh)


##
if __name__== "__main__":
  main()   
