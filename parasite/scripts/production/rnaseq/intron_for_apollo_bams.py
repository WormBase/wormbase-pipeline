#!/hps/software/users/wormbase/parasite/shared/.pyenv/versions/p395/bin/python

import sys
import os
import csv
import re
import pysam
import sysrsync
import gzip
import glob
import subprocess
from optparse import OptionParser
from pathlib import Path
from tools import *
from dependencies import *

EMBASSY_RNASEQER_PATH = os.environ["EMBASSY_APOLLO_PATH"]

EMBASSY_URL = os.environ["EMBASSY_URL"]
EMBASSY_BUCKET = os.environ["EMBASSY_BUCKET"]
ENA_RUN_API = 'https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=run_accession%3D%22{0}%22&fields=secondary_study_accession%2Cstudy_accession%2Crun_accession%2Csample_accession%2Csample_title%2Csample_alias%2Ccenter_name%2Cstudy_alias%2Cstudy_title&format=json'
JBROWSE_INSTALL_DIR="/hps/software/users/wormbase/parasite/software/jbrowse"
JBROWSE_OUT_DIR="$PARASITE_SCRATCH/jbrowse/WBPS${PARASITE_VERSION}/out"
TEMP_DIR="$PARASITE_SCRATCH/temp"
studies_metadata_tsv = "/nfs/production/flicek/wormbase/parasite/data/releases/release17/schistosoma_mansoni_prjea36577/schistosoma_mansoni_jbrowse_tracks.txt"
out_studies_json = os.environ["PARASITE_SCRATCH"]+"/jbrowse/WBPS"+os.environ["PARASITE_VERSION"]+"/WbpsExpression/schistosoma_mansoni_prjea36577/schistosoma_mansoni.studies.json"

parser = OptionParser(usage='usage: %prog [options] arguments')
parser.add_option("-a", "--apollo_instance_name", dest="apname", help="Name of the apollo instance that will be deployed.")
parser.add_option("-b", "--bam_file", dest="bams", help="Specific bam file on apollo which will be used. If mutliple, use comma-separated list.")

(options, args) = parser.parse_args()

if not options.apname:
    parser.print_help()
    sys.exit(1)

apollo_name = options.apname

if options.bams:
    selected_bams = [x.strip() for x in options.bams.split(",")]
else:
    selected_bams = []

eaws = os.environ["EMBASSY_COMMAND"]
APOLLO_INSTANCE = apollo_name
APOLLO_PATH = EMBASSY_RNASEQER_PATH + "/" + APOLLO_INSTANCE
OUTPUT_DIR = os.path.join(os.environ["PARASITE_SCRATCH"],"rnaseq",APOLLO_INSTANCE,"introns")
TOAPOLLO_DIR = os.path.join(OUTPUT_DIR,"introns_to_apollo")
SAM2INTRONS_software = os.path.join(os.environ["PARASITE_REPOSITORIES"],"perl_ditties","SAM_extract_junctions.pl")
INT_COV = "10"

bams_eaws_run = subprocess.run([eaws+" s3 ls "+APOLLO_PATH+"/ --recursive | awk '{print $4}' | grep '.bam$'"], shell=True, stdout=subprocess.PIPE)
bams_eaws = list(set([x for x in bams_eaws_run.stdout.decode('utf-8').split('\n') if x != '']))

bams = bams_eaws

if selected_bams != []:
    bams = [y for x in selected_bams for y in bams_eaws if x in y]

bams_urls = [os.path.join(EMBASSY_URL,EMBASSY_BUCKET,x) for x in bams]
print(bams_urls)

if not os.path.exists(TOAPOLLO_DIR):
    os.makedirs(TOAPOLLO_DIR)



for bam_url in bams_urls:
    bam_name = os.path.basename(bam_url)
    bam_name_no_ext = ".".join(bam_name.split(".")[0:-1])
    out_dir = os.path.join(OUTPUT_DIR,bam_name_no_ext)
    if os.path.exists(out_dir):
        print_info("Path already exists: "+out_dir+". Duplicate?")
    os.makedirs(out_dir)
    bam_path = os.path.join(out_dir,bam_name)
    bam_path_no_ext = os.path.join(out_dir,bam_name_no_ext)
    bash_command = curl_fastqs({bam_url:bam_path})
    bash_command += "ANC_LEN=$(samtools view "+bam_path+" | " \
                    "awk '{print length($10)}' | " \
                    " head -1000 | sort -u | awk '{s+=$1}END{print s/NR}' | " \
                    "xargs printf '%.*f\\n' \"$p\" | awk '{print $1*0.1}');\n\n"
    bash_command += "if (( $ANC_LEN<3 )); then ANC_LEN=3; fi;\n\n"
    bash_command += "if (( $ANC_LEN>12 )); then ANC_LEN=10; fi;\n\n"
    bash_command += "INT_COV="+INT_COV+";\n\n"
    bash_command += "perl "+SAM2INTRONS_software+" "+bam_path+" $INT_COV $ANC_LEN;\n\n"
    bash_command += "rm "+bam_path+";\n"
    bash_command += "rm " + bam_path + ".${INT_COV}.${ANC_LEN}.sam;\n\n"
    bash_command += "rm -rf " + os.path.join(out_dir,"SAM_extract_junctions")+";\n\n"
    bash_command += "mv " + bam_path + ".${INT_COV}.${ANC_LEN}.gff " + os.path.join(TOAPOLLO_DIR,bam_name_no_ext+".introns.gff") + ";\n\n"
    #print(bash_command)

    lsf_submit(bash_command,
               jobprefix=bam_name_no_ext,
               to_wait_id="",
               cpu=1, mem="2gb", cwd=out_dir, queue="production",
               job_name=bam_name_no_ext,
               only_write=False)




