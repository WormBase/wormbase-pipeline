import os
import sys
import time
import gzip
import datetime
import glob
import subprocess
import re
import random
import sqlalchemy
import string
import csv
import pandas as pd

def dtnow():
    ppdtnow = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return(str(ppdtnow))

def exit_with_error(error_message):
    print(dtnow() + ": ERROR - "+error_message+".")
    raise ValueError

def print_info(info_message):
    print(dtnow() + ": INFO - " + info_message + ".")

def print_warning(info_message):
    print(dtnow() + ": INFO - " + info_message + ".")

def print_w_indent(message):
    print("\t" + message)

def pnl():
    print("\n")

def check_if_file_exists(xfile):
    if os.path.exists(xfile):
        pass
    else:
        exit_with_error("File: "+xfile+" does not exist. Exiting.")

def bash_check_if_success_or_loop(bcommand,times):
    if times<2:
        times=2
    if not bcommand.endswith("\n"):
        bcommand += "\n"
    bash_command = "for run in {1.."+str(times)+"}; do" + "\n" + \
                   "\t"   + bcommand + \
                   "\t"   + 'if [ "$?" = "0" ] ; then' + "\n" + \
                   "\t\t" + 'echo "rsync completed normally";' + "\n" + \
                   "\t\t" + "break;" + "\n" + \
                   "\t"   + "else" + "\n" + \
                   "\t\t" + 'echo "Rsync failed. Retrying...";' + "\n" + \
                   "\t\t" + "sleep 60;" + "\n" + \
                   "\t"   + "fi" + "\n" + \
                   "done;\n"
    return(bash_command)

def decompress(infile, tofile):
    with open(infile, 'rb') as inf, open(tofile, 'w', encoding='utf8') as tof:
        decom_str = gzip.decompress(inf.read()).decode('utf-8')
        tof.write(decom_str)

def md5_check_command(original_file, copied_file):
    bash_command="\n\noriginal_md5=$(md5sum "+original_file+" | cut -d' ' -f1);\n" + \
                 "copied_md5=$(md5sum "+copied_file+" | cut -d' ' -f1);\n" + \
                 'if [ "$original_md5" != "$copied_md5" ];\n' + \
                 "\tthen echo 'File copy was unsuccessful. Exiting!'; exit 1;" + \
                 "fi;\n\n\n"
    return(bash_command)


def coordlist2vcf(coordlist, pgrange_df, col2get="file", chr="chr", start="start", end="end"):
    """Function that takes a list with genomic coordinates: ["chr5","3432423","3243342"]
    (if one of the positions miss you can just replace it with "") and returns a list
    of vcf files corresponding to this coordinates (read from a map text file)."""
    c1=coordlist[0]
    c2=coordlist[1]
    c3=coordlist[2]
    cindex=[]
    if c3=="":
        if c2!="":
            cindex+=pgrange_df.index[(pgrange_df[start]<int(c2)) & (pgrange_df[end]>int(c2)) & (pgrange_df[chr]==str(c1))].tolist()
            if len(cindex)==0:
                print("The coordinates you've entered are either not correct or there are not in the processed vcf files.")
                raise IndexError
            else:
                cfiles = pgrange_df.iloc[cindex][col2get]
            #cfiles=cfiles+pgrange_df[(pgrange_df[start]<int(c3)) & (pgrange_df[end]>int(c3)) & (pgrange_df[chr]==str(c1))][col2get].tolist()
        elif c2=="":
            cindex+=pgrange_df.index[(pgrange_df[chr]==str(c1))].tolist()
            if len(cindex)==0:
                print("The coordinates you've entered are either not correct or there are not in the processed vcf files.")
                raise IndexError
            else:
                cfiles = pgrange_df.iloc[cindex][col2get].tolist()
    elif c2!="" and c3!="":
        cfiles=[]
        cindex+=pgrange_df.index[(pgrange_df[start]<int(c2)) & (pgrange_df[end]>int(c2)) & (pgrange_df[chr]==str(c1))].tolist()
        cindex+=pgrange_df.index[(pgrange_df[start]<int(c3)) & (pgrange_df[end]>int(c3)) & (pgrange_df[chr]==str(c1))].tolist()
        try:
            cfiles+=[pgrange_df.iloc[cindex[0]][col2get]]
        except IndexError:
            print("The coordinates you've entered are not correct or there are not in the processed vcf files.")
        try:
            cfiles+=[pgrange_df.iloc[cindex[1]][col2get]]
        except IndexError:
            print("The coordinates you've entered are not correct or there are not in the processed vcf files.")
        try:
            cfiles+=pgrange_df.iloc[cindex[0]:cindex[0]][col2get].tolist()
        except IndexError:
            print("The coordinates you've entered are not correct or there are not in the processed vcf files.")
    cfiles = list(set(cfiles))
    return(cfiles)

def csvlines2list(csv_path):
    """Function that takes a path for a csv file, it opens it and returns each line as a
    python list element."""
    csv_in = csv_path.strip()
    with open(csv_in, mode='r') as csv_file:
        csv_reader = csv.reader(csv_file)
        outlist=[]
        for row in csv_reader:
            row = [x.strip() for x in row]
            row = [x.split("\xef\xbb\xbf")[1] if x.startswith("\xef\xbb\xbf") else x for x in row]
            outlist+=row
    return(outlist)


def getnum(text):
    retnnum = int(text.split("/")[-1].split("_")[5])
    return(retnnum)


def tidydirs(project, outdir, rundir):
    dirlist=[]
    dirlist.append(outdir+"/"+project) #project_out
    dirlist.append(rundir+"/"+project) #project_run
    dirlist.append(rundir+"/"+project+"/"+"jobscripts") #project_js
    dirlist.append(rundir+"/"+project+"/"+"tmp") #project_tmp
    dirlist.append(outdir+"/"+project+"/genomic_data") #genomic_data
    dirlist.append(outdir+"/"+project+"/interval_data") #interval_data
    for dir in dirlist:
        if not os.path.exists(dir):
            print(dir+" "+"has been created.")
            os.makedirs(dir)
    return(dirlist)


def randomString(stringLength=8):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(stringLength))

def lsf_submit(command, jobprefix, to_wait_id="", cpu=1, mem="4gb", cwd="./", queue="production"):
    """Function which takes a pipeline command, runtime arguments, jobname, job dependencies
    as input, submits the job to the cluster and returns the job id"""
    if (type(to_wait_id) == str):
        to_wait_id = [to_wait_id]
    to_wait_id = [x for x in to_wait_id if x!=""]
    #print("debug:")
    #print(to_wait_id)
    #print(":".join(to_wait_id))
    if cwd.endswith("/") !=True:
        cwd = cwd + "/"
    out_filename = cwd + jobprefix + ".jobscript"
    to_wait_id_to_write = ["done("+x+")" for x in to_wait_id]
    #whours = int(wtime.split(":")[0])
    mem = str(int(mem.split("gb")[0])*1000)
    #if whours < 4: wqueue="short"
    #elif whours >=4 and whours < 24: wqueue="medium"
    #elif whours >= 24: wqueue="long"
    #print(out_filename)
    with open(out_filename, 'w+') as out_file:
        out_file.write('#!/bin/bash')
        out_file.write('\n')
        out_file.write('\n')
        out_file.write('#BSUB -q ' + queue)
        out_file.write('\n')
        #out_file.write('#BSUB -P ' + project_name)
        #out_file.write('\n')
        if len(to_wait_id_to_write)!=0:
            out_file.write('#BSUB -w "'+" && ".join(to_wait_id_to_write)+'"')
            out_file.write('\n')
        #out_file.write('\n')
        out_file.write("#BSUB -n " + str(cpu))
        out_file.write('\n')
        out_file.write("#BSUB -M " + str(mem))
        out_file.write('\n')
        out_file.write("#BSUB -R rusage[mem=" + str(mem) + "]")
        out_file.write('\n')
        out_file.write("#BSUB -oo " + str(cwd + "/" + jobprefix + ".stdout"))
        out_file.write('\n')
        out_file.write("#BSUB -outdir " + cwd)
        out_file.write('\n')
        out_file.write("#BSUB -eo " + str(cwd + "/" + jobprefix + ".stderr"))
        out_file.write('\n')
        out_file.write("#BSUB -cwd " + str(cwd))
        out_file.write('\n')
        out_file.write('\n')
        out_file.write(command)
        out_file.write('\n')
        out_file.write('sleep 2')
        out_file.write('\n')
        out_file.write('exit')
    out_file.close()
    time.sleep(0.1)
    p = subprocess.Popen(["bsub"],stdout=subprocess.PIPE, stderr=subprocess.PIPE,stdin=open(out_filename, 'r'))
    out, err = p.communicate()
    jout = out
    jerr = err
    mout = re.findall(r"<(\d+)>",jout.decode("utf-8"))
    if mout[0]:
        mout = mout[0]
    else:
        print("No job identifier detected")
        raise()
    if "Job not submitted" in jerr.decode("utf-8"):
        print(jerr.decode("utf-8"))
        raise()
    return mout

def flatten(t):
    return [item for sublist in t for item in sublist]

def getnum(text):
    retnnum = int(text.split("/")[-1].split("_")[5])
    return(retnnum)

def rsync_command(origin, destination):
    bash_command="rsync -avz --checksum " + \
                 origin + " " + \
                 destination + ";\n"
    return(bash_command)
def cram2bam_command(cram, bam, fasta):
    bash_command="samtools view -b " + \
                 "-T " + fasta + " " + \
                 "-o " + bam + " " + \
                 cram + ";\n\n"
    return (bash_command)
def namesort_bam_command(inbam, outbam):
    bash_command="samtools sort -n " + \
                 "-o " + outbam + " " + \
                 inbam + ";\n\n"
    return bash_command
def positionsort_bam_command(inbam, outbam):
    bash_command="samtools sort " + \
                 "-o " + outbam + " " + \
                 inbam + ";\n\n"
    return bash_command
def fixmate_bam_command(inbam, outbam):
    bash_command="samtools fixmate -m " + \
                 inbam + " " + \
                 outbam + ";\n\n"
    return bash_command
def markdup_bam_command(inbam, outbam):
    bash_command="samtools markdup -r " + \
                 inbam + " " + \
                 outbam + ";\n\n"
    return bash_command
def index_bam_command(inbam):
    bash_command="samtools index " + \
                 inbam + ";\n\n"
    return bash_command
def merge_bams_command(inbams, outbam, threads):
    if (type(inbams) == str):
        inbams = [inbams]
    bash_command="samtools merge -f " + \
        "--threads " + str(threads) + " " + \
        outbam + " " + \
        " ".join(inbams) + ";\n\n"
    return(bash_command)
def sortsamrefname_command(inbam, outbam, sortsamrefname):
    bash_command=sortsamrefname + " " + \
                 "--samoutputformat BAM " + \
                 "-o " + outbam + " " + \
                 inbam + ";\n\n"
    return(bash_command)
def capbam_command(inbam, outbam, cap_reads, biostar154220):
    bash_command=biostar154220 + " " + \
                 "-n " + str(cap_reads) + "  " + \
                 "--samoutputformat BAM " + \
                 "-o " + outbam + " " + \
                 inbam + ";\n\n"
    return bash_command
def delete_if_ok(todelete):
    bash_command="if [ $? -eq 0 ]; then\n" + \
    '\trm ' + todelete + '\n' + \
    'else\n' + \
    '\techo 1>&2 "command failed";\n' \
    '\texit 1;\n' + \
    'fi;\n\n'
    return bash_command

def regex_match_dbs(pattern, databases):
    r = re.compile(".*" + pattern + ".*")
    filtdb_list = list(filter(r.match, databases))
    return (filtdb_list)

def regex_match_one_db(pattern, databases):
    r = re.compile(".*" + pattern + ".*")
    filtdb_list = list(filter(r.match, databases))
    if len(filtdb_list) == 0:
        exit_with_error("no db for: " + pattern)
    elif len(filtdb_list) > 1:
        exit_with_error("multiple_dbs_for: " + pattern)
    else:
        return (filtdb_list[0])