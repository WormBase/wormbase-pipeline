import subprocess
import time
import re
from ProductionUtils import print_info

def lsf_submit(command, jobprefix, to_wait_id="", cpu=1, mem="4gb", cwd="./", queue="production", job_name=None, only_write=False):
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
        if job_name!=None:
            out_file.write('#BSUB -J ' + job_name)
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
        out_file.write("#BSUB -oo " + str(cwd + "/" + jobprefix + ".%J.stdout"))
        out_file.write('\n')
        out_file.write("#BSUB -outdir " + cwd)
        out_file.write('\n')
        out_file.write("#BSUB -eo " + str(cwd + "/" + jobprefix + ".%J.stderr"))
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
    if only_write:
        return out_filename
    else:
        p = subprocess.Popen(["bsub"],stdout=subprocess.PIPE, stderr=subprocess.PIPE,stdin=open(out_filename, 'r'))
        out, err = p.communicate()
        jout = out
        jerr = err
        mout = re.findall(r"<(\d+)>" , jout.decode("utf-8"))
        if mout[0]:
            mout = mout[0]
        else:
            print("No job identifier detected")
            raise()
        if "Job not submitted" in jerr.decode("utf-8"):
            print(jerr.decode("utf-8"))
            raise()
        return mout

def lsf_job_status(job_id):
    p = subprocess.Popen(["bjobs", "-l", job_id], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    jout = out
    jerr = err
    mout = re.findall(r", Status <(\w+)>," , jout.decode("utf-8"))
    try:
        mout = mout[0]
    except IndexError:
        mout = "DONE"
    if "Job not submitted" in jerr.decode("utf-8"):
        mout = "NOT SUBMITTED"
    return mout

def lsf_bjobs():
    p = subprocess.Popen(["bjobs", "-w"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    jout = out
    jerr = err
    return(jout.decode("utf-8"))

def lsf_submitted_jobnames():
    lsf_split = [x.split()[6] for x in lsf_bjobs().split("\n") if x != '']
    return([x for x in lsf_split if x!="JOB_NAME"])

def lsf_running_jobnames():
    lsf_split = [x.split()[6] for x in lsf_bjobs().split("\n") if x != '' and x.split()[2]=="RUN"]
    return ([x for x in lsf_split if x != "JOB_NAME"])

def batch_jobs_check_status_periodically(batch_jobs_dict):
    """Expects a dictionary as its input. Keys are references (i.e. a genome) and values are job_ids."""
    all_jobs = list(batch_jobs_dict.values())
    while len(all_jobs)>0:
        for ref in batch_jobs_dict:
            job_id = batch_jobs_dict[ref]
            if lsf_job_status(job_id)!="RUN" and lsf_job_status(job_id)!="PEND":
                if job_id in all_jobs: all_jobs.remove(job_id)
                # print_info(job_id +" is DONE")
                if lsf_job_status(job_id) == "EXIT":
                    print_info(ref+" - "+job_id + " has FAILED.")
        print_info("Number of jobs still running: " + str(len(all_jobs)))
        time.sleep(15)