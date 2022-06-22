import os
from LSF import *
from ProductionUtils import *

PARASITE_SCRATCH = os.environ["PARASITE_SCRATCH"]
PARASITE_VERSION = os.environ["PARASITE_VERSION"]
PARASITE_FTP = os.environ["PARASITE_FTP"]

job_title = "wdjdwbps_ftp"+PARASITE_VERSION
parallel_jobs_limit = 10
jbrowse_species_dir = os.path.join(PARASITE_SCRATCH,"jbrowse","WBPS"+PARASITE_VERSION,"out")
jbrowse_species = [x for x in os.listdir(jbrowse_species_dir) if os.path.isdir(jbrowse_species_dir+"/"+x)]
jbrowse_log_dir = os.path.join(PARASITE_SCRATCH,"log","jbrowse_webdata_deploy_to_ftp")

if not os.path.exists(jbrowse_log_dir):
    os.makedirs(jbrowse_log_dir)

for js in jbrowse_species:
    dcommand = "sudo -u wormbase mkdir -p "+ \
               os.path.join(PARASITE_FTP,"web_data/jbrowse/releases","release-"+PARASITE_VERSION,js,"data") + " " + \
               "&& sudo -u wormbase rsync -avz " + \
               os.path.join(PARASITE_SCRATCH,"jbrowse","WBPS"+PARASITE_VERSION,"out",js+"/") + " " + \
               os.path.join(PARASITE_FTP,"web_data/jbrowse/releases","release-"+PARASITE_VERSION,js,"data/")
    running_jobs = len([x for x in lsf_submitted_jobnames() if job_title in x])

    if running_jobs < parallel_jobs_limit:
        print_info("Submitting job for "+js)
        lsf_submit(command=dcommand,
                   jobprefix="do_jbrowse_deploy." + js,
                   to_wait_id="",
                   cpu=1,
                   mem="1gb",
                   cwd=jbrowse_log_dir,
                   queue="datamover", job_name=job_title + "_" + js,
                   only_write=False)
    else:
        while running_jobs >= parallel_jobs_limit:
            print_info("Reached jobs limit. Checking again in a minute...")
            time.sleep(60)
            running_jobs = len([x for x in lsf_submitted_jobnames() if job_title in x])
        print_info("Submitting job for " + js)
        lsf_submit(command=dcommand,
                   jobprefix="do_jbrowse_deploy." + js,
                   to_wait_id="",
                   cpu=1,
                   mem="1gb",
                   cwd=jbrowse_log_dir,
                   queue="datamover", job_name=job_title + "_" + js,
                   only_write=False)















