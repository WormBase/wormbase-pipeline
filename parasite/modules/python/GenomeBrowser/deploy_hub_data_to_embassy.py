import os
from LSF import *
from ProductionUtils import *

EMBASSY_COMMAND = os.environ["EMBASSY_COMMAND"]
PARASITE_SCRATCH = os.environ["PARASITE_SCRATCH"]
PARASITE_VERSION = os.environ["PARASITE_VERSION"]
EMBASSY_PATH = os.environ["EMBASSY_PATH"]

job_title = "wdjdwbps"+PARASITE_VERSION
parallel_jobs_limit = 10
jbrowse_species_dir = os.path.join(PARASITE_SCRATCH,"jbrowse","WBPS"+PARASITE_VERSION,"hub")
jbrowse_species = [x for x in os.listdir(jbrowse_species_dir) if os.path.isdir(jbrowse_species_dir+"/"+x)]
jbrowse_log_dir = os.path.join(PARASITE_SCRATCH,"log","jbrowse_webdata_hub_deploy")

if not os.path.exists(jbrowse_log_dir):
    os.makedirs(jbrowse_log_dir)

for js in jbrowse_species:
    if js=="Schistosoma_mansoni_prjea36577":
        continue
    dcommand = "module load embassy;\n\n" + \
               EMBASSY_COMMAND + " s3 sync " + os.path.join(PARASITE_SCRATCH,"jbrowse","WBPS"+PARASITE_VERSION,"hub",js) + " " +\
               os.path.join(EMBASSY_PATH,"web_data/track_hubs/rnaseq/releases/release-"+PARASITE_VERSION,js) +\
               " --cli-connect-timeout 0;"
    running_jobs = len([x for x in lsf_submitted_jobnames() if job_title in x])

    if running_jobs < parallel_jobs_limit:
        print_info("Submitting job for "+js)
        lsf_submit(command=dcommand,
                   jobprefix="do_jbrowse_deploy." + js,
                   to_wait_id="",
                   cpu=1,
                   mem="1gb",
                   cwd=jbrowse_log_dir,
                   queue="production", job_name=job_title + "_" + js,
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
                   queue="production", job_name=job_title + "_" + js,
                   only_write=False)
















