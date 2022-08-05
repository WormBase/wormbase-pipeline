import os
from ProductionMysql import *
from LSF import *
from optparse import OptionParser
import subprocess
import re

BRC4PIPELINE_DIR = os.environ["RNASEQ_SCRATCH"]
LIVE_WEBDATA_DIR = os.path.join(os.environ["PARASITE_FTP"], "web_data")
LIVE_JBROWSE_DIR = os.path.join(LIVE_WEBDATA_DIR, "jbrowse", "releases", "release-"+os.environ["PARASITE_VERSION"])
EMBASSY_RNASEQER_PATH = os.environ["EMBASSY_RNASEQER_PATH"]
EMBASSY_PATH = os.environ["EMBASSY_PATH"]
EMBASSY_BUCKET = os.environ["EMBASSY_BUCKET"]
eaws = os.environ["EMBASSY_COMMAND"]

parser = OptionParser(usage='usage: %prog [options] arguments')
parser.add_option("-g", "--genome", dest="GENOME",
                  help="Optional: Genome for which .bw files will be deployed (in the spe_cies_bp format). Example "
                       "schistosoma_mansoni_prjea36577. For multiple entries use a comma-separated list.")
(options, args) = parser.parse_args()

inputs = {x:globals()[x] for x in globals() if x.endswith("_DIR") or x.endswith("_PATH")
          or x.endswith("aws") or x.endswith("_BUCKET") }

top_level_dirs = next(os.walk(BRC4PIPELINE_DIR))[1]
genomes = [x for x in top_level_dirs if len(staging.core_dbs(x+"_core"))==1]

if options.GENOME:
    GENOMES = [x for x in options.GENOME.strip().split(",") if x in genomes]
    if len(GENOMES)>=1:
        print_info("User defined genomes: "+", ".join(GENOMES))
    else:
        exit_with_info("Wrong genome(s) defined. please try again")
else:
    GENOMES = genomes
    print_info("No genomes were defined by the user, all genomes be used instead:\n"+"\n".join(GENOMES))

class Genome:
    def __init__(self, genome, inputs):
        self.genome_dir = os.path.join(inputs["BRC4PIPELINE_DIR"], genome)
        self.studies_dir = os.path.join(self.genome_dir, "component", genome)
        self.studies = next(os.walk(self.studies_dir))[1]
        #self.live_jbrowse_dir = os.path.join(inputs["LIVE_JBROWSE_DIR"], genome)
        #self.live_trackList = os.path.join(self.live_jbrowse_dir, "data", "trackList.json")
        self.genome_embassy_s3_path = os.path.join(inputs["EMBASSY_RNASEQER_PATH"], genome)
        self.genome_embassy_path = self.genome_embassy_s3_path.replace(inputs["EMBASSY_PATH"]+"/","")

class Study(Genome):
    def __init__(self, genome, inputs, study_id):
        super().__init__(genome, inputs)
        self.id = study_id
        self.study_dir = os.path.join(self.studies_dir, self.id)
        self.to_embassy_dir = os.path.join(self.study_dir, "to_embassy")
        self.runs = [x for x in next(os.walk(self.study_dir))[1] if x != "to_embassy"]
    def create_dirs(self):
        run_once = 0
        for dir in self.__dict__.keys():
            pot_path = self.__dict__[dir]
            if dir.endswith("dir") and not check_dir_exists(pot_path):
                if run_once == 0:
                    print_info("Creating Required paths:")
                print_w_indent("Creating " + pot_path)
                os.makedirs(pot_path)
                run_once += 1

class Run(Study):
    def __init__(self, genome, inputs, study_id, run_id):
        super().__init__(genome, inputs, study_id)
        self.inputs = inputs
        self.id = run_id
        self.study_id = study_id
        self.run_dir = os.path.join(self.study_dir, self.id)
        self.run_bw_original_name = os.path.join(self.run_dir, "results.bw")
        self.run_bw = os.path.join(self.to_embassy_dir, self.id+".bw")
        self.run_bw_filename = self.id+".bw"

    def rename_and_move_bw_command(self, to_wait_id=""):
        bash_command = "cp " + self.run_bw_original_name + " " + self.run_bw + ";\n\n"
        return(bash_command)

    def is_bw_deployed_in_embassy(self):
        command = os.environ["EMBASSY_COMMAND"] + ' s3api list-objects ' + \
            '--bucket ' + self.inputs["EMBASSY_BUCKET"] + ' ' + \
            '--prefix ' + self.genome_embassy_path + ' ' + \
            '--query "Contents[?contains(Key, \''+self.id+'.bw\')]"'
        result = subprocess.run(command, stdout=subprocess.PIPE, shell=True)
        if len(str(result.stdout).replace("\\n", "").split("},")) == 1:
            entry = str(result.stdout).replace("\\n", "").split("},")[0]
            remote_bw_path = re.search('\{\s+"Key":\s+"((\S+\/){1,}\w+.bw)",\s+', entry)[1]
            return(remote_bw_path.strip())
        elif len(str(result.stdout).replace("\\n", "").split("},")) > 1:
            return("Multi")
        else:
            return(False)

    def embassy_target_s3_path(self):
        target_bigwig_pre = self.is_bw_deployed_in_embassy()
        if target_bigwig_pre == False:
            target_bigwig_s3_path = "Couldn't find bigwig in embassy for " + run_id + ". Skipping..."
            raise FileNotFoundError
        elif target_bigwig_pre == "Multi":
            target_bigwig_s3_path = "Found multiple bigwigs in embassy for " + run_id + ". Skipping..."
            raise OSError
        else:
            target_bigwig_s3_path = EMBASSY_PATH + "/" + target_bigwig_pre
        return(target_bigwig_s3_path)

    def deploy_bw_command(self):
        bash_command = os.environ["EMBASSY_COMMAND"] + ' s3 cp ' + \
            self.run_bw + " " + self.embassy_target_s3_path() + ";\n\n"
        return(bash_command)


for genome_name in GENOMES:
    genome = Genome(genome_name, inputs)
    studies = genome.studies
    print_info("Processing: "+genome_name)
    for study_id in studies:
        if study_id == "SRP122521":
            continue
        study = Study(genome_name, inputs, study_id)
        study.create_dirs()
        runs = study.runs
        print_info("Study: " + study_id)
        for run_id in runs:
            run = Run(genome_name, inputs, study_id, run_id)
            print_info("Run: " + run_id)
            try:
                bash_command = run.rename_and_move_bw_command() + \
                               run.deploy_bw_command()
            except FileNotFoundError:
                print("Couldn't find bigwig in embassy for " + run_id + ". Skipping...")
                continue
            except OSError:
                print("Found multiple bigwigs in embassy for " + run_id + ". Skipping...")
                continue
            job_id = lsf_submit(bash_command,
                                jobprefix=run_id+"_bigwig_rename",
                                to_wait_id="",
                                cpu=1, mem="1gb",
                                cwd=run.run_dir,
                                queue="production",
                                only_write=False)



