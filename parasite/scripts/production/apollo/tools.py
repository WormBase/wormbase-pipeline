import os
import pysam
import glob
import datetime
import sysrsync
import gzip
from pathlib import Path
import dependencies
from dependencies import *
from diotools import *

def check_dependent_dirs():
    dep_dirs_names=[item for item in dir(dependencies) if item.endswith('dir')]
    for dir_name in dep_dirs_names:
        dir_path = globals()[dir_name]
        if not os.path.exists(dir_path):
            exit_with_error("The required path: "+dir_path+" does not exist. Please check"
                                                           "what is in the dependencies.py"
                                                           "script and try again.")

class Species:
    def __init__(self, species_name):
        self.species = species_name
        self.apollo_dir = apollo_dir+"/"+species_name
        self.permanent_apollo_dir = permanent_apollo_dir+"/"+species_name
        self.bam_dir = apollo_dir+"/"+species_name+"/bam"
        self.cram_dir = apollo_dir+"/"+species_name+"/cram"
        self.bigwig_dir = apollo_dir+"/"+species_name+"/bigwig"
        self.se_dir = expression_dir+"/"+species_name
        self.reference_dir = os.path.join(self.apollo_dir,"reference")
        self.log_dir = os.path.join(self.apollo_dir,"log")
    def create_dirs(self):
        run_once=0
        for dir in self.__dict__.keys():
            pot_path=self.__dict__[dir]
            if dir.endswith("dir") and not os.path.exists(pot_path):
                if run_once==0:
                    print_info("Creating Required paths:")
                print_w_indent("Creating "+pot_path)
                os.makedirs(pot_path)
                run_once+=1
    def report_input(self):
        dep_dirs_names = [item for item in dir(dependencies) if item.endswith('dir')]
        print_info("Given Input (from user and dependencies.py):")
        print_w_indent("Species: "+self.species)
        for depdir in self.__dict__.keys():
            pot_path=self.__dict__[depdir]
            if depdir.endswith("dir"):
                print_w_indent(depdir + ": " + pot_path)
    def studies(self):
        study_ids = [x for x in next(os.walk(self.se_dir))[1]
                      if os.path.isfile(self.se_dir + "/" + x + "/" + x + ".design.tsv")]
        return study_ids


class Study(Species):
    def __init__(self, species_name, study_id):
        super().__init__(species_name)
        self.study_id = study_id
        self.spath = os.path.join(self.se_dir, self.study_id)
        self.design_file = os.path.join(self.spath, self.study_id + ".design.tsv")
        self.study_bam_dir = os.path.join(self.bam_dir,self.study_id)
        self.study_bigwig_dir = os.path.join(self.bigwig_dir,self.study_id)
        self.ena_api_url = ena_api_url.format(self.study_id)
        self.genome_fa_prepared = 0
        self.ref_fasta = os.path.join(self.reference_dir, self.species+"."+self.study_id+".genomic.fa")
        self.merged_dir = os.path.join(self.bam_dir,"merged")
        self.merged_bam = os.path.join(self.merged_dir, self.study_id+".merged.bam")
        self.merged_bam_noext = str(Path(self.merged_bam).with_suffix(''))
        self.merged_bam_sortrefnamed = self.merged_bam_noext + ".sortrefnamed.bam"
        self.merged_bam_capped = self.merged_bam_noext + ".TEMP.capped.bam"
        self.merged_bam_final = self.merged_bam_noext + ".capped.bam"
    def design_dict(self):
        design_reader = csv.DictReader(open(self.design_file), delimiter="\t")
        desdict = {self.study_id + '_' + re.sub(r'[^\w]', '', row['Condition']): row['Run'] for row in design_reader}
        return desdict
    def design_keys(self):
        return list(self.design_dict().keys())
    def samples(self):
        return list(self.design_dict.values())
    def create_study_dirs(self):
        run_once = 0
        for dir in self.__dict__.keys():
            pot_path = self.__dict__[dir]
            if dir.endswith("dir") and not os.path.exists(pot_path):
                if run_once == 0:
                    print_info("Creating Required paths:")
                print_w_indent("Creating " + pot_path)
                os.makedirs(pot_path)
                run_once += 1
    def ref_fasta_ftpfile(self):
        for dek in self.design_keys():
            ssa = Sample(self.species, self.study_id, dek)
            fasta = ssa.ref_fasta_ftpfile_sample()
            break
        return fasta
    def is_fasta_gzipped(self):
        return self.ref_fasta_ftpfile().endswith(".gz")
    def ref_fasta_copied(self):
        if self.is_fasta_gzipped():
            return self.ref_fasta+".gz"
        else:
            return self.ref_fasta
    def copy_ref_fasta(self):
        if not os.path.exists(self.ref_fasta):
            print_info("Copying reference fasta file "+self.ref_fasta_ftpfile()+" for study "+self.study_id+".")
            sysrsync.run(
                source=self.ref_fasta_ftpfile(),
                destination=self.ref_fasta_copied(),
                options=['-avz', '--checksum'])
            if self.is_fasta_gzipped():
                self.decompress_ref_fasta()
    def decompress_ref_fasta(self):
        decompress(self.ref_fasta_copied(), self.ref_fasta)
        if os.path.isfile(self.ref_fasta):
            os.remove(self.ref_fasta_copied())
    def bigwigs(self):
        bigwigs = []
        for dek in self.design_keys():
            ssa = Sample(self.species, self.study_id, dek)
            bigwig = ssa.bigwig
            bigwigs.append(bigwig)
        return bigwigs
    def bams(self):
        bams = []
        for dek in self.design_keys():
            ssa = Sample(self.species, self.study_id, dek)
            bam = ssa.bam_deduped
            bams.append(bam)
        return bams
    def merge_and_cap_bams_command(self):
        bash_command = merge_bams_command(self.bams(), self.merged_bam, 4) + \
                       index_bam_command(self.merged_bam) + \
                       sortsamrefname_command(self.merged_bam, self.merged_bam_sortrefnamed, sortsamrefname) + \
                       delete_if_ok(self.merged_bam) + \
                       capbam_command(self.merged_bam_sortrefnamed, self.merged_bam_capped, 50, biostar154220) + \
                       delete_if_ok(self.merged_bam_sortrefnamed) + \
                       positionsort_bam_command(self.merged_bam_capped, self.merged_bam_final) + \
                       delete_if_ok(self.merged_bam_capped) + \
                       index_bam_command(self.merged_bam_final)
        return bash_command
    def merge_and_cap_bams_command_submit(self, to_wait_id=""):
        job_id = lsf_submit(self.merge_and_cap_bams_command(),
                            jobprefix=self.study_id+"_01_merge_bams",
                            to_wait_id=to_wait_id,
                            cpu=4,
                            mem="12gb",
                            cwd=self.log_dir,
                            queue="production")
        return(job_id)




class Sample(Study):
    def __init__(self, species_name, study_id, design_key):
        super().__init__(species_name, study_id)
        self.sample_name = design_key
        self.sample_id = self.design_dict()[self.sample_name]
        self.shortid = self.sample_id[0:6]
        self.ftp = os.path.join(rnaseq_ftp_dir,self.shortid,self.sample_id)
        self.sample_bam_dir = os.path.join(self.study_bam_dir,self.sample_id)
        self.sample_bigwig_dir = os.path.join(self.study_bigwig_dir, self.sample_id)
        self.bam = os.path.join(self.sample_bam_dir, self.sample_name+".bam")
        self.bam_noext = str(Path(self.bam).with_suffix(''))
        self.bam_nsorted = self.bam_noext + ".nsorted.bam"
        self.bam_psorted = self.bam_noext + ".psorted.bam"
        self.bam_fixmated = self.bam_noext + ".fixmated.bam"
        self.bam_deduped = self.bam_noext + ".deduped.bam"
        self.bigwig = os.path.join(self.sample_bigwig_dir, self.sample_name+".bigwig")
        self.fasta = self.ref_fasta
    def create_sample_dirs(self):
        for dir_name in [self.sample_bam_dir, self.sample_bigwig_dir]:
            if not os.path.exists(dir_name):
                print_info("Creating Required paths:")
                print_w_indent("Creating "+dir_name)
                os.makedirs(dir_name)
    def bam_ftpfile(self):
        bamfiles = glob.glob(os.path.join(self.ftp,self.sample_id+".bam"))
        cramfiles = glob.glob(os.path.join(self.ftp,self.sample_id+".cram"))
        if len(bamfiles)>0:
            alfile = bamfiles[0]
        elif len(bamfiles)<1 and len(cramfiles)>0:
            alfile = cramfiles[0]
        else:
            exit_with_error("Could not find BAM/CRAM file in FTP"
                            "directory: " + self.ftp + " for " + self.sample_id)
        return(alfile)
    def bigwig_ftpfile(self):
        bigwigfiles = glob.glob(os.path.join(self.ftp,self.sample_id+".nospliced.bw"))
        if len(bigwigfiles)<0:
            exit_with_error("Could not find a bigwig file in FTP"
                            "directory: " + self.ftp + " for " + self.sample_id)
        return(bigwigfiles[0])
    def bam_copied(self):
        if self.bam_ftpfile().endswith(".cram"):
            return str(Path(self.bam).with_suffix('.cram'))
        else:
            return self.bam
    def is_cram(self):
        return self.bam_ftpfile().endswith(".cram")
    def bioproject(self):
        pyal = pysam.AlignmentFile(self.bam_ftpfile(), "r")
        fasta = list(set([x['UR'] for x in pyal.header.to_dict()['SQ']]))[0]
        fasta_basename = os.path.basename(fasta)
        bioproject = fasta_basename.split(".")[1].strip()
        return(bioproject)
    def ref_fasta_ftpfile_sample(self):
        fasta = glob.glob(os.path.join(ps_ftp_dir,self.species,self.bioproject().upper(),"*.genomic.fa.*"))[0]
        check_if_file_exists(fasta)
        #print(os.path.join(ps_ftp_dir,self.species,self.bioproject().upper(),"*.genomic.fa.*"))
        return(fasta)
    def copy_bam_from_ftp_command(self):
        bash_command=rsync_command(self.bam_ftpfile(), self.bam_copied())
        return(bash_command)
    def copy_bigwig_from_ftp_command(self):
        bash_command = rsync_command(self.bigwig_ftpfile(), self.bigwig)
        return(bash_command)
    def copy_command(self):
        bash_command = bash_check_if_success_or_loop(self.copy_bam_from_ftp_command(),times=3) + \
                       bash_check_if_success_or_loop(self.copy_bigwig_from_ftp_command(),times=3)
        return(bash_command)
    def bam_processing_command(self):
        bash_command = "" + \
                       (cram2bam_command(self.bam_copied(), self.bam, self.fasta) +
                            delete_if_ok(self.bam_copied()) if self.is_cram() else "") + \
                       namesort_bam_command(self.bam, self.bam_nsorted) + \
                       delete_if_ok(self.bam) + \
                       fixmate_bam_command(self.bam_nsorted, self.bam_fixmated) + \
                       delete_if_ok(self.bam_nsorted) + \
                       positionsort_bam_command(self.bam_fixmated, self.bam_psorted) + \
                       delete_if_ok(self.bam_fixmated) + \
                       markdup_bam_command(self.bam_psorted, self.bam_deduped) + \
                       delete_if_ok(self.bam_psorted) + \
                       index_bam_command(self.bam_deduped)
        return bash_command
    def copy_command_submit(self, to_wait_id=""):
        job_id = lsf_submit(self.copy_command(),
                            jobprefix=self.sample_name+"_01_bam_copy",
                            to_wait_id=to_wait_id,
                            cpu=1,
                            mem="1gb",
                            cwd=self.log_dir,
                            queue="datamover")
        return(job_id)
    def bam_processing_command_submit(self, to_wait_id=""):
        job_id = lsf_submit(self.bam_processing_command(),
                            jobprefix=self.sample_name+"_02_bam_processing",
                            to_wait_id=to_wait_id,
                            cpu=1,
                            mem="4gb",
                            cwd=self.log_dir,
                            queue="production")
        return(job_id)










