#!/hps/software/users/wormbase/parasite/shared/.pyenv/versions/p395/bin/python

import os
import pysam
import glob
import time
import datetime
import sysrsync
import gzip
from optparse import OptionParser
from pathlib import Path
import dependencies
from dependencies import *
from ProductionMysql import *
from bash_commands import *
from LSF import *
from Deployment import *
from Alignments.bash_commands import *
from Reads.bash_commands import *
from Utils import *


class ui:
    def __init__(self, parser_object):
        (options, args) = parser_object.parse_args()
        if options.COREDB:
            self.core_db = options.COREDB
        else:
            parser.error('core_db not given')
        self.explore = options.explore
        self.apdep = options.apdep
        self.apname = options.apname
        self.cpext = options.cpext
        if options.genome_ref_dir:
            check_dir_exists(options.genome_ref_dir, to_raise=True)
            self.genome_ref_dir = options.genome_ref_dir
        else:
            if check_file_exists(os.path.join(reference_genomes_dir, self.core_db, "Genome")):
                self.genome_ref_dir = os.path.join(reference_genomes_dir, self.core_db)
            else:
                exit_with_error(
                    "No STAR aligner genome reference has been found. Please create one (check here: ) in " + \
                    os.path.join(reference_genomes_dir, self.core_db) + " " + \
                    "and then re-run")
        if options.genome_gtf:
            check_dir_exists(options.genome_ref_dir, to_raise=True)
            self.genome_gtf = options.genome_gtf
        else:
            candidate_reference_dir = os.path.join(reference_genomes_dir, self.core_db)
            if check_dir_exists(candidate_reference_dir):
                candidate_gtfs = glob.glob(os.path.join(candidate_reference_dir, "*.gtf"))
                if len(candidate_gtfs) > 0:
                    self.genome_gtf = candidate_gtfs[0]
                else:
                    exit_with_error(
                        "No gene models gtf for the genome reference has been detected. Please place one in " + \
                        os.path.join(reference_genomes_dir, self.core_db) + " " + \
                        "and then re-run")

        if options.selects:
            self.selects = [x.strip() for x in options.selects.split(",")]
            if len(self.selects) < 1:
                exit_with_error("Wrong select_studies (-t) input given. Please review and try again.")
            for sel in self.selects:
                if not validate_selects_format(sel):
                    exit_with_error(sel + " is not a valid format for -t. Please try again.")
        if self.apdep:
            if not options.apname: exit_with_error(
                "The apollo name instance has not been given. Please use the -a flag and provide one.")
            self.apname = options.apname
        if self.explore or options.selects == []:
            self.selected_studies = []

    def print_user_input(self):
        print_info("User Specified input:")
        for uo in self.__dict__.keys():
            actual_value = self.__dict__[uo]
            print_w_indent(uo + ": " + str(actual_value))


def check_dependent_dirs():
    dep_dirs_names = [item for item in dir(dependencies) if item.endswith('dir')]
    for dir_name in dep_dirs_names:
        dir_path = globals()[dir_name]
        check_dir_exists(dir_path, to_raise=True)


def check_dependent_software():
    dep_sw_paths = [item for item in dir(dependencies) if item.startwith('sw') and item.endswith('path')]
    for sw_path in dep_sw_paths:
        dir_path = globals()[dir_name]
        check_dir_exists(dir_path, to_raise=True)


class Species:
    def __init__(self, ui, studies_dict=None):
        self.ui = ui
        self.cdb = self.ui.core_db
        self.species = "_".join(self.cdb.split("_")[0:2])
        self.apname = self.ui.apname
        self.species_dir = os.path.join(rnaseq_dir, self.cdb)
        self.studies_dir = os.path.join(self.species_dir, "studies")
        self.to_move_dir = os.path.join(self.species_dir, "to_move")
        self.taxon_id = Core(staging.host, self.cdb).meta_value("species.taxonomy_id")
        self.species_ftp_path = Core(staging.host, self.cdb).ftp_filename_n_filename()
        self.species_ftp_genome = self.species_ftp_path + ".genomic.fa.gz"
        self.species_ftp_gtf = self.species_ftp_path + ".canonical_geneset.gtf.gz"
        self.get_rna_seq_studies_by_taxon_api_url = ena_rnaseq_by_taxon_url.format(self.taxon_id)
        self.aligner_reference = self.ui.genome_ref_dir
        self.aligner_gtf = self.ui.genome_gtf
        if self.ui.selects:
            self.selects = self.ui.selects
        if studies_dict:
            self.studies_dict = studies_dict
            self.studies = list(self.studies_dict.keys())
            self.samples = flatten([[y["run_accession"] for y in self.studies_dict[x]] for x in self.studies])
            self.studies_to_samples_dict = {x: [y["run_accession"] for y in self.studies_dict[x]] for x in self.studies}

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




class Study(Species):

    def __init__(self, ui, studies_dict, study_id):
        super().__init__(ui, studies_dict)
        self.study_id = study_id
        self.study_dir = os.path.join(self.studies_dir, self.study_id)
        self.merged_dir = os.path.join(self.study_dir, "merged")
        self.to_apollo_dir = os.path.join(self.study_dir, "to_apollo")
        self.log_dir = os.path.join(self.study_dir, "log")

        self.dict = {self.study_id: self.studies_dict[self.study_id]}
        self.sample_ids = [x["run_accession"] for x in self.dict[self.study_id]]
        self.samples_dict = {x["run_accession"]: x for x in self.dict[self.study_id]}

        self.merged_bam = os.path.join(self.merged_dir, self.study_id + ".merged.bam")
        self.merged_bam_noext = str(Path(self.merged_bam).with_suffix(''))
        self.merged_bam_sortrefnamed = os.path.join(self.merged_bam_noext, ".sortrefnamed.bam")
        self.merged_bam_capped = os.path.join(self.merged_bam_noext, ".TEMP.capped.bam")
        self.merged_bam_final = os.path.join(self.merged_bam_noext, ".capped.bam")

        self.embassy_apollo_rsync_path = os.path.join(embassy_apollo_path, self.apname, self.study_id)

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

    def bigwigs(self):
        bigwigs = []
        for sample_id in self.sample_ids:
            ssa = Sample(self.ui, self.studies_dict, self.study_id, sample_id)
            bigwig = ssa.bigwig
            bigwigs.append(bigwig)
        return bigwigs

    def bams(self):
        bigwigs = []
        for sample_id in self.sample_ids:
            ssa = Sample(self.ui, self.studies_dict, self.study_id, sample_id)
            bam = ssa.bam_deduped
            bams.append(bam)
        return bams

    def files_to_apollo(self):
        files_to_apollo = []
        for sample_id in self.sample_ids:
            ssa = Sample(self.ui, self.studies_dict, self.study_id, sample_id)
            files_to_apollo.append(ssa.bam_for_apollo)
            files_to_apollo.append(ssa.bigwig_for_apollo)
        return files_to_apollo

    def move_all_finals_to_final_dir_command(self):
        bash_commands = [do_if_exist(x, "\tmv " + \
                                     x + "* " + \
                                     self.to_apollo_dir) for x in self.files_to_apollo()]
        bash_command = "\n".join(bash_commands)
        return (bash_command)

    def deploy_to_embassy_command(self):
        bash_command = rsync_to_external_location_command(directory_from=self.to_apollo_dir,
                                                          remote_path=self.embassy_apollo_rsync_path,
                                                          method="embassy")
        return (bash_command)


class Sample(Study):
    def __init__(self, ui, studies_dict, study_id, sample_id):
        super().__init__(ui, studies_dict, study_id)
        # Sample Ids
        self.sample_id = sample_id
        self.sample_dir = os.path.join(self.study_dir, self.sample_id)
        self.sample_dict = self.samples_dict[self.sample_id]
        self.tmp_dir = os.path.join(self.sample_dir, "tmp")
        self.short_sample_id = self.sample_id[0:6]

        # FASTQ Options
        self.fastqs_ftp_paths = self.sample_dict["fastq_ftp"].split(";")
        self.fastqs_fire_paths = [sra_ftp_to_fire(fastq) for fastq in self.fastqs_ftp_paths]
        self.fastqs_dict = {fastq: os.path.join(self.sample_dir, self.sample_id + ".r" + str(cnt) + ".fastq.gz") for
                            cnt, fastq in
                            enumerate(self.fastqs_fire_paths, 1)}
        self.fastqs = [os.path.join(self.sample_dir, self.sample_id + ".r" + str(cnt) + ".fastq.gz") for cnt, fastq in
                       enumerate(self.fastqs_fire_paths, 1)]
        self.trimmed_fastqs = [os.path.join(self.sample_dir, self.sample_id + ".r" + str(cnt) + ".trimmed.fastq.gz") for
                               cnt, fastq in enumerate(self.fastqs_fire_paths, 1)]
        self.fastp_html_report = os.path.join(self.sample_dir, self.sample_id + ".fastp_report.html")
        self.fastp_json_report = os.path.join(self.sample_dir, self.sample_id + ".fastp_report.json")
        self.no_of_fastqs = len(self.fastqs_ftp_paths)
        self.fastp_threads = fastp_threads
        self.fastp_memory = fastp_memory

        # STAR ALIGNMENT OPTIONS
        self.star_align_threads = star_align_threads
        self.star_align_memory = star_align_memory
        self.star_align_outFileNamePrefix = os.path.join(self.sample_dir, self.sample_id) + "."
        self.star_align_outFileNameSuffix = star_align_outFileNameSuffix
        self.star_align_limitBAMsortRAM = star_align_limitBAMsortRAM
        self.star_align_sjdbOverhang = star_align_sjdbOverhang
        self.star_align_extra_params = star_align_extra_params
        self.star_align_outbam = self.star_align_outFileNamePrefix + self.star_align_outFileNameSuffix
        self.star_align_outbai = self.star_align_outbam + ".bai"
        self.star_align_readFilesCommand = star_align_readFilesCommand
        self.star_align_limitBAMsortRAM = star_align_limitBAMsortRAM
        self.star_align_sjdbOverhang = star_align_sjdbOverhang

        # BAM PROCESSING OPTIONS
        self.bam_refname_sorted = os.path.join(self.sample_dir, self.sample_id + ".refname_sorted.bam")
        self.bam_capped = os.path.join(self.sample_dir, self.sample_id + ".capped.bam")
        self.bam_capped_sorted = os.path.join(self.sample_dir, self.sample_id + ".capped.sorted.bam")
        self.bigwig = os.path.join(self.sample_dir, self.sample_id + ".bigwig")
        self.bam2bigwig_binSize = bam2bigwig_binSize

        self.bam_for_apollo = self.bam_capped_sorted
        self.bigwig_for_apollo = self.bigwig

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

    def download_fastqs_command(self):
        bash_command = curl_fastqs(self.fastqs_dict)
        return (bash_command)

    def trim_fastqs_command(self, previous_files=None, delete_previous=False):
        if not previous_files: previous_files = self.fastqs
        if delete_previous and not isinstance(previous_files, list): previous_files = [previous_files]
        bash_command = trim_paired_fqs(fastp=fastp, input_fqs=self.fastqs,
                                       output_fqs=self.trimmed_fastqs,
                                       fastp_html_report=self.fastp_html_report,
                                       fastp_json_report=self.fastp_json_report,
                                       threads=self.fastp_threads) + \
                                       (delete_if_ok(" ".join(previous_files)) if delete_previous else "")
        return (bash_command)

    def star_alignment_command(self, previous_files=None, delete_previous=False):
        if not previous_files: previous_files = self.trimmed_fastqs
        if delete_previous and not isinstance(previous_files, list): previous_files = [previous_files]
        bash_command = star_alignment(star=star,
                                      reference=self.aligner_reference,
                                      threads=self.star_align_threads,
                                      fastqs=self.trimmed_fastqs,
                                      outFileNamePrefix=self.star_align_outFileNamePrefix,
                                      limitBAMsortRAM=self.star_align_limitBAMsortRAM,
                                      sjdbOverhang=self.star_align_sjdbOverhang,
                                      gtf=self.aligner_gtf,
                                      extra_params=self.star_align_extra_params) + \
                                      (delete_if_ok(" ".join(previous_files)) if delete_previous else "")
        return (bash_command)
    
    def index_aligned_bam_command(self, previous_files=None, delete_previous=False):
        if not previous_files: previous_files = self.trimmed_fastqs
        if delete_previous and not isinstance(previous_files, list): previous_files = [previous_files]
        bash_command = index_bam(inbam=self.star_align_outbam)
        return (bash_command)

    def sortrefname_bam_command(self, previous_files=None, delete_previous=False):
        if not previous_files: previous_files = self.star_align_outbam
        if delete_previous and not isinstance(previous_files, list): previous_files = [previous_files]
        bash_command = sortrefname_bam(self.star_align_outbam,
                                       self.bam_refname_sorted,
                                       self.tmp_dir,
                                       sortsamrefname=sortsamrefname) + \
                                       (delete_if_ok(" ".join(previous_files)) if delete_previous else "")
        return (bash_command)

    def cap_bam_command(self, previous_files=None, delete_previous=False):
        if not previous_files: previous_files = self.bam_refname_sorted
        if delete_previous and not isinstance(previous_files, list): previous_files = [previous_files]
        bash_command = cap_bam(self.bam_refname_sorted,
                               self.bam_capped,
                               cap_reads=cap_reads,
                               biostar154220=biostar154220) + \
                               (delete_if_ok(" ".join(previous_files)) if delete_previous else "")
        return (bash_command)
    
    def sort_capped_bam_command(self, previous_files=None, delete_previous=False):
        if not previous_files: previous_files = self.bam_capped
        if delete_previous and not isinstance(previous_files, list): previous_files = [previous_files]
        bash_command = positionsort_bam(inbam=self.bam_capped,
                                        outbam=self.bam_capped_sorted) + \
                                        (delete_if_ok(" ".join(previous_files)) if delete_previous else "")
        return (bash_command)

    def index_final_bam_command(self, previous_files=None, delete_previous=False):
        if not previous_files: previous_files = self.trimmed_fastqs
        if delete_previous and not isinstance(previous_files, list): previous_files = [previous_files]
        bash_command = index_bam(inbam=self.bam_for_apollo)
        return (bash_command)

    def bam2bigwig_command(self, previous_files=None, delete_previous=False):
        if not previous_files: previous_files = self.star_align_outbam
        if delete_previous and not isinstance(previous_files, list): previous_files = [previous_files]
        bash_command = bam2bigwig(self.star_align_outbam,
                                  self.bigwig,
                                  self.bam2bigwig_binSize,
                                  bamCoverage=bamCoverage) + \
                                  (delete_if_ok(" ".join(previous_files)) if delete_previous else "")
        return (bash_command)

    def move_final_bam_to_apollo_dir_command(self):
        to_move_file = self.bam_for_apollo
        bash_command = do_if_exist(to_move_file + ".bai", "\tmv " + \
                                   to_move_file + "* " + \
                                   self.to_apollo_dir)
        return (bash_command)

    def move_final_bigwig_to_apollo_dir_command(self):
        to_move_file = self.bigwig_for_apollo
        bash_command = do_if_exist(to_move_file, "\tmv " + \
                                   to_move_file + " " + \
                                   self.to_apollo_dir)
        return (bash_command)


#     def check_if_studies(self):
#         if not os.path.exists(self.se_dir):
#             exit_with_error(
#                 "No expression data found for " + self.species + " in " + self.se_dir + ".\n No more work to do. Exiting.")
#     def report_input(self):
#         dep_dirs_names = [item for item in dir(dependencies) if item.endswith('dir')]
#         print_info("Given Input (from user and dependencies.py):")
#         print_w_indent("Species: "+self.species)
#         for depdir in self.__dict__.keys():
#             pot_path=self.__dict__[depdir]
#             if depdir.endswith("dir"):
#                 print_w_indent(depdir + ": " + pot_path)
#     def studies(self):
#         self.check_if_studies()
#         study_ids = [x for x in next(os.walk(self.se_dir))[1]
#                       if os.path.isfile(self.se_dir + "/" + x + "/" + x + ".design.tsv")]
#         return study_ids
#     def samples(self):
#         return_samples = []
#         study_ids = self.studies()
#         for study_id in study_ids:
#             sst = Study(self.ui, study_id)
#             return_samples += list(set(sst.design_dict().values()))
#         return(return_samples)
#     def studies_to_samples(self):
#         stu2sam_dict = {}
#         study_ids = self.studies()
#         for study_id in study_ids:
#             sst = Study(self.ui, study_id)
#             stu2sam_dict[study_id] = list(set(sst.design_dict().values()))
#         return(stu2sam_dict)
#     def all_bams(self):
#         return(flatten([Study(self.ui, study_id).bams() for study_id in self.studies()]))
#     def all_final_merged_bams(self):
#         return([Study(self.ui, study_id).merged_bam_final for study_id in self.studies()])
#     def all_bigwigs(self):
#         return(flatten([Study(self.ui, study_id).bigwigs() for study_id in self.studies()]))
#     def move_all_finals_to_final_dir_command(self, to_wait_id=""):
#         bam_files_to_move = self.all_bams() + self.all_final_merged_bams()
#         bigwigs_to_move = self.all_bigwigs()
#         bash_commands = [do_if_exist(x, "\tmv " + \
#                          x + "* " + \
#                          self.toapollo_dir) for x in bam_files_to_move]
#         bash_commands += [do_if_exist(x, "\tmv " + \
#                          x + " " + \
#                          self.toapollo_dir) for x in bigwigs_to_move]
#         bash_command = "\n".join(bash_commands)
#         return(bash_command)
#     def move_all_finals_to_final_dir_command_submit(self, to_wait_id=""):
#         job_id = lsf_submit(self.move_all_finals_to_final_dir_command(),
#                             jobprefix=self.species + "_01_move_to_final_dir",
#                             to_wait_id=to_wait_id,
#                             cpu=1,
#                             mem="1gb",
#                             cwd=self.log_dir,
#                             queue="production")
#         return (job_id)
#     def rsync_to_external_location_command_submit(self, to_wait_id=""):
#         job_id = lsf_submit(self.rsync_to_external_location_command(),
#                             jobprefix=self.species + "_02_copy_to_embassy",
#                             to_wait_id=to_wait_id,
#                             cpu=1,
#                             mem="1gb",
#                             cwd=self.log_dir,
#                             queue="production")
#         return(job_id)
#
# class Study(Species):
#     def __init__(self, ui, study_id):
#         super().__init__(ui)
#         self.study_id = study_id
#         self.spath = os.path.join(self.se_dir, self.study_id)
#         self.design_file = os.path.join(self.spath, self.study_id + ".design.tsv")
#         self.study_bam_dir = os.path.join(self.bam_dir,self.study_id)
#         self.study_bigwig_dir = os.path.join(self.bigwig_dir,self.study_id)
#         self.ena_api_url = ena_api_url.format(self.study_id)
#         self.genome_fa_prepared = 0
#         self.ref_fasta = os.path.join(self.reference_dir, self.species+"."+self.study_id+".genomic.fa")
#         self.merged_dir = os.path.join(self.bam_dir,"merged")
#         self.merged_bam = os.path.join(self.merged_dir, self.study_id+".merged.bam")
#         self.merged_bam_noext = str(Path(self.merged_bam).with_suffix(''))
#         self.merged_bam_sortrefnamed = self.merged_bam_noext + ".sortrefnamed.bam"
#         self.merged_bam_capped = self.merged_bam_noext + ".TEMP.capped.bam"
#         self.merged_bam_final = self.merged_bam_noext + ".capped.bam"
#     def design_dict(self):
#         design_reader = csv.DictReader(open(self.design_file), delimiter="\t")
#         desdict = {self.study_id + '_' + re.sub(r'[^\w]', '', row['Condition']): row['Run'] for row in design_reader}
#         return desdict
#     def design_keys(self):
#         return list(self.design_dict().keys())
#     def samples(self):
#         return list(self.design_dict.values())
#     def create_study_dirs(self):
#         run_once = 0
#         for dir in self.__dict__.keys():
#             pot_path = self.__dict__[dir]
#             if dir.endswith("dir") and not os.path.exists(pot_path):
#                 if run_once == 0:
#                     print_info("Creating Required paths:")
#                 print_w_indent("Creating " + pot_path)
#                 os.makedirs(pot_path)
#                 run_once += 1
#     def ref_fasta_ftpfile(self):
#         for dek in self.design_keys():
#             ssa = Sample(self.ui, self.study_id, dek)
#             fasta = ssa.ref_fasta_ftpfile_sample()
#             break
#         return fasta
#     def is_fasta_gzipped(self):
#         return self.ref_fasta_ftpfile().endswith(".gz")
#     def ref_fasta_copied(self):
#         if self.is_fasta_gzipped():
#             return self.ref_fasta+".gz"
#         else:
#             return self.ref_fasta
#     def copy_ref_fasta(self):
#         if not os.path.exists(self.ref_fasta):
#             print_info("Copying reference fasta file "+self.ref_fasta_ftpfile()+" for study "+self.study_id+".")
#             sysrsync.run(
#                 source=self.ref_fasta_ftpfile(),
#                 destination=self.ref_fasta_copied(),
#                 options=['-avz', '--checksum'])
#             if self.is_fasta_gzipped():
#                 self.decompress_ref_fasta()
#     def decompress_ref_fasta(self):
#         decompress(self.ref_fasta_copied(), self.ref_fasta)
#         if os.path.isfile(self.ref_fasta):
#             os.remove(self.ref_fasta_copied())
#     def bigwigs(self):
#         bigwigs = []
#         for dek in self.design_keys():
#             ssa = Sample(self.ui, self.study_id, dek)
#             bigwig = ssa.bigwig
#             bigwigs.append(bigwig)
#         return bigwigs
#     def bams(self):
#         bams = []
#         for dek in self.design_keys():
#             ssa = Sample(self.ui, self.study_id, dek)
#             bam = ssa.bam_deduped
#             bams.append(bam)
#         return bams
#     def merge_and_cap_bams_command(self):
#         bash_command = merge_bams_command(self.bams(), self.merged_bam, 4) + \
#                        index_bam_command(self.merged_bam) + \
#                        sortsamrefname_command(self.merged_bam, self.merged_bam_sortrefnamed, sortsamrefname) + \
#                        delete_if_ok(self.merged_bam) + \
#                        capbam_command(self.merged_bam_sortrefnamed, self.merged_bam_capped, 50, biostar154220) + \
#                        delete_if_ok(self.merged_bam_sortrefnamed) + \
#                        positionsort_bam_command(self.merged_bam_capped, self.merged_bam_final) + \
#                        delete_if_ok(self.merged_bam_capped) + \
#                        index_bam_command(self.merged_bam_final)
#         return bash_command
#     def merge_and_cap_bams_command_submit(self, to_wait_id=""):
#         job_id = lsf_submit(self.merge_and_cap_bams_command(),
#                             jobprefix=self.study_id+"_01_merge_bams",
#                             to_wait_id=to_wait_id,
#                             cpu=4,
#                             mem="24gb",
#                             cwd=self.log_dir,
#                             queue="production")
#         return(job_id)
#
#
#
#
#
#
# class Sample(Study):
#     def __init__(self, ui, study_id, design_key):
#         super().__init__(ui, study_id)
#         self.sample_name = design_key
#         self.sample_id = self.design_dict()[self.sample_name]
#         self.shortid = self.sample_id[0:6]
#         self.sample_bam_dir = os.path.join(self.study_bam_dir,self.sample_id)
#         self.sample_bigwig_dir = os.path.join(self.study_bigwig_dir, self.sample_id)
#         self.bam = os.path.join(self.sample_bam_dir, self.sample_name+".bam")
#         self.bam_noext = str(Path(self.bam).with_suffix(''))
#         self.bam_nsorted = self.bam_noext + ".nsorted.bam"
#         self.bam_psorted = self.bam_noext + ".psorted.bam"
#         self.bam_fixmated = self.bam_noext + ".fixmated.bam"
#         self.bam_deduped = self.bam_noext + ".deduped.bam"
#         self.bigwig = os.path.join(self.sample_bigwig_dir, self.sample_name+".bigwig")
#         self.fasta = self.ref_fasta
#     def ftp(self):
#         ftp = os.path.join(rnaseq_ftp_dir, self.shortid, self.sample_id)
#         if not os.path.isdir(ftp):
#             alt_ftp_path = os.path.join(rnaseq_ftp_dir, self.shortid, "*", self.sample_id)
#             alt_ftp_dirs = glob.glob(alt_ftp_path)
#             if len(alt_ftp_dirs) > 1:
#                 exit_with_error("There are many FTP directories for "+alt_ftp_path+". Exiting.")
#             elif len(alt_ftp_dirs) == 0:
#                 exit_with_error("Cannot find an FTP directory in "+ftp+" or "+alt_ftp_path)
#             else:
#                 return(alt_ftp_dirs[0])
#         return(ftp)
#     def create_sample_dirs(self):
#         for dir_name in [self.sample_bam_dir, self.sample_bigwig_dir]:
#             if not os.path.exists(dir_name):
#                 print_info("Creating Required paths:")
#                 print_w_indent("Creating "+dir_name)
#                 os.makedirs(dir_name)
#     def bam_ftpfile(self):
#         bamfiles = glob.glob(os.path.join(self.ftp(),self.sample_id+".bam"))
#         cramfiles = glob.glob(os.path.join(self.ftp(),self.sample_id+".cram"))
#         if len(bamfiles)>0:
#             alfile = bamfiles[0]
#         elif len(bamfiles)<1 and len(cramfiles)>0:
#             alfile = cramfiles[0]
#         else:
#             exit_with_error("Could not find BAM/CRAM file in FTP"
#                             "directory: " + self.ftp() + " for " + self.sample_id)
#         return(alfile)
#     def bigwig_ftpfile(self):
#         bigwigfiles = glob.glob(os.path.join(self.ftp(),self.sample_id+".nospliced.bw"))
#         if len(bigwigfiles)<0:
#             exit_with_error("Could not find a bigwig file in FTP"
#                             "directory: " + self.ftp() + " for " + self.sample_id)
#         return(bigwigfiles[0])
#     def bam_copied(self):
#         if self.bam_ftpfile().endswith(".cram"):
#             return str(Path(self.bam).with_suffix('.cram'))
#         else:
#             return self.bam
#     def is_cram(self):
#         return self.bam_ftpfile().endswith(".cram")
#     def bioproject(self):
#         pyal = pysam.AlignmentFile(self.bam_ftpfile(), "r")
#         fasta = list(set([x['UR'] for x in pyal.header.to_dict()['SQ']]))[0]
#         fasta_basename = os.path.basename(fasta)
#         bioproject = fasta_basename.split(".")[1].strip()
#         return(bioproject)
#     def ref_fasta_ftpfile_sample(self):
#         fasta = glob.glob(os.path.join(ps_ftp_dir,self.species,self.bioproject().upper(),"*.genomic.fa.*"))[0]
#         check_if_file_exists(fasta)
#         #print(os.path.join(ps_ftp_dir,self.species,self.bioproject().upper(),"*.genomic.fa.*"))
#         return(fasta)
#     def copy_bam_from_ftp_command(self):
#         bash_command=rsync_command(self.bam_ftpfile(), self.bam_copied())
#         return(bash_command)
#     def copy_bigwig_from_ftp_command(self):
#         bash_command = rsync_command(self.bigwig_ftpfile(), self.bigwig)
#         return(bash_command)
#     def copy_command(self):
#         bash_command = bash_check_if_success_or_loop(self.copy_bam_from_ftp_command(),times=3) + \
#                        bash_check_if_success_or_loop(self.copy_bigwig_from_ftp_command(),times=3)
#         return(bash_command)
#     def bam_processing_command(self):
#         bash_command = "" + \
#                        (cram2bam_command(self.bam_copied(), self.bam, self.fasta) +
#                             delete_if_ok(self.bam_copied()) if self.is_cram() else "") + \
#                        namesort_bam_command(self.bam, self.bam_nsorted) + \
#                        delete_if_ok(self.bam) + \
#                        fixmate_bam_command(self.bam_nsorted, self.bam_fixmated) + \
#                        delete_if_ok(self.bam_nsorted) + \
#                        positionsort_bam_command(self.bam_fixmated, self.bam_psorted) + \
#                        delete_if_ok(self.bam_fixmated) + \
#                        markdup_bam_command(self.bam_psorted, self.bam_deduped) + \
#                        delete_if_ok(self.bam_psorted) + \
#                        index_bam_command(self.bam_deduped)
#         return bash_command
#     def copy_command_submit(self, to_wait_id=""):
#         job_id = lsf_submit(self.copy_command(),
#                             jobprefix=self.sample_name+"_01_bam_copy",
#                             to_wait_id=to_wait_id,
#                             cpu=1,
#                             mem="1gb",
#                             cwd=self.log_dir,
#                             queue="datamover")
#         return(job_id)
#     def bam_processing_command_submit(self, to_wait_id=""):
#         job_id = lsf_submit(self.bam_processing_command(),
#                             jobprefix=self.sample_name+"_02_bam_processing",
#                             to_wait_id=to_wait_id,
#                             cpu=1,
#                             mem="12gb",
#                             cwd=self.log_dir,
#                             queue="production")
#         return(job_id)


def explore_print(ui):
    print("Species" + "\t" + "Study_id" + "\t" + "Sample_id" + "\t" + "sample_id_condition" + "\t" + "ftp_location")
    for spec_study in species.studies():
        sst = Study(ui, spec_study)
        for dek in sst.design_keys():
            ssa = Sample(ui, spec_study, dek)
            print(ssa.species + "\t" + ssa.study_id + "\t" + ssa.sample_id + "\t" + ssa.sample_name + "\t" + ssa.ftp())

def all_studies_dict_for_species(get_rna_seq_studies_by_taxon_api_url):
    url = get_rna_seq_studies_by_taxon_api_url
    return (ena_rnaseq_by_taxon_id_to_dict(url))

def selected_studies_dict_for_species(ui, all_studies_dict):
    if not ui.selects:
        print_warning("No selected studies - All studies will be processed.")
        time.sleep(30)
        return (self.all_studies_dict)
    selected_studies_dict = {}
    for select in ui.selects:
        study_id = select.split(":")[0]
        if study_id in all_studies_dict.keys():
            selected_studies_dict[study_id] = all_studies_dict[study_id]
            all_sample_ids = [x["run_accession"] for x in selected_studies_dict[study_id]]
            try:
                sample_ids = select.split(":")[1].split(";")
            except IndexError:
                print_info("No samples were identified for study: " + study_id
                           + ". Processing with all samples.")
                sample_ids = all_sample_ids
            for sample_id in sample_ids:
                if sample_id not in all_sample_ids:
                    exit_with_error("The given sample ID: " + sample_id + " is not a valid "
                                                                          "sample id from the study " + study_id + ". Please double-check. Exiting.")
            study_dict_values = [x for x in selected_studies_dict[study_id] if x["run_accession"] in sample_ids]
            selected_studies_dict[study_id] = study_dict_values
        else:
            exit_with_error("The given study ID: " + study_id + " is not a valid "
                                                                "study id for that species")
    return (selected_studies_dict)