import os


expression_dir = "/hps/software/users/wormbase/parasite/repositories/wbps-expression/studies"
#rnaseq_ftp_dir = "/nfs/ftp/public/databases/arrayexpress/data/atlas/rnaseq"
#ps_ftp_dir = "/nfs/ftp/public/databases/wormbase/parasite/releases/current/species/"
rnaseq_dir = os.getenv("PARASITE_SCRATCH")+"/rnaseq"
ena_api_url = "https://www.ebi.ac.uk/ena/portal/api/search?result=study&query=secondary_study_accession={0}&fields=study_accession&format=json"
ena_rnaseq_by_taxon_url = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=tax_eq({0})%20AND%20(library_strategy%3D%22RNA-Seq%22%20OR%20library_strategy%3D%22OTHER%22)&fields=study_accession%2Caccession%2Cfastq_ftp%2Csecondary_study_accession%2Csecondary_sample_accession&format=json"
ena_rnaseq_by_taxon_url_onlyrnaseq = "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=tax_eq({0})%20AND%20(library_strategy%3D%22RNA-Seq%22)&fields=study_accession%2Caccession%2Cfastq_ftp%2Csecondary_study_accession%2Csecondary_sample_accession&format=json"
ena_secondary_study_id_count_url = "https://www.ebi.ac.uk/ena/portal/api/count?dataPortal=ena&query=secondary_study_accession%3D{0}&result=study"
ena_run_accession_count_url = "https://www.ebi.ac.uk/ena/portal/api/count?dataPortal=ena&query=run_accession%3D{0}&result=read_run"
permanent_apollo_dir="/nfs/production/flicek/wormbase/parasite/apollo"
reference_genomes_dir="/nfs/production/flicek/wormbase/parasite/data/reference_genomes/star/"
reference_genomes_ftp="https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/species/"
sra_ftp_path_prefix="ftp.sra.ebi.ac.uk/vol1/fastq/"
sra_fire_path_prefix="hh.fire.sdo.ebi.ac.uk/fire/public/era/fastq/"

#RUN PARAMETERS
fastp_memory = "10000"
fastp_threads = "6"

star_align_memory = "12000"
star_align_threads = "6"
star_align_readFilesCommand = "zcat"

star_align_outSAMtype = "BAM SortedByCoordinate"
star_align_outFileNameSuffix = "Aligned.sortedByCoord.out.bam"

star_align_limitBAMsortRAM = "18486355837"
star_align_sjdbOverhang = "99"
star_align_sjdbGTFtagExonParentGene = "gene_id"
star_align_quantMode = "GeneCounts"
star_align_outFilterMultimapNmax = "1"
star_align_extra_params = "--outSAMtype " + star_align_outSAMtype + " " + \
                          "--readFilesCommand " + star_align_readFilesCommand + " " + \
                          "--sjdbGTFtagExonParentGene " + star_align_sjdbGTFtagExonParentGene + " " + \
                          "--quantMode " + star_align_quantMode + " " #+ \
                          #"--outFilterMultimapNmax " + star_align_outFilterMultimapNmax + " "
                    
cap_reads = "30"

bam2bigwig_binSize = "10"

#CONNECTIONS
external_ssh = "sangerngs"
external_ssh_path = "/data/production/parasites/apollo"
ssh_host = external_ssh + ":" + external_ssh_path
try:
    embassy_bucket = os.environ["EMBASSY_BUCKET"]
    embassy_apollo_path = os.environ["EMBASSY_APOLLO_PATH"]
    embassy_access_url_apollo = os.environ["EMBASSY_ACCESS_URL_APOLLO"]
except KeyError:
    print("The embassy module has not been loaded. Maybe you want to load it and re-run?")


#SOFTWARE
sw_wgtools_image_path = "/hps/software/users/wormbase/parasite/images/wiggletools.sif"
wgtools = "singularity run "+sw_wgtools_image_path
sw_biostar154220_path = "/hps/software/users/wormbase/parasite/software/jvarkit/dist/biostar154220.jar"
sw_sortsamrefname_path = "/hps/software/users/wormbase/parasite/software/jvarkit/dist/sortsamrefname.jar"
biostar154220 = "java -jar " + sw_biostar154220_path
sortsamrefname = "java -jar " + sw_sortsamrefname_path
sw_fastp_path = os.getenv("PARASITE_SOFTWARE") + "/fastp"
fastp = sw_fastp_path
star = "STAR"
sw_path_path = "/hps/software/users/wormbase/parasite/shared/.pyenv/versions/p395/bin"
bamCoverage = 'export PATH='+sw_path_path+':$PATH\nexport PYTHONPATH=""\n\nbamCoverage'