import os


expression_dir = "/hps/software/users/wormbase/parasite/repositories/wbps-expression/studies"
rnaseq_ftp_dir = "/nfs/ftp/public/databases/arrayexpress/data/atlas/rnaseq"
ps_ftp_dir = "/nfs/ftp/public/databases/wormbase/parasite/releases/current/species/"
apollo_dir = os.getenv("PARASITE_SCRATCH")+"/apollo"
ena_api_url = "https://www.ebi.ac.uk/ena/portal/api/search?result=study&query=secondary_study_accession={0}&fields=study_accession&format=json"
permanent_apollo_dir="/nfs/production/flicek/wormbase/parasite/apollo"


#CONNECTIONS
external_ssh = "sangerngs"
external_ssh_path = "/data/production/parasites/apollo"
ssh_host = external_ssh + ":" + external_ssh_path

#SOFTWARE
wgtools_image = "/hps/software/users/wormbase/parasite/images/wiggletools.sif"
wgtools = "singularity run "+wgtools_image
biostar154220_path = "/hps/software/users/wormbase/parasite/software/jvarkit/dist/biostar154220.jar"
sortsamrefname_path = "/hps/software/users/wormbase/parasite/software/jvarkit/dist/sortsamrefname.jar"
biostar154220 = "java -jar " + biostar154220_path
sortsamrefname = "java -jar " + sortsamrefname_path