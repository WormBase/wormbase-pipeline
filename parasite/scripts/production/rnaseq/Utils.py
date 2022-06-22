import requests
from ProductionUtils import *
from dependencies import *

def ena_rnaseq_by_taxon_id_to_dict(url):
    r = requests.get(url)
    r_json = r.json()
    study_dict = {}
    for sample_accession in r_json:
        if sample_accession["secondary_study_accession"] not in study_dict.keys():
            study_dict[sample_accession["secondary_study_accession"]] = []
        study_dict[sample_accession["secondary_study_accession"]].append(sample_accession)
    return(study_dict)

def sra_ftp_to_fire(fastq, sra_ftp_path_prefix=sra_ftp_path_prefix, sra_fire_path_prefix=sra_fire_path_prefix):
    fastq = fastq.strip()
    fire_fastq = fastq.replace(sra_ftp_path_prefix, sra_fire_path_prefix)
    return(fire_fastq)
    # if url_file_exists("http://" + fire_fastq):
    #     return (fire_fastq)
    # elif url_file_exists("http://" + fastq):
    #     return (fastq)
    # else:
    #     exit_with_error("The remote fastq file does not exist in any of its forms:\n" +
    #                     "\n".join([fire_fastq, fastq]))

def validate_selects_format(user_select):
    validation_regex = "^\w+$|^\w+:\w+$|^\w+:\w+(;\w+){1,}$"
    return (bool(re.search(validation_regex, user_select)))

def species_input(ui_species):
    if bool(re.match("^[a-z0-9]+_[a-z0-9]+_[a-z0-9]+$", ui_species)):
        return ("spe_cies_bp")
    elif bool(re.match("^[a-z0-9]+_[a-z0-9]+$", ui_species)):
        return ("spe_cies")
    else:
        exit_with_error("Wrong species input. Please try again.")