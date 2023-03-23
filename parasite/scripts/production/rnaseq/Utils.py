import requests
from ProductionUtils import *
from dependencies import *
import time

def ena_rnaseq_by_taxon_id_to_dict(url):
    r = requests.get(url)
    r_json = r.json()
    study_dict = {}
    for sample_accession in r_json:
        if sample_accession["secondary_study_accession"] not in study_dict.keys():
            study_dict[sample_accession["secondary_study_accession"]] = []
        study_dict[sample_accession["secondary_study_accession"]].append(sample_accession)
    return(study_dict)

def is_ena_secondary_study(study_id, ena_secondary_study_id_count_url=ena_secondary_study_id_count_url):
    url = ena_secondary_study_id_count_url.format(study_id)
    r = requests.get(url)
    counter = 1
    while counter < 5:
        try:
            r_json = r.json()
            break
        except ValueError:
            counter += 1
            time.sleep(10)
            r_json=0
    if r_json >= 1:
        return True
    else:
        print_warning(study_id + " is not a valid ENA secondary_study_id.")
        return False

def is_ena_run_accession(run_accession, ena_run_accession_count_url=ena_run_accession_count_url):
    url = ena_run_accession_count_url.format(run_accession)
    r = requests.get(url)
    counter = 1
    while counter < 5:
        try:
            r_json = r.json()
            break
        except ValueError:
            counter += 1
            time.sleep(10)
            r_json=0
    if r_json >= 1:
        return True
    else:
        print_warning(run_accession + " is not a valid ENA run_accession.")
        return False

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
    validation_regex = "^\w+$|^\w+:\w+$|^\w+:\[?\w+(]?;\[?\w+){1,}]?|^\w+$"
    if "[" in user_select:
        bracket_error=False
        starting_positions = findOccurrences(user_select, "[")
        ending_positions = findOccurrences(user_select, "]")
        for count, value in enumerate(starting_positions):
            if starting_positions[count] >= ending_positions[count]:
                bracket_error=True
                break
    if bool(re.search(validation_regex, user_select))==False or bracket_error==True:
        return(False)
    else:
        return(True)

def species_input(ui_species):
    if bool(re.match("^[a-z0-9]+_[a-z0-9]+_[a-z0-9]+$", ui_species)):
        return ("spe_cies_bp")
    elif bool(re.match("^[a-z0-9]+_[a-z0-9]+$", ui_species)):
        return ("spe_cies")
    else:
        exit_with_error("Wrong species input. Please try again.")

def outlier_is_stranded_runs(production_host, db):

    production = Production(PRODUCTION_HOST=production_host, pattern=db)
    get_msg_res = [x.msg for x in production.connect().execute("SELECT DISTINCT(msg) FROM msg WHERE msg LIKE \"%Could not create a consensus:%\";")]
    msg_list = [{"metrics": list(y["groups"].groups()),
                 "runsdict": y["runsdict"]}
                for y in
                    [{"groups":re.search("metadata '(is_stranded)' \((\w+)=(\w+), (\w+)=(\w+)\)",x),
                                         "runsdict": json.loads(x.split("All metadata: ")[1])} for x in get_msg_res]
                if y["groups"] is not None]

    msg_list = [ {"metrics": x["metrics"],
                  "minority_is_str": x["metrics"][3] if x["metrics"][4] < x["metrics"][2] else x["metrics"][1],
                  "runsdict": x["runsdict"]} for x in msg_list]

    sep_runs = [[x for x in y["runsdict"] if str(y["runsdict"][x]["is_stranded"])==str(y["minority_is_str"])] for y in msg_list]

    return(sep_runs)

def all_runs_with_amgiguous_strand_inference(production_host, db):

    production = Production(PRODUCTION_HOST=production_host, pattern=db)
    get_msg_res = [x.msg for x in production.connect().execute("SELECT DISTINCT(msg) FROM msg WHERE msg LIKE \"%Could not create a consensus:%\";")]
    msg_list = [{"metrics": list(y["groups"].groups()),
                 "runsdict": y["runsdict"]}
                for y in
                    [{"groups":re.search("metadata '(is_stranded)' \((\w+)=(\w+), (\w+)=(\w+)\)",x),
                                         "runsdict": json.loads(x.split("All metadata: ")[1])} for x in get_msg_res]
                if y["groups"] is not None]

    msg_list = [ {"metrics": x["metrics"],
                  "minority_is_str": x["metrics"][3] if x["metrics"][4] < x["metrics"][2] else x["metrics"][1],
                  "runsdict": x["runsdict"]} for x in msg_list]

    sep_runs = [[x for x in y["runsdict"]] for y in msg_list]

    return(sep_runs)

def all_runs_with_failed_strand_inference(production_host, db):

    production = Production(PRODUCTION_HOST=production_host, pattern=db)
    get_msg_res = [x.msg for x in production.connect().execute("SELECT DISTINCT(msg) FROM msg WHERE msg LIKE \"%Ambiguous strand inference%\";")]

    runs = list(set([list(y["groups"].groups())[1]
                for y in
                    [{"groups":re.search("(\w+)_(\w+)_\w+_all", x),
                      "text":x} for x in get_msg_res]
                if y["groups"] is not None]))

    return(runs)