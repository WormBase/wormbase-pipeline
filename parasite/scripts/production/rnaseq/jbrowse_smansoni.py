import os
import subprocess
import json
import requests
from ProductionUtils import *

EMBASSY_RNASEQER_PATH = os.environ["EMBASSY_APOLLO_PATH"]
eaws = os.environ["EMBASSY_COMMAND"]
APOLLO_INSTANCE = "schistosoma_mansoni_v9_ENA"
APOLLO_PATH = EMBASSY_RNASEQER_PATH + "/" + APOLLO_INSTANCE
EMBASSY_URL = os.environ["EMBASSY_URL"]
EMBASSY_BUCKET = os.environ["EMBASSY_BUCKET"]
ENA_RUN_API = 'https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=run_accession%3D%22{0}%22&fields=secondary_study_accession%2Cstudy_accession%2Crun_accession%2Csample_accession%2Csample_title%2Csample_alias%2Ccenter_name%2Cstudy_alias%2Cstudy_title&format=json'
JBROWSE_INSTALL_DIR="/hps/software/users/wormbase/parasite/software/jbrowse"
JBROWSE_OUT_DIR="$PARASITE_SCRATCH/jbrowse/WBPS${PARASITE_VERSION}/out"
TEMP_DIR="$PARASITE_SCRATCH/temp"
studies_metadata_tsv = "/nfs/production/flicek/wormbase/parasite/data/releases/release17/schistosoma_mansoni_prjea36577/schistosoma_mansoni_jbrowse_tracks.txt"
out_studies_json = os.environ["PARASITE_SCRATCH"]+"/jbrowse/WBPS"+os.environ["PARASITE_VERSION"]+"/WbpsExpression/schistosoma_mansoni_prjea36577/schistosoma_mansoni.studies.json"

bigwigs_eaws = subprocess.run([eaws+" s3 ls "+APOLLO_PATH+"/ --recursive | awk '{print $4}' | grep '.bigwig'"], shell=True, stdout=subprocess.PIPE)
bigwigs = list(set([x for x in bigwigs_eaws.stdout.decode('utf-8').split('\n') if x != '']))

studies_metadata = csvlines2list(studies_metadata_tsv, delimiter="\t", skiplines=1)
studies_metadata_dict = {}
studies_metadata_dict = {x[0]:{} for x in studies_metadata}
for row in studies_metadata:
    studies_metadata_dict[row[0]].update({row[1]:{"curated_label":row[4], "developmental_stage":row[5], "pubmed":row[6], "sex":row[7], "strain":row[8], "center":row[9]}})

#print(studies_metadata_dict)
#
studies_dict = {}
for bw in bigwigs:
    bw_url = EMBASSY_URL + "/" + EMBASSY_BUCKET + "/" + bw
    study = bw.split("/")[-2]
    run = bw.split("/")[-1].split(".")[0]
    ena_run_api = ENA_RUN_API.format(run)
    r = requests.get(ena_run_api)
    r_json = r.json()
    sample_accession = r_json[0]["sample_accession"]
    sample_alias = r_json[0]["sample_alias"]
    sample_title = r_json[0]["sample_title"]
    study_title = r_json[0]["study_title"]
    submitting_center = r_json[0]["center_name"]
    bioproject = r_json[0]["study_accession"]
    #print("\t".join([study, run, sample_title, sample_alias, study_title, submitting_center]))
    curated_label = studies_metadata_dict[study][run]["curated_label"]
    developmental_stage = studies_metadata_dict[study][run]["developmental_stage"]
    pubmed_id = studies_metadata_dict[study][run]["pubmed"]
    pubmed = '<a href=\"https://www.ncbi.nlm.nih.gov/pubmed/'+pubmed_id+'\">'+pubmed_id+'</a>'
    sex = studies_metadata_dict[study][run]["sex"]
    strain = studies_metadata_dict[study][run]["strain"]
    center = studies_metadata_dict[study][run]["center"]
    if study not in studies_dict.keys():
        studies_dict[study] = {"study_id":study,
                                 "runs":[],
                                 "attributes":{
                                     "ENA BioProject":"<a href=\"https://www.ebi.ac.uk/ena/data/view/"+bioproject+"\">Study page: "+bioproject+"</a>",
                                     "ENA Study": "<a href=\"https://www.ebi.ac.uk/ena/data/view/"+study+"\">Study page: "+study+"</a>",
                                     "Study description": study_title,

                                 }}
        if pubmed!="NA": studies_dict[study]["attributes"].update({"pubmed":pubmed})
        if center!="NA": studies_dict[study]["attributes"].update({"submitting_centre":center})
        if strain!="NA": studies_dict[study]["attributes"].update({"strain":strain})
    run_dict = {"replicate":sample_accession, "bigwig":bw_url, "attributes":{}, "run_id":run, "condition":curated_label}
    if strain != "NA": run_dict["attributes"].update({"strain": strain})
    if sex != "NA": run_dict["attributes"].update({"sex": sex})
    if developmental_stage != "NA": run_dict["attributes"].update({"developmental_stage": developmental_stage})
    studies_dict[study]["runs"].append(run_dict)

studies_dict_list = [studies_dict[x] for x in studies_dict]

with open(out_studies_json, "w") as final:
    json.dump(studies_dict_list, final)

    # bw_cmd = "perl "+JBROWSE_INSTALL_DIR+"/bin/add-bw-track.pl " \
    #           "--in "+"/homes/digri/trackList.json " \
    #           "--out "+"/homes/digri/trackList.json " \
    #           "--bw_url " + bw_url + " " + \
    #           "--label " + run + ":" + curated_label + " " + \
    #           "--config '" + '{"metadata": { ' + \
    #                             '"Track": "'+run+':'+curated_label+'",' + \
    #                             '"Category": "RNASeq",' + \
    #                             '"Developmental stage": "'+developmental_stage+'",' \
    #                             '"ENA Study": "<a href=https://www.ebi.ac.uk/ena/browser/view/'+study+'>Study page: '+study+'</a>"' + \
    #                             (','+'"sex": "'+sex+'"' if sex !="NA" else "") + \
    #                             (',' + '"strain": "' + strain + '"' if strain != "NA" else "") + \
    #                             (','+'"pubmed": "'+pubmed+'"' if pubmed !="NA" else "") + \
    #                             (',' + '"submitting_centre": "' + center + '"' if center != "NA" else "") + \
    #                             ' }}' + "'"
    # print(bw_command)
    # bw_command = subprocess.run([bw_cmd], shell=True, stdout=subprocess.PIPE)

