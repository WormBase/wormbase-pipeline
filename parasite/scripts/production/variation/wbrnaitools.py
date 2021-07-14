from __future__ import print_function
from intermine.webservice import Servicedd
import sys
import os
import pandas as pd
import re
import requests
import traceback
import datetime



def dtnow():
    ppdtnow = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return(str(ppdtnow))

def wpaXref_parser(file_url):
    """This function reads the WormBase WpaXref Paper ID mapping file as a pandas dataframe
    process the data so the dataframe will have one row per wormbase paper_id and one column
    with a list with all the identifiers linked to each paper. It then returns the pandas d
    dataframe object.
    :param file_url: String. Url of the WpaXref file. Defaults to WPAXREF_URL varialbe set in
    the dependencies file.
    :return: Pandas dataframe. Processed WormBase's Paper ID mapping file table.
    """
    #Read url file
    try:
        wpaxref_df = pd.read_csv(file_url, header=None, names=['wbp', 'map_id'], sep='\t')
    except urllib.error.HTTPError:
        print(dtnow() + ": ERROR - HTTPError Returned.\nReason: Cannot parse WormBase WpaXref file from"
              "here {0}. Please make sure the url in the dependencies.py file is "
              "correct and try again.".format(file_url))
        raise ValueError
    #Check the format of the df
    if list(set(wpaxref_df['wbp'].str.contains('^WBPaper.+').to_list())) != [True]:
        print(dtnow() + ": ERROR - File format Error.\nThe file in this address {0} has not been loaded"
              "correctly: Not all first column values follow the WBPaper ID regex.\n"
              "Please check the URL manually and try again.".format(file_url))
        raise ValueError
    #Process and return the WpaXref df
    wpaxref_df['map_id'] = wpaxref_df['map_id'].str.replace("<BR>", "")
    wpaxref_df = wpaxref_df.groupby('wbp')['map_id'].apply(list).reset_index(name='ids')
    return(wpaxref_df)

def wbphenoftp_parser(file, taxon_id, VAR_URL, RNAI_URL, names=None):
    """
    This function reads the WormBase phenotypes_association file
    from the given WormBase ftp path using the provided column
    names, it then processes the file to make it compatible with the
    ImportWORMBASE.pm of the phenotype import pipeline of
    ensembl-variation team.
    :param file: String. File path to the WB FTP phenotypes dump file.
    :param names: List. Column names to be used to name the output pandas df.
                  Defaults to ['gene_id','gene_symbol',
                             'pheno_type','pheno_id','source_id',
                             'var_id','object_type',
                             'taxon'].
    :return: A Pandas df object.
    """
    #Default behaviour
    if names==None: names = ['gene_id','gene_symbol',
                             'pheno_type','pheno_id','source_id',
                             'var_id','object_type',
                             'taxon']
    wbpheno_df = pd.read_csv(file, skiprows=3, header=None, sep='\t')[[1,2,3,4,5,7,11,12]]
    wbpheno_df.columns = names
    wbpheno_df = wbpheno_df[wbpheno_df['taxon']=='taxon:'+taxon_id]
    wbpheno_df['source_id'] = wbpheno_df['source_id'].apply(lambda x: x.split(':')[1])
    wbpheno_df['source_id'] = wbpheno_df['source_id'].apply(lambda x: check_source_id_format(x))
    wbpheno_df['new_var_id'] = wbpheno_df['var_id'].apply(lambda z: [x.replace('WB:', '') for x in str(z).split('|')])
    wbpheno_df_extended = wbpheno_df.explode('new_var_id').drop_duplicates()
    wbpheno_df_extended['source_type'] = wbpheno_df_extended.apply(lambda x: decide_source(x), axis=1)
    wbpheno_df_extended['final_source_id'] = wbpheno_df_extended.apply(lambda x: decide_source_id(x), axis=1)
    return(wbpheno_df_extended)

def wbrnaidump_parser(file,
                      names=None):
    """
    This function reads the WormBase dump RNAi phenotypes dump
    file from the given path using the provided column
    names, it then processes the file to make it compatible with the phenotypes table of the
    ensembl-variation database and returns it as a pandas df.
    :param file: String. File path to the WB RNAi phenotypes dump file. Defaults to RNAI_FILE.
    :param names: List. Column names to be used to name the output pandas df.
                  Defaules to ['gene','sequence','phenotype_name','phenotype_id','rnai:study'].
    :return: A Pandas df object.
    """
    #Default behaviour
    if names==None: names = ['gene','sequence','phenotype_name','phenotype_id','rnai:study']
    wbrnai_df = pd.read_csv(file, names=names, header=None, sep='\t')
    wbrnai_df['papers'] = wbrnai_df['rnai:study'].apply(lambda z: list(set([x.split('|')[1] for x in z.split(' ')])))
    wbrnai_df_extended = wbrnai_df.explode('papers').drop_duplicates()
    wbrnai_df_extended = wbrnai_df_extended[wbrnai_df_extended['papers'] != '']
    return(wbrnai_df_extended)

def wbpaperurl2json(wbpid, wbpaperapiurl, wbpmethods=None):
    """
    This function returns the results of specific methods of the WormBase Paper API results
    (https://wormbase.org/about/userguide/for_developers/API-REST/Paper#0--10)
     for a given WormBase Paper ID
    Args:
        :param wbpid (required): String. WormBase Paper ID.
        :param wbmethods (optional): String or List. Method(s) of the WormBase Paper API which will be
                          used in the search. Default: ['title','authors','year','pmid']
    :return: Dictionary with the requested data as keys and the corresponding API returned in json format
             as values.
    """
    wbp_json_dict = {}
    if wbpmethods is None:
        wbpmethods = ['title','authors','year','pmid','doi']
    if not isinstance(wbpmethods, list):
        wbpmethods = [wbpmethods]
    for method in wbpmethods:
        wbpurl = wbpaperapiurl.format(wbpid, method)
        response = requests.get(wbpurl)
        try:
            wbp_json = response.json()
        except json.decoder.JSONDecodeError:
            print(dtnow() + ": ERROR - ERROR: json.decoder.JSONDecodeError:")
            print("Error parsing {0} for {1}.".format(wbpurl, wbpid))
            print("Please Try again.\n")
            raise ValueError
        wbp_json_dict[method] = wbp_json[method]
    return(wbp_json_dict)

def wormmine_pheno_parser(url):
    """
    This function queries WormBase's WormMine for phenoytpes and returns a pandas dataframe
    with phenotype ids, names and descriptions
    :param url (string): WormMine query URL
    :return: pandas dataframe with three columns: WormBase Phenotype ID, Phenotype Name, Phenotype ID
    """
    service = Service(url)
    query = service.new_query("Phenotype")
    query.add_view("identifier", "name", "description")
    wm_pheno_dict = {}
    for row in query.rows():
        wm_pheno_dict[row["identifier"]] = [row["identifier"], row["name"], row["description"]]
    wm_pheno_df = pd.DataFrame.from_dict(wm_pheno_dict, orient='index')
    wm_pheno_df.columns = ['id','name','description']
    wm_pheno_df = wm_pheno_df.reset_index(drop=True)
    return(wm_pheno_df)


def paperid2varstudydict(wbpid, wbpaperapiurl, ncbi_url, wbpmethods=None):
    """
    This function returns the results of methods title,authors,years and pmid from the WormBase Paper API
    (https://wormbase.org/about/userguide/for_developers/API-REST/Paper#0--10) for a given WormBase Paper ID
    in a format suitable for each row of the study-file which is used as an input for the  import_wormbase class
    of the PhenotypeAnnotation pipeline in the ensembl-variation repository
    (https://github.com/Ensembl/ensembl-variation/tree/release/104/modules/Bio/EnsEMBL/Variation/Pipeline/PhenotypeAnnotation).

    Args:
        :param wbpid (required): String. WormBase Paper ID.
        :param wbmethods (optional): String or List. Method(s) of the WormBase Paper API which will be
                          used in the search. Default: ['title','authors','year','pmid']

    :return: Dictionary with the requested data as keys and the corresponding API returned in json format
             as values.
    """
    wb_study_dict = {}
    wb_json_dict = wbpaperurl2json(wbpid, wbpaperapiurl)
    try:
        wb_study_dict['title'] = wb_json_dict['title']['data']
        wb_study_dict['name'] = "{0} et al., {1}".format(wb_json_dict['authors']['data'][0]['label'], wb_json_dict['year']['data'])
        wb_study_dict['pmid'] = wb_json_dict['pmid']['data']
        if wb_study_dict['pmid'] == None:
            wb_study_dict['pmid'] = papertitle2pmidncbi(wb_study_dict['title'], ncbi_url)
        print(dtnow() + ": PROGRESS - Got PMID:{0} for Paper_ID:{1}".format(wb_study_dict['pmid'], wbpid))
    except KeyError:
        print(dtnow() + ": ERROR - Some key is not in the API object for {0}\n".format(wb_json_dict))
        raise KeyError
    return(wb_study_dict)

def varid2varstudydict(wbvid, wbvmethods=None):
    """This Function takes a WormBase Variant ID as its input and returns a dictionary
       similar to the paperid2varstudydict output."""
    wb_study_dict = {}
    wb_study_dict['title'] = wbvid
    wb_study_dict['name'] = wbvid
    wb_study_dict['pmid'] = wbvid
    print(dtnow() + ": PROGRESS - Got PMID:{0} for VarID:{0}".format(wbvid))
    return(wb_study_dict)

def sourceid2varstudydict(source_id, wbpaperapiurl, ncbi_url):
    if source_id.startswith('WBPaper'):
        return(paperid2varstudydict(source_id, wbpaperapiurl, ncbi_url))
    elif source_id.startswith('WBVar'):
        return(varid2varstudydict(source_id))
    else:
        print(dtnow() + ": ERROR:")
        print("WB ID {0} doesn't have a proper WB ID format".format(source_id))
        print("Please Try again.\n")
        raise ValueError

def sourceids2studyfiledf(unique_paper_ids, wbpaperapiurl, ncbi_url):
    """
    Function that takes a list of source IDs in the form WBPaperXXXXXXXX/WBVarXXXXXXXX as its input
    and returns a pandas df with the necessary
    columns required to populate a variation db.
    :param unique_paper_ids: List. List of unique WB Paper IDs in the form of WBPaperXXXXXXXX.
    :return: pandas df to create the phenotype_studies.txt file.
    """
    study_file_dict = {x: sourceid2varstudydict(x, wbpaperapiurl, ncbi_url) for x in unique_paper_ids}
    study_file_df = pd.DataFrame.from_dict(study_file_dict, orient='index')
    study_file_df['wbp'] = study_file_df.index
    study_file_df = study_file_df[['wbp','pmid','name','title']]
    return(study_file_df)

def wbphenourl2json(wbpid, wbphenoapiurl, wbpmethods=None):
    """
    This function returns the results of specific methods of the WormBase Paper API results
    (https://wormbase.org/about/userguide/for_developers/API-REST/Paper#0--10)
     for a given WormBase Paper ID
    Args:
        :param wbpid (required): String. WormBase Paper ID.
        :param wbmethods (optional): String or List. Method(s) of the WormBase Paper API which will be
                          used in the search. Default: ['title','authors','year','pmid']
    :return: Dictionary with the requested data as keys and the corresponding API returned in json format
             as values.
    """
    wbp_json_dict = {}
    if wbpmethods is None:
        wbpmethods = ['title','authors','year','pmid','doi']
    if not isinstance(wbpmethods, list):
        wbpmethods = [wbpmethods]
    for method in wbpmethods:
        wbpurl = wbpaperapiurl.format(wbpid, method)
        response = requests.get(wbpurl)
        try:
            wbp_json = response.json()
        except json.decoder.JSONDecodeError:
            print(dtnow() + ": ERROR - ERROR: json.decoder.JSONDecodeError:")
            print("Error parsing {0} for {1}.".format(wbpurl, wbpid))
            print("Please Try again.\n")
            raise ValueError
        wbp_json_dict[method] = wbp_json[method]
    return(wbp_json_dict)

def paperid2varstudydict(wbpid, wbpaperapiurl, ncbi_url, wbpmethods=None):
    """
    This function returns the results of methods title,authors,years and pmid from the WormBase Paper API
    (https://wormbase.org/about/userguide/for_developers/API-REST/Paper#0--10) for a given WormBase Paper ID
    in a format suitable for each row of the study-file which is used as an input for the  import_wormbase class
    of the PhenotypeAnnotation pipeline in the ensembl-variation repository
    (https://github.com/Ensembl/ensembl-variation/tree/release/104/modules/Bio/EnsEMBL/Variation/Pipeline/PhenotypeAnnotation).

    Args:
        :param wbpid (required): String. WormBase Paper ID.
        :param wbmethods (optional): String or List. Method(s) of the WormBase Paper API which will be
                          used in the search. Default: ['title','authors','year','pmid']

    :return: Dictionary with the requested data as keys and the corresponding API returned in json format
             as values.
    """
    wb_study_dict = {}
    wb_json_dict = wbpaperurl2json(wbpid, wbpaperapiurl)
    try:
        wb_study_dict['title'] = wb_json_dict['title']['data']
        wb_study_dict['name'] = "{0} et al., {1}".format(wb_json_dict['authors']['data'][0]['label'], wb_json_dict['year']['data'])
        wb_study_dict['pmid'] = wb_json_dict['pmid']['data']
        if wb_study_dict['pmid'] == None:
            wb_study_dict['pmid'] = papertitle2pmidncbi(wb_study_dict['title'], ncbi_url)
        print(dtnow() + ": PROGRESS - Got PMID:{0} for Paper_ID:{1}".format(wb_study_dict['pmid'], wbpid))
    except KeyError:
        print(dtnow() + ": ERROR - Some key is not in the API object for {0}\n".format(wb_json_dict))
        raise KeyError
    return(wb_study_dict)

def papertitle2pmidncbi(paper_title, ncbi_url):
    """
    Function that takes a string with a paper title as its argument and an ncbi title API url as its inputs and returns a string
    with the PMID for the input paper title.
    :param ncbi_url: String. String with the NCBI API URL.
    :param paper_tile: String. String with a paper title.
    """
    paper_ncbi_url = ncbi_url.format(paper_title)
    response = requests.get(paper_ncbi_url)
    ncbi_paper_page = requests.get(paper_ncbi_url).content
    try:
        pmid = re.findall('<Id>(\d{7,8})<\/Id>',str(ncbi_paper_page))[0]
        return (pmid)
    except IndexError:
        print(dtnow() + ": WARNING - Could not find PMID for this paper title: {0} on NCBI API URL: {1}.\n".format(paper_title, ncbi_url))
        return(None)

def species2taxonidncbi(species, ncbi_url):
    """
    Function that takes a string with a WBPS species_name_bioproject as its argument and
    an ncbi taxonomy API url as its inputs and returns a string
    with the NCBI TAXON ID for the input paper title.
    :param ncbi_url: String. String with the NCBI API URL.
    :param ncbi_url: String. String with a species name (e.g. caenorhabditis_elegans_prjna13758).
    """
    ncbi_species = " ".join(species.split('_')[0:2])
    taxon_ncbi_url = ncbi_url.format(ncbi_species)
    ncbi_taxon_page = requests.get(taxon_ncbi_url).content
    try:
        taxonid = re.findall('<Id>(\d+)<\/Id>',str(ncbi_taxon_page))[0]
    except IndexError:
        print(dtnow() + ": ERROR - Could not find taxon ID for this species: {0} on NCBI API URL: {1}.\n".format(ncbi_species, taxon_ncbi_url))
        raise IndexError
    return(taxonid)

def ifpmidnull_runncbiapi(x, ncbi_url):
    """Function that runs papertitle2pmidncbi on a dataframe if pmid is None"""
    if x['pmid']==None:
        return(papertitle2pmidncbi(x['title'],ncbi_url))
    else:
        return(x['pmid'])

def check_source_id_format(wbsourceid):
    """Function that takes a WormBase Paper ID and ensures it has the
    correct WBPaperXXXXXXXXX format."""
    if re.search('wbpaper', wbsourceid, re.IGNORECASE):
        match_nu = re.findall('\d+', wbsourceid)
        if len(match_nu) != 1 or len(match_nu[0])!=8:
            print(dtnow() + ": ERROR - WBPaper id {0} is not following the format WBPaperXXXXXXXX.\n".format(wbsourceid))
            raise ValueError
        source_numid = match_nu[0]
        correct_source_id = 'WBPaper' + source_numid
    elif re.search('wbvar', wbsourceid, re.IGNORECASE):
        match_nu = re.findall('\d+', wbsourceid)
        if len(match_nu) != 1 or len(match_nu[0])!=8:
            print(dtnow() + ": ERROR - WBVar id {0} is not following the format WBVarXXXXXXXX.\n".format(wbsourceid))
            raise ValueError
        source_numid = match_nu[0]
        correct_source_id = 'WBVar' + source_numid
    else:
        correct_source_id = wbsourceid
    return(correct_source_id)

def add_nopheno(x):
    if x['pheno_type']=='NOT':
        return(x['pheno_name'] + ', no phenotype observed')
    else:
        return(x['pheno_name'])

def fix_var_id(x, VAR_URL, RNAI_URL):
    listofvars = str(x).split('|')
    listofnewvars = []
    for oldvar in listofvars:
        newvar = oldvar.replace('WB:', '')
        if newvar.startswith('WBVar'):
            newvar = '' + newvar
        elif newvar.startswith('WBRNAi'):
            newvar = '' + newvar
        listofnewvars.append(newvar)
    return(listofnewvars)

def decide_source(daf):
    if daf['source_id'].startswith('WBVar'):
        srctype = 'Var'
    elif daf['source_id'].startswith('WBPaper'):
        if daf['new_var_id'].startswith('WBVar'):
            srctype = 'Var'
        if daf['new_var_id'].startswith('WBRNAi'):
            srctype = 'RNAi'
    else:
        srctype = None
    return(srctype)

def decide_source_id(daf):
    if daf['source_id'].startswith('WBVar'):
        srcid = daf['source_id']
    elif daf['source_id'].startswith('WBPaper'):
        srcid = daf['new_var_id']
    else:
        srcid = None
    return(srcid)

def get_source_url(daf, VAR_URL, RNAI_URL):
    if daf['source_type']=='Var':
        newurl = VAR_URL + daf['final_source_id']
    elif daf['source_type']=='RNAi':
        newurl = RNAI_URL + daf['final_source_id']
    else:
        newurl = None
    return(newurl)