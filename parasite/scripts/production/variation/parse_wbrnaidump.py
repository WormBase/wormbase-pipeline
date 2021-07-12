#!/hps/software/users/wormbase/parasite/shared/.pyenv/versions/p395/bin/python
from __future__ import print_function
from wbrnaitools import *
import pandas as pd
import urllib
import re
import os
import requests
import traceback
import sys
import datetime
from intermine.webservice import Service

#Check Input
if len(sys.argv) != 5:
    print(dtnow()+": ERROR - Usage: python parse_wbdump_rnai.py SPECIES WORMBASE_VERSION"
                  "WORMBASE_DATABASES PHENOTYPE_HOME")
    print("\tExample: ./parse_wbdump_rnai.py caenorhabditis_elegans_prjna13758 $WORMBASE_VERSION "
          "/nfs/production/flicek/wormbase/wb/DATABASES $PHENOTYPE_HOME/WORMBASE/WBPS$PARASITE_VERSION")
    raise ValueError

SPECIES = sys.argv[1]
WORMBASE_VERSION = sys.argv[2]
PARASITE_VERSION = sys.argv[3]
PHENOTYPE_HOME = sys.argv[4]

SPECIES = "caenorhabditis_elegans_prjna13758"
WORMBASE_VERSION = os.getenv('WORMBASE_VERSION')
PARASITE_VERSION = os.getenv('PARASITE_VERSION')
PHENOTYPE_HOME = os.getenv('PHENOTYPE_HOME')
WORKDIR = '{0}/WBPS{1}/WORMBASE_PHENOTYPE/{2}'.format(PHENOTYPE_HOME,PARASITE_VERSION,SPECIES)
if os.path.exists(WORKDIR) is False:
    print(dtnow() + ': INFO - {0} could not be found. The directory will be created..'.format(WORKDIR))
    os.makedirs(WORKDIR)

#INPUT
WPAXREF_URL = 'http://tazendra.caltech.edu/~azurebrd/cgi-bin/forms/generic.cgi?action=WpaXref'
NCBI_PMID_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&field=title&term="{}"'
NCBI_TAXON_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term="{}"'
WB_FTP_PHENO_FILE = 'http://ftp.ebi.ac.uk/pub/databases/wormbase/releases/WS{0}/ONTOLOGY/phenotype_association.WS{0}.wb'.format(WORMBASE_VERSION)
#WB_FTP_PHENO_FILE = '~/phenotype_association.WS279.wb'
WB_FTP_RNAi_FILE = 'http://ftp.ebi.ac.uk/pub/databases/wormbase/releases/WS{0}/ONTOLOGY/phenotype_association.WS{0}.wb'.format(WORMBASE_VERSION)
WB_PAPER_API_URL = 'http://api.wormbase.org/rest/field/paper/{0}/{1}'
WB_WORMMINE_URL = 'http://intermine.wormbase.org/tools/wormmine/service'
WB_VAR_URL = '/species/c_elegans/variation/'
WB_RNAi_URL = '/species/c_elegans/rnai/'

#OUTPUT
PHENOTYPES_FILE = '{0}/wb_phenotypes.txt'.format(WORKDIR)
PHENOTYPES_FILE_VAR = '{0}/var_phenotypes.txt'.format(WORKDIR)
PHENOTYPES_FILE_RNAI = '{0}/rnai_phenotypes.txt'.format(WORKDIR)

print(dtnow() + ': INFO - Arguments given to parse_wbrnaidump.py:')
print('SPECIES = '+SPECIES)
print('WORMBASE_VERSION = '+WORMBASE_VERSION)
print('WORKDIR = '+WORKDIR)
print('WPAXREF_URL = '+WPAXREF_URL)
print('WORMBASE_FILE = '+WB_FTP_PHENO_FILE)
print('PHENOTYPES_FILE = '+PHENOTYPES_FILE)

print(dtnow() + ': INFO - Parsing Input Files')
#Read/Parse and process the WormBase phenotypes FTP file
wbpheno_df = wbphenoftp_parser(file=WB_FTP_PHENO_FILE,
                               taxon_id=species2taxonidncbi(SPECIES,NCBI_TAXON_URL),
                               VAR_URL = WB_VAR_URL,
                               RNAI_URL = WB_RNAi_URL,
                               names=None)

#Get PMIDs for the WBPapers or WBVar
print(dtnow() + ': INFO - Create a dataframe with study information.')
unique_paper_ids = list(wbpheno_df['source_id'].unique())
study_file_df = sourceids2studyfiledf(unique_paper_ids, wbpaperapiurl=WB_PAPER_API_URL, ncbi_url=NCBI_PMID_URL)

#Print WBPaper IDs without a PMID:
if study_file_df['pmid'].isnull().sum()==0:
    pass
else:
    failed_wbps = ', '.join(study_file_df[study_file_df['pmid'].isnull()]['wbp'])
    print(dtnow() + ": WARNING - Cannot get PMID neither from WB API nor from NCBI API for papers: {0}\n".format(failed_wbps))


#Get phenotype descriptions:
print(dtnow() + ': INFO - Get phenotype descriptions from WormMine.')
pheno_df = wormmine_pheno_parser(WB_WORMMINE_URL)
pheno_df.columns = ['pheno_id','pheno_name','pheno_desc']


print(dtnow() + ': INFO - Create the final phenotype association dataframe '
                'with the study and phenotype details.')


print(dtnow() + ': INFO - Adding study information to the final dataframe.')
merged_df_studies = pd.merge(wbpheno_df, study_file_df[['wbp','pmid','name','title']], left_on=wbpheno_df['source_id'], right_on=study_file_df['wbp'])
if len(merged_df_studies) == len(wbpheno_df): #check if all IDs matched
    merged_df_studies = merged_df_studies.drop('key_0', axis=1)
    merged_df_studies = merged_df_studies.drop_duplicates()
    if merged_df_studies[merged_df_studies.columns.difference(['pheno_type', 'var_id', 'pmid'])].isnull().values.any():
        print(dtnow() + ": ERROR - There are NULL values in the study merged df. NULL columns sum:")
        print(merged_df_studies[merged_df_studies.columns.difference(['pheno_type','var_id','pmid'])].isnull().sum())
        raise ValueError
else:
    print(dtnow() + ": ERROR - FORMAT ERROR: Phenotype IDs of the WB phenotype FTP file {0} "
                    "do not match 100% with the IDs parsed from : {1}.\n".format(WB_PHENO_WORMMINE_URL, WB_PAPER_API_URL))
    raise ValueError

print(dtnow() + ': INFO - Adding phenotype information to the final dataframe.')
merged_df_phenos = pd.merge(merged_df_studies, pheno_df[['pheno_id','pheno_name']], left_on=merged_df_studies['pheno_id'], right_on=pheno_df['pheno_id'])
if len(merged_df_phenos) == len(wbpheno_df): #check if all IDs matched
    merged_df_phenos = merged_df_phenos.drop('key_0', axis=1)
    merged_df_phenos = merged_df_phenos.drop_duplicates()
    if merged_df_phenos[merged_df_phenos.columns.difference(['pheno_type', 'var_id', 'pmid'])].isnull().values.any():
        print(dtnow() + ": ERROR - There are NULL values in the phenotypes merged df. NULL columns sum:")
        print(merged_df_phenos[merged_df_phenos.columns.difference(['pheno_type','var_id','pmid'])].isnull().sum())
        print('Failed IDs: ')
        print(merged_df_phenos[merged_df_phenos[merged_df_phenos.columns.difference(['pheno_type','var_id'])].isnull().any(axis=1)][['gene_id','pheno_id','source_id']].drop_duplicates())
        raise ValueError
else:
    print(dtnow() + ": ERROR - FORMAT ERROR: Source IDs of the WB phenotype FTP file {0} "
                    "do not match 100% with the IDs parsed from : {1}.\n".format(WB_FTP_PHENO_FILE, WB_PAPER_API_URL))
    raise ValueError

print(dtnow() + ': INFO - Processing final dataframe.')
merged_df_phenos['pheno_name'] = merged_df_phenos.apply(lambda x: add_nopheno(x), axis=1)
merged_df_phenos_to_write = merged_df_phenos[['gene_id', 'gene_symbol', 'pheno_name', 'pheno_id_y', 'final_source_id',
                                              'object_type', 'source_type', 'wbp', 'pmid', 'name', 'title']]
merged_df_phenos_to_write['final_source_url'] = merged_df_phenos_to_write.apply(lambda x: get_source_url(x,
                                                                                                       VAR_URL=WB_VAR_URL,
                                                                                                         RNAI_URL=WB_RNAi_URL), axis=1)
merged_df_phenos_to_write = merged_df_phenos_to_write.drop_duplicates()

print(dtnow() + ': INFO - Write output to {0} and {1}.'.format(PHENOTYPES_FILE_VAR, PHENOTYPES_FILE_RNAI))
#Write df to files
merged_df_phenos_to_write.to_csv(PHENOTYPES_FILE,sep="\t",header=None, index=False,na_rep='NULL')


# paper_id_consistency = wbpheno_df[wbpheno_df['source_id'].str.contains("(?i)wbpaper")]['source_id'].str.lower().isin(wpaxref_df['wbp'].str.lower())
#
# if paper_id_consistency.all() != true:
#     print(dtnow() + ": error - format error: rnai dump phenotype file format error.\n"
#           "reason: the wb paper ids loaded from file: {0} "
#           "do not seem to be consistent with the wb paper ids in the"
#           "wb paper mapping file here: {1}.\n"
#           "the failing ids are:".format(rnai_file, wpaxref_url))
#     print(wbpheno_df[wbpheno_df['source_id'].str.contains("(?i)wbpaper")][~paper_id_consistency]['source_id'])
#     raise valueerror

#Change the paper IDs from the WB RNAi dump file with their corresponding IDs from the
#WB Paper Mapping file to make sure the IDs are correct (there are some cases which
#paper id is WBpaperXXXXXXXX instead of WBPaperXXXXXXXX) and we don't want that.
# merged_df = pd.merge(wbrnai_df, wpaxref_df['wbp'], left_on=wbrnai_df["papers"].str.lower(), right_on=wpaxref_df['wbp'].str.lower())
# if len(merged_df) == len(wbrnai_df): #check if all IDs matched
#     wbrnai_df = merged_df.drop('key_0', axis=1)
# else:
#     print(dtnow() + ": ERROR - FORMAT ERROR: RNAi Dump Phenotype file FORMAT ERROR.\n"
#           "Reason: The WB Paper IDs loaded from file: {0} "
#           "do not seem to be consistent with the WB Paper IDs in the"
#           "WB Paper mapping file here: {1}.\n"
#           "The Failing IDs are:".format(RNAI_FILE, WPAXREF_URL))
#     print(wbrnai_df[~wbrnai_df['papers'].isin(merged_df['wbp'])])
#     raise ValueError