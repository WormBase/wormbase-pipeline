#!/hps/software/users/wormbase/parasite/shared/.pyenv/versions/p395/bin/python
from __future__ import print_function
from helpers import *
import pandas as pd
import urllib
import re
import os
import requests
import traceback
import sys
import datetime
from intermine.webservice import Service

#USAGE
def usage():
    return('USAGE: python cele_pheno_datachecks phenotype_file.txt')


if len(sys.argv) > 1:
    PHENOTYPES_FILE = sys.argv[1]
else:
    print(usage())
    raise ValueError

#INPUT
SPECIES = "caenorhabditis_elegans_prjna13758"
WORMBASE_VERSION = os.getenv('WORMBASE_VERSION')
PARASITE_VERSION = os.getenv('PARASITE_VERSION')
PHENOTYPE_HOME = os.getenv('PHENOTYPE_HOME')
WORKDIR = '{0}/WBPS{1}/WORMBASE_PHENOTYPE/{2}'.format(PHENOTYPE_HOME,PARASITE_VERSION,SPECIES)

wbpheno_df_columns = {'gene_id':0,
                      'gene_symbol':1,
                      'pheno_type':2,
                      'pheno_id':3,
                      'var_id':4,
                      'object_type':5,
                      'source_type':6,
                      'source_wb_id':7,
                      'pmid':8,
                      'author':9,
                      'phenotype_description':10,
                      'source_url':11
                      }

phefile = pd.read_csv(PHENOTYPES_FILE, header=None, sep='\t')
print(dtnow() + ': INFO - Loaded '+PHENOTYPES_FILE)

print(dtnow() + ': INFO - Starting '+PHENOTYPES_FILE+' datachecks.')
error_messages = [z for z in [final_file_datacheck(phefile, x, wbpheno_df_columns[x]) for x in wbpheno_df_columns] if z != ""]

if len(error_messages) > 0:
    print(dtnow() + ': ERROR - Datachecks failed:')
    print('\n'.join(error_messages))
    print(
        'Check whether the ' + file + ' format has been changed and adapt to these changes but tweaking the wbpheno_df_columns in the main script. Exiting.')
    raise ValueError
print(dtnow() + ': INFO - Datachecks passed!.')


