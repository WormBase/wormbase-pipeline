#!/hps/software/users/wormbase/parasite/shared/.pyenv/versions/production-tools/bin/python

import re
import sys
import os
import datetime
from ensembl.database import DBConnection, UnitTestDB
from id_to_external_db_name import idregex_to_db
from sqlalchemy import insert

def dtnow():
    ppdtnow = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return(str(ppdtnow))

def exit_with_error(error_message):
    print(dtnow() + ": ERROR - "+error_message+".")
    raise ValueError

def id_to_db(id):
    regcounter = 0
    for regx in idregex_to_db:
        regcounter += 1
        if re.match(regx, id):
            return idregex_to_db[regx]
            break
        elif regcounter==3:
            exit_with_error("The ID: "+id+" is not following one of the regex in the id_to_external_db_name.py file."+str(idregex_to_db))

#Check Input
if len(sys.argv) != 2:
    print(dtnow()+": ERROR - Usage: python add_external_db_name.py <DATABASE_CONNECTION_URL>.")
    print("\tExample: python add_external_db_name.py "
          "mysql://ensro@mysql-ps-staging-1.ebi.ac.uk:4451/caenorhabditis_elegans_prjna13758_variation_17_105_282")
    raise ValueError

DB_URL = sys.argv[1]

dbc = DBConnection(DB_URL)

#Get all protein feature attrib IDs with external ID
GET_PF_SQL="SELECT phenotype_feature_id,value FROM phenotype_feature_attrib LEFT JOIN attrib_type USING(attrib_type_id) WHERE code='external_id'"
get_pf_sql_res = dbc.execute(GET_PF_SQL)
pf_count = len(dbc.execute(GET_PF_SQL).all())
pf_to_id = {pf_id[0]:pf_id[1] for pf_id in get_pf_sql_res}
if len(pf_to_id)!=pf_count:
    exit_with_error("The phenotype_feature_id -> value mapping of the phenotype_feature_attrib "
                    "is not 1-1. SQL Statement: "+GET_PF_SQL)

#Get the external_db_name attrib_id:
GET_EDB_ID_SQL="SELECT attrib_type_id FROM attrib_type WHERE code='external_db';"
get_ebd_id_res = dbc.execute(GET_EDB_ID_SQL)
ebd_id=[x[0] for x in get_ebd_id_res][0]


#Add the external db name in the phenotype_attrib table
for pfid in pf_to_id:
    db_name = id_to_db(pf_to_id[pfid])
    stmt = insert(dbc.tables['phenotype_feature_attrib']).values(phenotype_feature_id=pfid,
                                                                 attrib_type_id=ebd_id,
                                                                 value=db_name)
    dbc.execute(stmt)
    #print(stmt)





