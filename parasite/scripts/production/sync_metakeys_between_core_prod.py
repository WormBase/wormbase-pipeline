#!/hps/software/users/wormbase/parasite/shared/.pyenv/versions/p395/bin/python
import os
import re
import sys
import datetime
import subprocess

def dtnow():
    ppdtnow = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    return(str(ppdtnow))

def listify(subprocess_sql_run):
    listified = [x for x in subprocess_sql_run.stdout.decode('utf-8').split("\n") if x != '']
    return(listified)

print(sys.argv)
if len(sys.argv) != 5:
    print(dtnow() + ": ERROR - Usage: python sync_metakeys_between_core_prod.py STAGING_SERVER PRODUCTION_SERVER"
                    "PRODUCTION_DATABASE SQL_REGEX_TO_CAPTURE_CORE_DBS")
    print("\tExample: ./sync_metakeys_between_core_prod.py mysql-ps-staging-1 mysql-ps-prod-1 "
          "ensembl_production_parasite %_core_16_101_%")
    raise ValueError

STAGING_SERVER = sys.argv[1]
PRODUCTION_SERVER = sys.argv[2]
PRODUCTION_DB = sys.argv[3]
SREGEX = sys.argv[4]

#STAGING_SERVER=os.getenv("HANDOVER_STAGING_MYSQL")
#PRODUCTION_SERVER="mysql-ps-prod-1"
PRODUCTION_SERVER_W=PRODUCTION_SERVER+"-w"
#PRODUCTION_DB="ensembl_production_handover_16_105"
ENSEMBL_VERSION=os.getenv("ENSEMBL_VERSION")
#SREGEX="%_"+ENSEMBL_VERSION+"_%"

GET_CORE_SQL="SHOW DATABASES LIKE \""+SREGEX+"\";"
CORE_SQL="SELECT DISTINCT(meta_key) FROM meta"
PROD_SQL_REQ="SELECT name FROM meta_key WHERE FIND_IN_SET('core', db_type) AND is_current = 1 AND is_optional=0"
PROD_SQL_ALL="SELECT name FROM meta_key WHERE FIND_IN_SET('core', db_type) AND is_current = 1"

core_dbs = listify(subprocess.run([STAGING_SERVER,'-Ne',GET_CORE_SQL], stdout=subprocess.PIPE))

print(dtnow() + ': INFO - Arguments given:')
print('                     STAGING_SERVER: '+STAGING_SERVER)
print('                     PRODUCTION_SERVER = '+PRODUCTION_SERVER)
print('                     PRODUCTION_DB = '+PRODUCTION_DB)
print('                     ENSEMBL_VERSION = '+ENSEMBL_VERSION)
print('                     CORE_SQL = '+GET_CORE_SQL)

for this_core in core_dbs:
    core_meta_keys = list(set(listify(subprocess.run([STAGING_SERVER,this_core,'-Ne',CORE_SQL], stdout=subprocess.PIPE))))
    all_prod_meta_keys = list(set(listify(subprocess.run([PRODUCTION_SERVER,PRODUCTION_DB,'-Ne',PROD_SQL_ALL], stdout=subprocess.PIPE))))
    req_prod_meta_keys = list(set(listify(subprocess.run([PRODUCTION_SERVER,PRODUCTION_DB,'-Ne',PROD_SQL_REQ], stdout=subprocess.PIPE))))
    #Handle meta_keys in core_dbs not present in the production db
    print(dtnow() + " - INFO: Syncing meta_keys for "+this_core+".")
    extra_core_meta = [x for x in core_meta_keys if x not in all_prod_meta_keys]
    if len(extra_core_meta)==0:
        print("                      Meta keys in core db meta table and production db meta_key table are in sync. Skipping.")
        print("\n")
        pass
    else:
        print("                      Core db meta_keys not present in the production db:")
        print("                       " + "\n                       ".join(extra_core_meta))
        PROD_INSERT_SQL = "INSERT IGNORE INTO meta_key(name, is_optional, is_current, db_type, is_multi_value) VALUES"
        PROD_INSERT_SQL += " " + ",".join(['(' + "'" + x + "'" + ",'1','1','core','0')" for x in extra_core_meta])
        print(PROD_INSERT_SQL)
        #subprocess.run([PRODUCTION_SERVER_W, PRODUCTION_DB, '-Ne', PROD_INSERT_SQL], stdout=subprocess.PIPE)
    print("\n")


print("\n\n")
print(dtnow()+ " - INFO: Done!")
quit()
