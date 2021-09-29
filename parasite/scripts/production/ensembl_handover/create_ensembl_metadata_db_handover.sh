#!/usr/bin/bash

RELEASE=$1
EMAIL=$2
if [[ -z $RELEASE || -z $EMAIL ]]; then
   echo "Usage: $0 <current_release> <email>"
   echo "       current_release is passed to metadata_updater.pl as -current_release"
   echo "       email is passed to metadata_updater.pl as -email"
   exit 255
fi

for VAR in  HANDOVER_STAGING_MYSQL  \
            PARASITE_SCRATCH        \
            PARASITE_VERSION        \
            ENSEMBL_VERSION         \
            EG_VERSION              \
            ENSEMBL_CVS_ROOT_DIR
do
   if [ -z $(eval "echo \$$VAR") ]; then
      echo "$VAR must be set in your environment before running $0"
      exit 254
   fi
done

# bale out immediately on error
set -e

db_name="ensembl_metadata_${ENSEMBL_VERSION}"
datestamp=$(date +%Y-%m-%d)
timestamp=$(date +%Y-%m-%d_%H:%M:%S)
EM_repo="${ENSEMBL_CVS_ROOT_DIR}/ensembl-metadata"
EM_table_SQL="${EM_repo}/sql/table.sql"

echo "Dropping current ${db_name}"
${HANDOVER_STAGING_MYSQL}-w -N -e "drop database if exists ${db_name};"

echo "Creating new ${db_name} from ${EM_table_SQL}"
${HANDOVER_STAGING_MYSQL}-w -N -e "create database ${db_name};"
cat ${EM_table_SQL} | ${HANDOVER_STAGING_MYSQL}-w ${db_name}
for this_core_db in $(${HANDOVER_STAGING_MYSQL} -N -e "SHOW DATABASES LIKE \"%core_${EG_VERSION}_${ENSEMBL_VERSION}%\""); do
   echo "Updating ${db_name} with metadata from ${this_core_db}";
   perl ${ENSEMBL_CVS_ROOT_DIR}/ensembl-metadata/misc_scripts/metadata_updater.pl  -metadata_uri $($HANDOVER_STAGING_MYSQL-w details url)${db_name}      \
                        -database_uri     $($HANDOVER_STAGING_MYSQL-w details url)${this_core_db} \
                        -e_release        ${ENSEMBL_VERSION}                                          \
                        -eg_release       ${EG_VERSION}                                         \
                        -release_date     "${datestamp}"                                            \
                        -current_release  "${ENSEMBL_VERSION}"                                              \
                        -email            "'${EMAIL}'"                                                \
                        -comment          "'creating ${db_name}'"        \
                        -source           "'Manual'"
done

echo -n "Done at "
date
