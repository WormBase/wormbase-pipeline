#!/usr/bin/bash

RELEASE=$1
EMAIL=$2
if [[ -z $RELEASE || -z $EMAIL ]]; then
   echo "Usage: $0 <current_release> <email>"
   echo "       current_release is passed to metadata_updater.pl as -current_release"
   echo "       email is passed to metadata_updater.pl as -email"
   exit 255
fi

for VAR in  PARASITE_STAGING_MYSQ   \
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

db_name='ensembl_metadata'
datestamp=$(date +%Y-%m-%d)
timestamp=$(date +%Y-%m-%d_%H:%M:%S)
EM_backup_file="${PARASITE_SCRATCH}/backup/${db_name}_${timestamp}.sql"
EM_repo="${ENSEMBL_CVS_ROOT_DIR}/ensembl-metadata"
EM_table_SQL="${EM_repo}/sql/table.sql"

echo "Dumping up current ${db_name} to ${EM_backup_file}.gz"
$PARASITE_STAGING_MYSQL mysqldump ${db_name} > ${EM_backup_file} && gzip ${EM_backup_file}

echo "Dropping current ${db_name}"
${PARASITE_STAGING_MYSQL}-ensrw -N -e "drop database if exists ${db_name};"

echo "Creating new ${db_name} from ${EM_table_SQL}"
${PARASITE_STAGING_MYSQL}-ensrw -N -e "create database ${db_name};"
cat ${EM_table_SQL} | ${PARASITE_STAGING_MYSQL}-ensrw ${db_name}
for this_core_db in $(${PARASITE_STAGING_MYSQL} -N -e "SHOW DATABASES LIKE \"%core_${PARASITE_VERSION}_${ENSEMBL_VERSION}%\""); do
   echo "Updating ${db_name} with metadata from ${this_core_db}";
   metadata_updater.pl  -metadata_uri     $($PARASITE_STAGING_MYSQL-ensrw details url)${db_name}      \
                        -database_uri     $($PARASITE_STAGING_MYSQL-ensrw details url)${this_core_db} \
                        -e_release        ${ENSEMBL_VERSION}                                          \
                        -eg_release       ${EG_VERSION}                                               \
                        -parasite_release ${PARASITE_VERSION}                                         \
                        -release_date     "'${datestamp}'"                                            \
                        -current_release  "'${RELEASE}'"                                              \
                        -email            "'${EMAIL}'"                                                \
                        -comment          "'creating ${db_name} for WBPS ${PARASITE_VERSION}'"        \
                        -source           "'Manual'"
done

echo -n "Done at "
date
