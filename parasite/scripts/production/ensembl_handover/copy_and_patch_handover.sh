#!/usr/bin/bash
# Takes the core databases from previous staging server, and gets the current staging server to the "start of new release" situation
# - copies core databases from the previous release, upgrading them to the new schema
# USAGE
error_message="Usage: $0 handover_list.tsv [drop_old_db/no_drop]\nno_drop: The previous db will not be dropped.\ndrop_old_db: The previous db will be dropped\n."
if [ ! "$1" ] || [ ! "$2" ]; then printf "$error_message"; exit 1; fi
hanlist=$1
drop_option=$2
if [[ ! "$drop_option" =~ ^("drop_old_db"|"no_drop")$ ]]; then printf "$error_message"; exit 1; fi

while read -r DB; do

echo "Copying $DB"
NEWDB=$(echo $DB | cut -d'_' -f1,2)_$(echo $DB | rev | cut -d'_' -f4 | rev)_${EG_VERSION}_${ENSEMBL_VERSION}_$(echo $DB | rev | cut -d'_' -f1 | rev)
echo "Looking for previous versions of $NEWDB"
#DB_PAT=$( sed "s/_$EG_VERSION\_$ENSEMBL_VERSION/_$PREVIOUS_EG_VERSION\_$PREVIOUS_ENSEMBL_VERSION/" <<< $NEWDB )
DB_PAT=$(echo $NEWDB | grep 'core' |  sed "s/core_.*/core%/")
if [[ ${drop_option} == "drop_old_db" ]];
  then echo "drop_old_db option selected";
       ${HANDOVER_STAGING_MYSQL} -Ne "show databases like \"$DB_PAT\"" | grep -v "$DB_PAT" | while read -r OLD_DB ; do
       echo "Dropping database $OLD_DB"
       ${HANDOVER_FROM_STAGING_MYSQL}-ensrw -e "DROP DATABASE $OLD_DB"
       done;
elif [[ ${drop_option} == "no_drop" ]];
  then echo "no_drop option selected: Skipped dropping OLD_DB";
else printf "$error_message"; exit 1;
fi
echo  "Creating $NEWDB"
if ! ${HANDOVER_STAGING_MYSQL}-w -e "CREATE DATABASE $NEWDB"; then
  echo "Could not create $NEWDB - should this script be running? Bailing out."
  exit 1
fi
echo "Dumping $DB to $NEWDB"
${HANDOVER_FROM_STAGING_MYSQL} mysqldump $DB | ${HANDOVER_STAGING_MYSQL}-w $NEWDB
echo "Using schema_patcher.pl to patch $NEWDB"
perl $ENSEMBL_CVS_ROOT_DIR/ensembl/misc-scripts/schema_patcher.pl $(${HANDOVER_STAGING_MYSQL}-w details script) \
   --release $ENSEMBL_VERSION --nointeractive \
   --database "$NEWDB"
done < $hanlist


