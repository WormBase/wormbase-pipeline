#!/usr/bin/bash
# Takes the core databases from previous staging server, and gets the current staging server to the "start of new release" situation
# - copies core databases from the previous release, upgrading them to the new schema

perl -MProductionMysql -E '
   say for ProductionMysql->previous_staging->core_databases(@ARGV ? @ARGV : "core_$ENV{PREVIOUS_PARASITE_VERSION}_$ENV{PREVIOUS_ENSEMBL_VERSION}_");
   say for ProductionMysql->previous_staging->variation_databases(@ARGV ? @ARGV : "variation_$ENV{PREVIOUS_PARASITE_VERSION}_$ENV{PREVIOUS_ENSEMBL_VERSION}_");
' "$@" | while read -r DB; do
  echo "Copying $DB"
  NEWDB=$( sed "s/_$PREVIOUS_PARASITE_VERSION\_$PREVIOUS_ENSEMBL_VERSION/_$PARASITE_VERSION\_$ENSEMBL_VERSION/" <<< $DB )
  echo "Looking for previous versions of $NEWDB"
  DB_PAT=$(echo $NEWDB | grep 'core' |  sed "s/core_.*/core%/")
  if [ -z $DB_PAT ]; then
    DB_PAT=$(echo $NEWDB | grep 'variation' |  sed "s/variation_.*/variation%/")
  fi
  ${PARASITE_STAGING_MYSQL} -Ne "show databases like \"$DB_PAT\"" | grep -v "$NEWDB" | while read -r OLD_DB ; do
    echo "Dropping database $OLD_DB"
    ${PARASITE_STAGING_MYSQL}-w -e "DROP DATABASE $OLD_DB"
  done
  echo  "Creating $NEWDB"
  if ! ${PARASITE_STAGING_MYSQL}-w -e "CREATE DATABASE $NEWDB"; then
    echo "Could not create $NEWDB - should this script be running? Bailing out."
    #exit 1
    continue
  fi
  echo "Dumping $DB to $NEWDB"
  ${PREVIOUS_PARASITE_STAGING_MYSQL} mysqldump $DB | ${PARASITE_STAGING_MYSQL}-w $NEWDB
  echo "Using schema_patcher.pl to patch $NEWDB"
  perl $ENSEMBL_CVS_ROOT_DIR/ensembl/misc-scripts/schema_patcher.pl $(${PARASITE_STAGING_MYSQL}-w details script) \
     --release $ENSEMBL_VERSION --nointeractive \
     --database "$NEWDB"
done

