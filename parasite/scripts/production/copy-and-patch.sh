
## Script to copy core databases from previous staging site to current
## After copy, relevant patches are applied to bring databases up to the relevant Ensembl schema
## Before running, load the correct module as this creates the environment variables, e.g. "module load parasite_prod_rel6"

echo "These were the comparator databases: No refresh needed? "
$PREVIOUS_PARASITE_STAGING_MYSQL -e "show databases like \"%_core_%\"" | grep -v core_${PREVIOUS_PARASITE_VERSION}
echo "Begin database copying"

databaselist=$($PREVIOUS_PARASITE_STAGING_MYSQL -NB -e 'SHOW DATABASES LIKE "%\_core\_%"')
databases=( `echo ${databaselist}` )

for DB in "${databases[@]}"
do
  echo "Copying $DB"
  NEWDB=$(echo $DB | sed "s/core_$PREVIOUS_PARASITE_VERSION\_$PREVIOUS_ENSEMBL_VERSION/core_$PARASITE_VERSION\_$ENSEMBL_VERSION/" | sed "s/_$PREVIOUS_ENSEMBL_VERSION/_$ENSEMBL_VERSION/" | sed "s/_$PREVIOUS_WORMBASE_VERSION$/_$WORMBASE_VERSION/")
  echo "Looking for previous versions of $NEWDB"
  DB_PAT=$(echo $NEWDB | sed "s/core_.*/core%/")
  ${PARASITE_STAGING_MYSQL} -e "show databases like \"$DB_PAT\"" | grep -v Database | grep -v "$NEWDB" | while read -r OLD_DB ; do
    echo "Dropping database $OLD_DB"
    ${PARASITE_STAGING_MYSQL}-ensrw -e "DROP DATABASE $OLD_DB"
  done
  echo  "Creating $NEWDB"
  if ! ${PARASITE_STAGING_MYSQL}-ensrw -e "CREATE DATABASE $NEWDB"; then
    echo "Could not create $NEWDB - should this script be running? Bailing out."
    exit 1
  fi
  echo "Dumping $DB to $NEWDB"
  ${PREVIOUS_PARASITE_STAGING_MYSQL}-ensrw mysqldump $DB | ${PARASITE_STAGING_MYSQL}-ensrw $NEWDB
done

echo "Using schema_patcher.pl to patch databases"

perl $ENSEMBL_CVS_ROOT_DIR/ensembl/misc-scripts/schema_patcher.pl $(${PARASITE_STAGING_MYSQL}-ensrw details script) --release $ENSEMBL_VERSION --type core --verbose --nointeractive

