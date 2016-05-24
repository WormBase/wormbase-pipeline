
## Script to copy core databases from previous staging site to current
## After copy, relevant patches are applied to bring databases up to the relevant Ensembl schema
## Before running, load the correct module as this creates the environment variables, e.g. "module load parasite_prod_rel6"

echo "Begin database copying"

databaselist=$($PREVIOUS_PARASITE_STAGING_MYSQL -NB -e 'SHOW DATABASES LIKE "%\_core\_%"')
databases=( `echo ${databaselist}` )

for DB in "${databases[@]}"
do
  echo "Copying $DB"
  NEWDB=$(echo $DB | sed "s/core_$PREVIOUS_PARASITE_VERSION\_$PREVIOUS_ENSEMBL_VERSION/core_$PARASITE_VERSION\_$ENSEMBL_VERSION/" | sed "s/core_$PREVIOUS_ENSEMBL_VERSION/core_$ENSEMBL_VERSION/" | sed "s/_$PREVIOUS_WORMBASE_VERSION/_$WORMBASE_VERSION/")
  echo  "Creating $NEWDB"
  ${PARASITE_STAGING_MYSQL}-ensrw -e "CREATE DATABASE $NEWDB"
  echo "Dumping $DB to $NEWDB"
  ${PREVIOUS_PARASITE_STAGING_MYSQL}-ensrw mysqldump $DB | ${PARASITE_STAGING_MYSQL}-ensrw $NEWDB
done

echo "Using schema_patcher.pl to patch databases"

perl $ENSEMBL_CVS_ROOT_DIR/ensembl/misc-scripts/schema_patcher.pl $(${PARASITE_STAGING_MYSQL}-ensrw details script) --release $ENSEMBL_VERSION --type core --verbose --nointeractive

