#!/usr/bin/bash
# Takes the core databases from previous staging server, and gets the current staging server to the "start of new release" situation
# - copies core databases from the previous release, upgrading them to the new schema
# - pinches new comparator DBs from known Ensembl places

perl -MProductionMysql -E '
   say for ProductionMysql->previous_staging->core_databases(@ARGV ? @ARGV : "core_$ENV{PREVIOUS_PARASITE_VERSION}");
' "$@" | while read -r DB; do
  echo "Copying $DB"
  NEWDB=$( sed "s/core_$PREVIOUS_PARASITE_VERSION\_$PREVIOUS_ENSEMBL_VERSION/core_$PARASITE_VERSION\_$ENSEMBL_VERSION/" <<< $DB )
  echo "Looking for previous versions of $NEWDB"
  DB_PAT=$(echo $NEWDB | sed "s/core_.*/core%/")
  ${PARASITE_STAGING_MYSQL} -Ne "show databases like \"$DB_PAT\"" | grep -v "$NEWDB" | while read -r OLD_DB ; do
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
  echo "Using schema_patcher.pl to patch $NEWDB"
  perl $ENSEMBL_CVS_ROOT_DIR/ensembl/misc-scripts/schema_patcher.pl $(${PARASITE_STAGING_MYSQL}-ensrw details script) \
     --release $ENSEMBL_VERSION --nointeractive \
     --database "$NEWDB"
done

if [ $# -gt 0 ] ; then
  echo "Ran in partial mode. Will not drop old stuff or copy comparators"
  exit
fi
find_host_and_db(){
  species="$1"
  # if mysql-ens-sta-[1234] don't have the comparator dbs required as the ensembl version is too old
  # then mysql-ens-mirror-[1234] will hopefully provide the db
  # *but* check it's OK to dump from  the mirror servers before using them
  for cmd in mysql-ens-sta-1 mysql-ens-sta-2 mysql-ens-sta-3 mysql-ens-sta-4; do
    db=$( $cmd -Ne "show databases like \"%${species}_core%_${ENSEMBL_VERSION}_%\"" )
    if [ "$db" ]; then
      echo $cmd$'\t'$db
      return
    fi
  done
}

copy_comparator(){
  res=$(find_host_and_db "$@" )
  if [ ! "$res" ]; then
    echo "Could not find comparator: $@"
    exit 1
  fi
  read host db <<< $res
  echo "Copying: $host $db"
  $PARASITE_STAGING_MYSQL-ensrw -e "create database if not exists $db"
  $host mysqldump $db | $PARASITE_STAGING_MYSQL-ensrw $db
}

echo "Dropping old comparators and databases that aren't for this WPBS version"
$PARASITE_STAGING_MYSQL -Ne "show databases like \"%_core_%\"" | grep -v core_${PARASITE_VERSION} | while read -r core_db; do
  echo "Dropping: $core_db"
  ${PARASITE_STAGING_MYSQL}-ensrw -e "drop database $core_db"
done

$PREVIOUS_PARASITE_STAGING_MYSQL -Ne "show databases like \"%_core_%\"" | grep -v core_${PREVIOUS_PARASITE_VERSION} | perl -pe 's/_core.*//' \
  | sort | while read -r species; do
  copy_comparator $species
done
