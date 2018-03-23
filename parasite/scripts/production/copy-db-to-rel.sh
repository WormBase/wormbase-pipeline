#!/usr/bin/bash
diff_and_assert(){
  if [ $(diff_from_schema "$@" | wc -l) -gt 0 ]; then
    echo "WARNING - tables not the same: $@"
    return 1
  else
    return 0
  fi
}
diff_from_schema(){
 server_1=$1
 server_2=$2
 core_db=$3
 diff \
        <($(which $server_1) -e "select TABLE_NAME,TABLE_ROWS from information_schema.TABLES where TABLE_SCHEMA=\"$core_db\" order by TABLE_NAME") \
        <($(which $server_2) -e "select TABLE_NAME,TABLE_ROWS from information_schema.TABLES where TABLE_SCHEMA=\"$core_db\" order by TABLE_NAME")
}

DIR=/nfs/nobackup/ensemblgenomes/wormbase/parasite/parasite-release/WBPS${PARASITE_VERSION}
mkdir -p $DIR

stage_db() {
  DB=$1
  f=$DIR/$DB.sql
  echo $(date) " Dumping $DB to $f"
  $PARASITE_STAGING_MYSQL mysqldump $DB > $f
  echo $(date) "Loading $DB.sql to $DB mysql-ps-rel, mysql-ps-rest-rel, mysql-ps-intrel"
  mysql-ps-rel-ensrw $DB < $f &
  mysql-ps-rest-rel-ensrw $DB < $f &
  mysql-ps-intrel-ensrw $DB < $f &
  wait
  echo $(date) "Load done. Checking databases for integrity and optimizing"
  mysql-ps-rest-rel-ensrw mysqlcheck --silent --optimize $DB
  diff_and_assert $PARASITE_STAGING_MYSQL mysql-ps-rest-rel-ensrw $DB
  mysql-ps-rel-ensrw mysqlcheck --silent --optimize $DB
  diff_and_assert $PARASITE_STAGING_MYSQL mysql-ps-rel-ensrw $DB
  mysql-ps-intrel-ensrw mysqlcheck --silent --optimize $DB
  diff_and_assert $PARASITE_STAGING_MYSQL mysql-ps-intrel-ensrw $DB
}


$PARASITE_STAGING_MYSQL -NB -e 'SHOW DATABASES LIKE "%_core%_'${ENSEMBL_VERSION}'_%"' | while read -r core_db; do
  stage_db $core_db
done
stage_db ensembl_compara_parasite_${PARASITE_VERSION}_${ENSEMBL_VERSION}
stage_db ensembl_ontology_${ENSEMBL_VERSION}
stage_db ensembl_website_${ENSEMBL_VERSION}
stage_db ensemblgenomes_info_${EG_VERSION}
stage_db ensemblgenomes_stable_ids_${PARASITE_VERSION}_${ENSEMBL_VERSION}

# Compress the archive then push to the EBI Archive Freezer
echo $(date) "Creating release archive for the EBI Freezer"
tar -zcvf $DIR/release-${PARASITE_VERSION}.tar.gz $DIR/*.sql && ear-put $DIR/release-${PARASITE_VERSION}.tar.gz && rm $DIR/*.sql
