
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
mkdir /nfs/nobackup/ensemblgenomes/wormbase/parasite/parasite-release
cd /nfs/nobackup/ensemblgenomes/wormbase/parasite/parasite-release

# Dump all the databases into the tmp space
echo "Dumping and loading databases"
databaselist=$($PARASITE_STAGING_MYSQL -NB -e 'SHOW DATABASES LIKE "%_core%_'${ENSEMBL_VERSION}'_%"')
databases=( `echo ${databaselist}` ensembl_compara_parasite_${PARASITE_VERSION}_${ENSEMBL_VERSION} ensembl_ontology_${ENSEMBL_VERSION} ensembl_website_${ENSEMBL_VERSION} ensemblgenomes_info_${EG_VERSION} ensemblgenomes_stable_ids_${PARASITE_VERSION}_${ENSEMBL_VERSION} )
stagingurl=$(echo $($PARASITE_STAGING_MYSQL details url) | awk '{ sub(/^mysql:\/\//,""); print }' | awk '{ sub(/\/$/,""); print }')
for DB in "${databases[@]}"
do
  echo $DB
  $PARASITE_STAGING_MYSQL mysqldump $DB > $DB.sql
  echo "Loading to mysql-ps-rel"
  mysql-ps-rel-ensrw -e "CREATE DATABASE $DB"
  mysql-ps-rel-ensrw $DB < $DB.sql
  mysql-ps-rel-ensrw mysqlcheck -o $DB
  diff_and_assert $PARASITE_STAGING_MYSQL mysql-ps-rel-ensrw $DB
  echo "Loading to mysql-ps-rest-rel"
  mysql-ps-rest-rel-ensrw -e "CREATE DATABASE $DB"
  mysql-ps-rest-rel-ensrw $DB < $DB.sql
  mysql-ps-rest-rel-ensrw mysqlcheck -o $DB
  diff_and_assert $PARASITE_STAGING_MYSQL mysql-ps-rest-rel-ensrw $DB
  echo "Loading to mysql-ps-intrel"
  mysql-ps-intrel-ensrw -e "CREATE DATABASE $DB"
  mysql-ps-intrel-ensrw $DB < $DB.sql
  mysql-ps-intrel-ensrw mysqlcheck -o $DB
  diff_and_assert $PARASITE_STAGING_MYSQL mysql-ps-intrel-ensrw $DB
done


# Compress the archive then push to the EBI Archive Freezer
echo "Creating release archive for the EBI Freezer"
tar -zcvf release-${PARASITE_VERSION}.tar.gz *.sql
rm *.sql
ear-put release-${PARASITE_VERSION}.tar.gz

