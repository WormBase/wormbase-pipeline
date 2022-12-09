species=$1

core_db=$($PARASITE_STAGING_MYSQL -Ne "SHOW DATABASES" | grep "${species}_core")
how_many_core_dbs=$(echo $core_db | wc -l)

if [[ $how_many_core_dbs != 1 ]]; then
  echo "Multiple core dbs. Exiting"
  exit
fi

otherfeatures_db=$(echo $core_db | perl -pe "s/_core_/_otherfeatures_/");
echo creating $otherfeatures_db;
$PARASITE_STAGING_MYSQL-w -e "drop database if exists $otherfeatures_db; create database $otherfeatures_db";
echo loading schema from $core_db to $otherfeatures_db;
$PARASITE_STAGING_MYSQL mysqldump --no-data $core_db | $PARASITE_STAGING_MYSQL-w -D $otherfeatures_db;
echo copying tables;
$PARASITE_STAGING_MYSQL mysqldump $core_db meta assembly coord_system seq_region \
                                           seq_region_attrib attrib_type external_db \
                                           misc_set unmapped_reason biotype seq_region_synonym | $PARASITE_STAGING_MYSQL-w -D $otherfeatures_db;
echo "Done"
