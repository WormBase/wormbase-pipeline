#!/usr/bin/bash

for VAR in  HANDOVER_STAGING_MYSQL
do
   if [ -z $(eval "echo \$$VAR") ]; then
      echo "$VAR must be set in your environment before running $0"
      exit 254
   fi
done;

# Make sure the dbs have the correct species.production name in the db
for f in $($HANDOVER_STAGING_MYSQL -Ne "SHOW DATABASES LIKE \"%${EG_VERSION}_${ENSEMBL_VERSION}_%\";");

  do new_s_name=$(echo $f | cut -d'_' -f1,2);

  $HANDOVER_STAGING_MYSQL-w $f -Ne "UPDATE meta SET meta_value = '${new_s_name}'
  WHERE meta_key = 'species.production_name';";

  # Make sure the dbs have the correct species.division name in the db
  $HANDOVER_STAGING_MYSQL-w $f -Ne "UPDATE meta SET meta_value = 'EnsemblMetazoa'
  WHERE meta_key = 'species.division';";

  # Remove any repeat elements with repeat_start < 0 from the repeat features table
  $HANDOVER_STAGING_MYSQL-w $f -e "DELETE FROM repeat_feature WHERE repeat_start < 1;";

  # Remove gene versions:
  $HANDOVER_STAGING_MYSQL-w $f -Ne "UPDATE gene INNER JOIN
  seq_region sr USING (seq_region_id) INNER JOIN
  coord_system cs USING (coord_system_id)
  SET gene.version = NULL
  WHERE cs.species_id = 1
  AND gene.version IS NOT NULL;";

  #Correct meta values regarding assembly and annotation providers names and urls:
  $HANDOVER_STAGING_MYSQL-w $f -Ne "UPDATE meta SET meta_key='assembly.provider_name' WHERE meta_key='provider.name';";

  $HANDOVER_STAGING_MYSQL-w $f -Ne "UPDATE meta SET meta_key='assembly.provider_url' WHERE meta_key='provider.url';";

  $HANDOVER_STAGING_MYSQL-w $f -Ne "INSERT IGNORE INTO meta(species_id,meta_key,meta_value)
  VALUES ('1','annotation.provider_name',(SELECT meta_value FROM (SELECT meta_value FROM meta WHERE meta_key = 'assembly.provider_name' LIMIT 1) t));";

  $HANDOVER_STAGING_MYSQL-w $f -Ne "INSERT IGNORE INTO meta(species_id,meta_key,meta_value)
  VALUES ('1','annotation.provider_url',(SELECT meta_value FROM (SELECT meta_value FROM meta WHERE meta_key = 'assembly.provider_url' LIMIT 1) t));";

  #Insert genebuild.method if not exist in meta table
  $HANDOVER_STAGING_MYSQL-w $f -Ne "INSERT IGNORE INTO meta(species_id,meta_key,meta_value)
  VALUES ('1','genebuild.method','import');";

  #Drop old tables which do not belong to the schema
  $HANDOVER_STAGING_MYSQL-w $f -Ne "DROP TABLE IF EXISTS input_id_analysis;";
  $HANDOVER_STAGING_MYSQL-w $f -Ne "DROP TABLE IF EXISTS analysis_description_bak;";
  $HANDOVER_STAGING_MYSQL-w $f -Ne "DROP TABLE IF EXISTS input_id_type_analysis;";
  $HANDOVER_STAGING_MYSQL-w $f -Ne "DROP TABLE IF EXISTS rule_conditions;";
  $HANDOVER_STAGING_MYSQL-w $f -Ne "DROP TABLE IF EXISTS rule_goal;";

  trsamp=$($HANDOVER_STAGING_MYSQL $core_db -Ne \
  "SELECT meta_value FROM meta WHERE meta_key='sample.transcript_param';");
  bestsamp=$($HANDOVER_STAGING_MYSQL $core_db -Ne \
  "SELECT stable_id FROM transcript WHERE stable_id LIKE \"%${trsamp}%\" LIMIT 1;");
  echo $core_db $trsamp $bestsamp;
  $HANDOVER_STAGING_MYSQL-w $core_db -Ne "UPDATE meta SET meta_value='${bestsamp}' WHERE meta_key='sample.transcript_param';";

  $HANDOVER_STAGING_MYSQL-w $core_db -Ne \
  "UPDATE transcript t
  INNER JOIN xref x on t.display_xref_id = x.xref_id
  INNER JOIN object_xref USING (xref_id)
  INNER JOIN seq_region USING (seq_region_id)
  INNER JOIN coord_system USING (coord_system_id)
  SET t.display_xref_id = NULL
  WHERE x.dbprimary_acc regexp '-20[[:digit:]]$' AND
  x.dbprimary_acc = x.display_label AND
  ensembl_object_type = 'Transcript' AND
  species_id = 1;";
done
