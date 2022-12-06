#!/usr/bin/bash

for VAR in  HANDOVER_STAGING_MYSQL
do
   if [ -z $(eval "echo \$$VAR") ]; then
      echo "$VAR must be set in your environment before running $0"
      exit 254
   fi
done;

# Make sure the dbs have the correct species.production name in the db
for f in $($PARASITE_STAGING_MYSQL -Ne "SHOW DATABASES LIKE \"%_core_${PARASITE_VERSION}_${ENSEMBL_VERSION}_%\";");

  do 

  # Remove any repeat elements with repeat_start < 0 from the repeat features table
  $PARASITE_STAGING_MYSQL-w $f -e "DELETE FROM repeat_feature WHERE repeat_start < 1;";

  # Remove gene versions:
  $PARASITE_STAGING_MYSQL-w $f -Ne "UPDATE gene INNER JOIN
  seq_region sr USING (seq_region_id) INNER JOIN
  coord_system cs USING (coord_system_id)
  SET gene.version = NULL
  WHERE cs.species_id = 1
  AND gene.version IS NOT NULL;";

  #Insert genebuild.method if not exist in meta table
  $PARASITE_STAGING_MYSQL-w $f -Ne "INSERT IGNORE INTO meta(species_id,meta_key,meta_value)
  VALUES ('1','genebuild.method','import');";

  #Drop old tables which do not belong to the schema
  $PARASITE_STAGING_MYSQL-w $f -Ne "DROP TABLE IF EXISTS input_id_analysis;";
  $PARASITE_STAGING_MYSQL-w $f -Ne "DROP TABLE IF EXISTS analysis_description_bak;";
  $PARASITE_STAGING_MYSQL-w $f -Ne "DROP TABLE IF EXISTS input_id_type_analysis;";
  $PARASITE_STAGING_MYSQL-w $f -Ne "DROP TABLE IF EXISTS rule_conditions;";
  $PARASITE_STAGING_MYSQL-w $f -Ne "DROP TABLE IF EXISTS rule_goal;";
  $PARASITE_STAGING_MYSQL-w $f -Ne "DROP TABLE IF EXISTS job_status;";
  $PARASITE_STAGING_MYSQL-w $f -Ne "DROP TABLE IF EXISTS job;";

  trsamp=$($PARASITE_STAGING_MYSQL $f -Ne \
  "SELECT meta_value FROM meta WHERE meta_key='sample.transcript_param';");
  bestsamp=$($PARASITE_STAGING_MYSQL $f -Ne \
  "SELECT stable_id FROM transcript WHERE stable_id LIKE \"%${trsamp}%\" LIMIT 1;");

  $PARASITE_STAGING_MYSQL-w $f -Ne "UPDATE meta SET meta_value='${bestsamp}' WHERE meta_key='sample.transcript_param';";

  gsd=$($PARASITE_STAGING_MYSQL-w $f -Ne "SELECT meta_value FROM meta WHERE meta_key='genebuild.start_date';");
  gv=$($PARASITE_STAGING_MYSQL-w $f -Ne "SELECT meta_value FROM meta WHERE meta_key='genebuild.version';");

  if [[ "$gsd" == "$gv" ]]; then $PARASITE_STAGING_MYSQL-w $f -Ne "UPDATE meta SET meta_value='WBPS${PARASITE_VERSION}' WHERE meta_key='genebuild.version';"; fi;

  $PARASITE_STAGING_MYSQL-w $f -Ne "UPDATE meta SET meta_value='WBPS${PARASITE_VERSION}' WHERE meta_key='genebuild.version' AND meta_value='1';";

  $PARASITE_STAGING_MYSQL-w $f -Ne "UPDATE meta SET meta_value='repeatmask_customlib' WHERE meta_key='repeat.analysis' AND meta_value='repeatmask';";
  $PARASITE_STAGING_MYSQL-w $f -Ne "UPDATE meta SET meta_value='trf' WHERE meta_key='repeat.analysis' AND meta_value='repeatmask_repbase';";

  $PARASITE_STAGING_MYSQL-w $f -Ne "DELETE FROM meta WHERE meta_key NOT IN ('division', 'patch', 'schema_type', 'schema_version') AND species_id IS NULL;";

  #  $PARASITE_STAGING_MYSQL mysqldump $f mapping_session > $HANDOVER_SCRATCH/$f.mapping_session.sql;
  #  sed -i s,"\`old_assembly\` varchar(20) NOT NULL DEFAULT","\`old_assembly\` varchar(80) NOT NULL DEFAULT",g $HANDOVER_SCRATCH/$f.mapping_session.sql;
  #  sed -i s,"\`new_assembly\` varchar(20) NOT NULL DEFAULT","\`new_assembly\` varchar(80) NOT NULL DEFAULT",g $HANDOVER_SCRATCH/$f.mapping_session.sql;
  #  $PARASITE_STAGING_MYSQL-w $f -Ne "DROP TABLE mapping_session;";
  #  $PARASITE_STAGING_MYSQL-w $f < $HANDOVER_SCRATCH/$f.mapping_session.sql;
  #
  #  $PARASITE_STAGING_MYSQL mysqldump $f stable_id_event > $HANDOVER_SCRATCH/$f.stable_id_event.sql;
  #  sed -i s,"\`type\` enum('gene'\,'transcript'\,'translation') NOT NULL","\`type\` enum('gene'\,'transcript'\,'translation'\,'rnaproduct') NOT NULL",g $HANDOVER_SCRATCH/$f.stable_id_event.sql;
  #  $PARASITE_STAGING_MYSQL-w $f -Ne "DROP TABLE stable_id_event;";
  #  $PARASITE_STAGING_MYSQL-w $f < $HANDOVER_SCRATCH/$f.stable_id_event.sql;

  $PARASITE_STAGING_MYSQL-w $f -Ne "UPDATE transcript t
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
