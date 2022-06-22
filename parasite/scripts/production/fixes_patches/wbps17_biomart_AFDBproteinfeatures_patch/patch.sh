workdir=${PARASITE_SCRATCH}/biomart/WBPS${PARASITE_VERSION}/afdb_patch_v1

mkdir -p $workdir
echo "Workdir: ${workdir}"
echo ""

for species_path in $(ls -1 -d ${PARASITE_SCRATCH}/alphafold/WBPS${PARASITE_VERSION}/*_*_*);
    do species=$(basename ${species_path});
    if  perl -MProductionMysql -E 'say for ProductionMysql->staging->core_databases(@ARGV)' ""  | grep -q $species;
      then echo "Fixing:" $species;
      cdb=$(perl -MProductionMysql -E 'say for ProductionMysql->staging->core_databases(@ARGV)' "$species" | head -n1);
      sqlfile=$(ls -1 ~/martbuilder_wbps_${PARASITE_VERSION}/${cdb}/gene/*_gene__protein_feature_alphafold_import__dm.sql);
      sqlfile_basename=$(basename $sqlfile)
      tabname="${sqlfile_basename%.*}"
      maintabname=$(echo $tabname | sed -e s,"protein_feature_alphafold_import__dm","translation__main",g)
      if [ -f "$sqlfile" ]; then
        echo "    Found ${sqlfile} biomart sql.";
      else
        echo "File ${sqlfile} does not exist. Exiting...";
        exit 1;
      fi;
      speciesdir=${workdir}/${species}
      mkdir -p ${speciesdir}
      sqlfile_fixed=${speciesdir}/${sqlfile_basename}
      echo "    Copying to: ${sqlfile_fixed}"
      cat ${sqlfile} | \
      sed -r '/^alter table parasite_mart_\w+\.\w+_gene__translation__main add column \(protein_feature_alphafold_import_bool integer default 0\);$/d' | \
      sed -r '/^create index \w+ on parasite_mart_\w+\.\w+_gene__translation__main\(protein_feature_alphafold_import_bool\);$/d' > ${sqlfile_fixed}
      echo "    Datacheck: Make sure that only two lines have been removed from ${sqlfile}."
      count_pre_fixed=$(cat ${sqlfile} | wc -l)
      count_fixed=$(cat ${sqlfile_fixed} | wc -l)
      count_fixed_plus_two=$(expr $count_fixed + 2)
      if [ "${count_pre_fixed}" == "${count_fixed_plus_two}" ]; then
        echo "        Passed.";
      else
        echo "        FAILED. Exiting";
        exit 1;
      fi;
      if ${PARASITE_STAGING_MYSQL} parasite_mart_${PARASITE_VERSION} -Ne "SHOW TABLES;" | grep -q $tabname;
        then echo "    Found existing AFDB biomart table: ${tabname}";
      else
        echo "    Didn't find an existing AFDB biomart table. Not sure why. Exiting..."
        exit 1;
      fi;
      echo "    Removing existing AFDB biomart table: ${tabname}"
      $PARASITE_STAGING_MYSQL-w parasite_mart_${PARASITE_VERSION} -Ne "DROP TABLE ${tabname};"
      echo "    Creating the new correct AFDB biomart table: ${tabname}"
      $PARASITE_STAGING_MYSQL-w parasite_mart_${PARASITE_VERSION} < ${sqlfile_fixed}
      echo "    Datacheck: Make sure that the number of total NON-NULL AFDB entries in ${tabname} is the same with ${maintabname}:"
      count_afdb_table=$($PARASITE_STAGING_MYSQL parasite_mart_${PARASITE_VERSION} -Ne "SELECT COUNT(DISTINCT(translation_id_1068_key)) FROM ${tabname} WHERE hit_name_1048 IS NOT NULL;")
      count_main_table=$($PARASITE_STAGING_MYSQL parasite_mart_${PARASITE_VERSION} -Ne "SELECT COUNT(*) FROM ${maintabname} WHERE protein_feature_alphafold_import_bool IS NOT NULL;")
      if [ "${count_afdb_table}" == "${count_main_table}" ]; then
        echo "    Passed.";
      else
        echo ${count_afdb_table};
        echo ${count_main_table};
        echo "    FAILED. Exiting...";
        exit 1;
      fi;
      echo "Done"
      echo ""
    fi;
done
