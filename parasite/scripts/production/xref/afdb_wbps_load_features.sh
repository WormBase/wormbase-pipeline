duplicated_translations_file=${PARASITE_SCRATCH}/alphafold/WBPS${PARASITE_VERSION}/afdb_wbps_load_features.duplicates.tsv
rm $duplicated_translations_file
touch $duplicated_translations_file

for species_path in $(ls -1 -d ${PARASITE_SCRATCH}/alphafold/WBPS${PARASITE_VERSION}/*_*_*);
    do species=$(basename ${species_path});
    if  perl -MProductionMysql -E 'say for ProductionMysql->staging->core_databases(@ARGV)' ""  | grep -q $species;
      then echo "Found mappings for:" $species;
      cdb=$(perl -MProductionMysql -E 'say for ProductionMysql->staging->core_databases(@ARGV)' "$species" | head -n1)
      echo "Core db: $cdb";
      echo "Working...";
      #Check if afdb protein features have already been added for this core db
      afdb_pf_count=$($PARASITE_STAGING_MYSQL $cdb -Ne "SELECT COUNT(*) FROM protein_feature t1 JOIN analysis t2 ON t1.analysis_id=t2.analysis_id WHERE t2.logic_name='alphafold_import';")
      if (( $afdb_pf_count < 0 ));
        then echo "Found ${afdb_pf_count} AlphaFold protein features. Skipping $cdb.";
        else echo "No existing AlphaFold protein features in ${cdb}. Progressing...";
          #standaloneJob.pl ${WORM_CODE}/parasite/modules/AFDB_import/HiveLoadAlphaFoldDBProteinFeatures.pm \
          #-alpha_path ${species_path}/alpha_mappings.txt \
          #$(${PARASITE_STAGING_MYSQL}-w details script | sed s,'--','-core_db',g) \
          #-core_dbname $cdb \
          #-cs_version $($PARASITE_STAGING_MYSQL $cdb -Ne "SELECT meta_value FROM meta WHERE meta_key='assembly.name';") \
          #-species $species \
          #-rest_server https://www.ebi.ac.uk/gifts/api/;
          afdb_pf_count=$($PARASITE_STAGING_MYSQL $cdb -Ne "SELECT COUNT(*) FROM protein_feature t1 JOIN analysis t2 ON t1.analysis_id=t2.analysis_id WHERE t2.logic_name='alphafold_import';");
          am_count=$(cat ${species_path}/alpha_mappings.txt | wc -l);
          echo "${afdb_pf_count} protein features have been added from ${am_count} alpha mappings";
          echo "Storing duplicate translation->Uniprot mappings to ${duplicated_translations_file} for ${species}..."
          $PARASITE_STAGING_MYSQL $cdb -Ne "SELECT hit_description FROM protein_feature t1 JOIN analysis t2 USING(analysis_id) WHERE t2.logic_name='alphafold_import' GROUP BY translation_id HAVING COUNT(translation_id)>1;" | while read mres;
            do printf "${species}\t${mres}\n" >> $duplicated_translations_file;
          done;
          printf "Done\n\n";
      fi;
    fi;
  done