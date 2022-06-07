for species_path in $(ls -1 -d ${PARASITE_SCRATCH}/alphafold/WBPS${PARASITE_VERSION}/*_*_*);
    do species=$(basename ${species_path});
    afdb_pf_count=$($PARASITE_STAGING_MYSQL $cdb -Ne "SELECT COUNT(*) FROM protein_feature t1 JOIN analysis t2 ON t1.analysis_id=t2.analysis_id WHERE t2.logic_name='alphafold_import';")
    echo "Found ${afdb_pf_count} AlphaFold protein features in $cdb";
    $PARASITE_STAGING_MYSQL $cdb -Ne "SELECT hit_description FROM protein_feature t1 JOIN analysis t2 USING(analysis_id) WHERE t2.logic_name='alphafold_import' GROUP BY translation_id HAVING COUNT(translation_id)>1;" | while read mres;
        do printf "${species}\t${mres}\n" >> $duplicated_translations_file;
      done;
      printf "Done\n\n";
    fi;
  done