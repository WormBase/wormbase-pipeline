for DB in $($PARASITE_STAGING_MYSQL -Ne "SHOW DATABASES" | grep core_${PARASITE_VERSION}_${ENSEMBL_VERSION}_${PREVIOUS_WORMBASE_VERSION});
  do NEW_DB=$(sed "s/${PREVIOUS_WORMBASE_VERSION}$/${WORMBASE_VERSION}/g" <<< $DB)
    if [[ $($PARASITE_STAGING_MYSQL -Ne "SHOW DATABASES" | grep $NEW_DB) ]];
      then
        echo "database ${NEW_DB} already exists. Skipping";
      else
        echo "Creating ${NEW_DB}...";
        if ! ${PARASITE_STAGING_MYSQL}-w -e "CREATE DATABASE ${NEW_DB}"; then
          echo "Could not create ${NEW_DB} - should this script be running? Bailing out."
          exit 1
        fi
        echo "Dumping ${DB} to ${NEW_DB}";
        ${PARASITE_STAGING_MYSQL} mysqldump $DB | ${PARASITE_STAGING_MYSQL}-w ${NEW_DB}
        echo "Droping ${DB}"
        ${PARASITE_STAGING_MYSQL}-w -e "DROP DATABASE ${DB}"
        echo "Done."
    fi;
  done