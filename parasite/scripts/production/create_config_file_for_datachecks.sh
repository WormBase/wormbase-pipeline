#!/usr/bin/bash

for VAR in  PREVIOUS_PARASITE_STAGING_MYSQL    \
            PARASITE_CONF          \
            PARASITE_SCRATCH
do
   if [ -z $(eval "echo \$$VAR") ]; then
      echo "$VAR must be set in your environment before running $0"
      exit 254
   fi
done

# bale out immediately on error
set -e
printf "{
  \"datacheck_params\" : {
    \"registry_file\"  : null,
    \"server_uri\"     : null,
    \"old_server_uri\" : [\"$(${PREVIOUS_PARASITE_STAGING_MYSQL} details url)\"]
  },
  \"history_file\" : null,
  \"output_dir\"   : \"${PARASITE_SCRATCH}/datachecks/\",
  \"output_file\"  : null
}" > ${PARASITE_CONF}/datachecks_config.json