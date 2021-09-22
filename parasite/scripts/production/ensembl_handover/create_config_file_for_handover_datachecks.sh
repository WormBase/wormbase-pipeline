#!/usr/bin/bash

for VAR in  HANDOVER_STAGING_MYSQL    \
            PARASITE_STAGING_MYSQL    \
            PARASITE_SCRATCH          \
            ENSEMBL_VERSION           \
            EG_VERSION                \
            PARASITE_HANDOVER_VERSION \
            HANDOVER_CONF
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
    \"old_server_uri\" : [\"$(${PARASITE_STAGING_MYSQL} details url)\"]
  },
  \"history_file\" : null,
  \"output_dir\"   : \"${PARASITE_SCRATCH}/log/WBPS_HANDOVER_${PARASITE_HANDOVER_VERSION}/datachecks/\",
  \"output_file\"  : null
}" > ${HANDOVER_CONF}/datachecks_config.json


