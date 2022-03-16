# This script repairs the xref and the object_xref tables of all our core dbs
# as it has been observed that after xrefs some tables are getting corrupted
# causing GPAD load to run really slow. See ticket:
# https://www.ebi.ac.uk/panda/jira/browse/ENSINT-1136
dtn() {
  msg=$1
  dt=$(date '+%d/%m/%Y %H:%M:%S:');
  printf "$dt ${msg}\n"
}

${PARASITE_STAGING_MYSQL} -Ne "SHOW DATABASES" | grep _${PARASITE_VERSION}_${ENSEMBL_VERSION} \
| while read cdb;
  do
  dtn "Processing $cdb";
  dtn "Reparing $cdb object xref table...";
  ${PARASITE_STAGING_MYSQL}-w $cdb -Ne "REPAIR TABLE object_xref";
  dtn "Reparing $cdb xref table...";
  ${PARASITE_STAGING_MYSQL}-w $cdb -Ne "REPAIR TABLE xref";
  dtn "Done";
  printf "\n"
done

printf "\n"
dtn "Finished repairing. Happy GPAD Load :)";
###########################################################################