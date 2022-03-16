# WormBase ParaSite performs some customisations to the e! xrefs pipeline 
# to add wormbase-specific external references and it the populates the description
# display_xref_id fields of the gene and trascript table in a slight different way. 
# This script performs a quick datacheck to make sure that these fields have been populated.
dtn() {
  msg=$1
  dt=$(date '+%d/%m/%Y %H:%M:%S:');
  printf "$dt ${msg}\n"
}

fail_message() {
  dtn "ERROR: $cdb failed the $DC with count=${dccount}.";
  exit 1
}

evaluate_count() {
  dccount=$1
  cthreshold=$2
  if (( $dccount < $cthreshold ));
    then return 1
  fi;
}

perform_dc() {
  cdb=$1
  cthreshold=$2
  result_count=$3
  evaluate_count $result_count $cthreshold || fail_message $sqls $cdb $DC
  dtn "Success: $DC: got $result_count"
}

${PARASITE_STAGING_MYSQL} -Ne "SHOW DATABASES" | grep _${PARASITE_VERSION}_${ENSEMBL_VERSION}_${WORMBASE_VERSION} \
| while read cdb;
  do
  dtn "Processing $cdb";
  ###########################################################################
  DC="1/2: Gene display_xref_ids DC"
  SQLS="SELECT COUNT(*) FROM gene WHERE display_xref_id IS NOT NULL"
  threshold=10000
  dccount=$($PARASITE_STAGING_MYSQL $cdb -Ne "$SQLS")
  perform_dc $cdb $threshold $dccount $DC
  ###########################################################################
  DC="2/2: Gene description DC"
  SQLS="SELECT COUNT(*) FROM gene WHERE description IS NOT NULL"
  threshold=2000
  dccount=$($PARASITE_STAGING_MYSQL $cdb -Ne "$SQLS")
  perform_dc $cdb $threshold $dccount $DC
   ###########################################################################
  DC="Extra 1/2: Transcript display_xref_ids DC"
  SQLS="SELECT COUNT(*) FROM transcript WHERE display_xref_id IS NOT NULL"
  threshold=10000
  dccount=$($PARASITE_STAGING_MYSQL $cdb -Ne "$SQLS")
  dtn "Info: ${cdb} transcript display_xref_id count: ${dccount}"
  ###########################################################################
  DC="Extra 2/2: Transcript description Advisory DC"
  SQLS="SELECT COUNT(*) FROM transcript WHERE description IS NOT NULL"
  dccount=$($PARASITE_STAGING_MYSQL $cdb -Ne "$SQLS")
  dtn "Info: ${cdb} transcript description count: ${dccount}"
  printf "\n"
done
 