#!/usr/bin/env bash

# https://www.ebi.ac.uk/seqdb/confluence/display/WORMBASE/Creating+Core+Database#CreatingCoreDatabase-BUSCOandCEGMA
# 
# this is the code from the wiki page; avoids copying & pasting into a terminal,
# also allows species to be passed as arguments rather that iterating through
# all species in $PARASITE_DATA

submit_busco() {
   core_db=$1
   species=$(echo $core_db | cut -d'_' -f1,2,3)
   echo "Submitting BUSCO job for $core_db"
   MEM_MB=10000
   log=$PARASITE_SCRATCH/busco-annotation/WBPS${PARASITE_VERSION}/log
   mkdir -pv $log
   bsub \
      -o $log/$species.%J.out \
      -e $log/$species.%J.err \
      -n 8 -M${MEM_MB} -R "select[mem > ${MEM_MB}] rusage[mem=${MEM_MB}]" \
      bash -x $WORM_CODE/parasite/scripts/production/core-creation/busco-annotation.sh $core_db
}
 
get_db_for_species() {
    species=$1
    $PARASITE_STAGING_MYSQL -Ne "show databases like \"$species%\""
}

meta_table_has_key_pattern() {
    pattern=$1
    core_db=$2
    [ "$core_db" ] && [ 0 -lt $($PARASITE_STAGING_MYSQL $core_db -Ne "select count(*) from meta where meta_key like \"$pattern\" " ) ]
}

if [[ $1 =~ ^\-\-?h ]]; then
   echo "Usage: $0 [species1 [species2..speciesN] ]"
   exit
fi

# check required environment variables are set
if [ -z ${PARASITE_VERSION+x} ]; then
   echo "Error: the PARASITE_VERSION environment variable isn't set."
   exit 255
fi
if [ -z ${PARASITE_DATA+x} ]; then
   echo "Error: the PARASITE_DATA environment variable isn't set."
   exit 255
fi
if [ -z ${PARASITE_SCRATCH+x} ]; then
   echo "Error: the PARASITE_SCRATCH environment variable isn't set."
   exit 255
fi
if [ -z ${WORM_CODE+x} ]; then
   echo "Error: the PARASITE_SCRATCH environment variable isn't set."
   exit 255
fi

CORE_DBS="$@"

if [[ -z ${CORE_DBS} ]]; then
   # SPECIES=$( ls $PARASITE_DATA )
   CORE_DBS=$($PARASITE_STAGING_MYSQL -Ne "SHOW DATABASES LIKE \"%_core_${PARASITE_VERSION}_${ENSEMBL_VERSION}%\"")
fi

for this_core_db in ${CORE_DBS} ; do
   #meta_table_has_key_pattern %busco% this_core_db || submit_busco this_core_db
   submit_busco $this_core_db
done
