#!/usr/bin/env bash

# https://www.ebi.ac.uk/seqdb/confluence/display/WORMBASE/Creating+Core+Database#CreatingCoreDatabase-BUSCOandCEGMA
# 
# this is the code from the wiki page; avoids copying & pasting into a terminal,
# also allows species to be passed as arguments rather that iterating through
# all species in $PARASITE_DATA

submit_busco() {
   species=$1
   echo "Submitting BUSCO job for $species"
   MEM_MB=10000
   log=$PARASITE_SCRATCH/busco/WBPS${PARASITE_VERSION}/log
   mkdir -pv $log
   bsub \
      -o $log/$species.%J.out \
      -e $log/$species.%J.err \
      -n 8 -M${MEM_MB} -R "select[mem > ${MEM_MB}] rusage[mem=${MEM_MB}]" \
      bash -x $WORM_CODE/parasite/scripts/production/core-creation/busco.sh $PARASITE_DATA/$species/${species}.fa
}

submit_busco_annotation() {
   core_db=$1
   species=$(echo $core_db | cut -d'_' -f1,2,3)
   echo "Submitting BUSCO annotation job for $core_db"
   MEM_MB=10000
   log=$PARASITE_SCRATCH/busco-annotation/WBPS${PARASITE_VERSION}/log
   mkdir -pv $log
   bsub \
      -o $log/$species.%J.out \
      -e $log/$species.%J.err \
      -n 8 -M${MEM_MB} -R "select[mem > ${MEM_MB}] rusage[mem=${MEM_MB}]" \
      bash -x $WORM_CODE/parasite/scripts/production/core-creation/busco-annotation.sh $core_db
}

submit_cegma() {
   species=$1
   echo "Submitting CEGMA job for $species"
   MEM_MB=1500
   log=$PARASITE_SCRATCH/cegma/WBPS${PARASITE_VERSION}/log
   mkdir -pv $log
   bsub \
     -o $log/$species.%J.out \
     -e $log/$species.%J.err \
     -n 8 -M${MEM_MB} -R "select[mem > ${MEM_MB}] rusage[mem=${MEM_MB}]" \
     bash -x $WORM_CODE/parasite/scripts/production/core-creation/cegma.sh $PARASITE_DATA/$species/${species}.fa     
}
 
get_db_for_species() {
    species=$1
    $PARASITE_STAGING_MYSQL -Ne "show databases like \"$species%\""
}
 
meta_table_has_key_pattern() {
    pattern=$1
    species=$2
    core_db=$(get_db_for_species $species)
    echo $core_db
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

SPECIES="$@"

if [[ -z ${SPECIES} ]]; then
   # SPECIES=$( ls $PARASITE_DATA )
   SPECIES=$( find $PARASITE_DATA  -maxdepth 1 -mindepth 1 -type d -exec basename {} \; )
fi

for this_species in ${SPECIES} ; do
   core_db=$(get_db_for_species $species)
   meta_table_has_key_pattern %assembly.busco% $this_species || submit_busco $this_species
   meta_table_has_key_pattern %cegma% $this_species || submit_cegma $this_species
   meta_table_has_key_pattern %annotation.busco% $this_species || submit_busco_annotation $core_db
done
