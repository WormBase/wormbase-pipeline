#!/usr/bin/env bash

# https://www.ebi.ac.uk/seqdb/confluence/display/WORMBASE/Creating+Core+Database#CreatingCoreDatabase-BUSCOandCEGMA
# 
# this is the code from the wiki page; avoids copying & pasting into a terminal,
# also allows species to be passed as arguments rather that iterating through
# all species in $PARASITE_DATA

submit_busco() {
   core_db=$1
   species=$(echo $core_db | cut -d'_' -f1,2);
   full_bioproject=$(echo $core_db | cut -d'_' -f3);
   bioproject=$(${WORM_CODE}/parasite/scripts/production/get_bioproject.py $core_db);
   ftpbioproject=$(${WORM_CODE}/parasite/scripts/production/get_ftp_bioproject.py $core_db);
   fasta=${PARASITE_FTP}/releases/WBPS${PARASITE_VERSION}/species/${species}/${ftpbioproject^^}/${species}.${ftpbioproject^^}.WBPS${PARASITE_VERSION}.genomic.fa.gz
   if [[ -f "$fasta" ]];
    then
      :
    else
      fasta=$PARASITE_DATA/${species}_${full_bioproject}/${species}_${full_bioproject}.fa
      if [[ ! -f $fasta ]];
        then
          echo "Could not find a fasta file. Exiting.";
      fi;
    fi;
   #echo "Submitting BUSCO job for $species"
   MEM_MB=10000
   log=$PARASITE_SCRATCH/busco/WBPS${PARASITE_VERSION}/log
   #echo $fasta
   mkdir -pv $log
    bsub \
    -o $log/$species.%J.out \
    -e $log/$species.%J.err \
    -n 8 -M${MEM_MB} -R "select[mem > ${MEM_MB}] rusage[mem=${MEM_MB}]" \
        bash -x $WORM_CODE/parasite/scripts/production/core-creation/busco-assembly.sh ${species}_${bioproject} ${fasta} ${core_db}
}

submit_busco_update_assembly() {
   core_db=$1
   species=$(echo $core_db | cut -d'_' -f1,2);
   bioproject=$(${WORM_CODE}/parasite/scripts/production/get_bioproject.py $core_db);
   ftpbioproject=$(${WORM_CODE}/parasite/scripts/production/get_ftp_bioproject.py $core_db);
   fasta=${PARASITE_FTP}/releases/WBPS${PARASITE_VERSION}/species/${species}/${ftpbioproject^^}/${species}.${ftpbioproject^^}.WBPS${PARASITE_VERSION}.genomic.fa.gz
   #echo "Submitting BUSCO job for ${species} with FASTA:${fasta}"
   MEM_MB=10000
   log=$PARASITE_SCRATCH/busco/WBPS${PARASITE_VERSION}/log
   mkdir -pv $log
   bsub \
   -o $log/$species.%J.out \
   -e $log/$species.%J.err \
   -n 8 -M${MEM_MB} -R "select[mem > ${MEM_MB}] rusage[mem=${MEM_MB}]" \
   bash -x $WORM_CODE/parasite/scripts/production/core-creation/busco-assembly.sh ${species}_${bioproject} $fasta $core_db
}

submit_busco_annotation() {
   core_db=$1
   species=$(echo $core_db | cut -d'_' -f1,2,3)
   echo "Submitting BUSCO annotation job for $core_db"
   MEM_MB=10000
   log=$PARASITE_SCRATCH/busco-annotation/WBPS${PARASITE_VERSION}/popdb_log
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
    [ "$core_db" ] && [ 0 -lt $($PARASITE_STAGING_MYSQL $core_db -Ne "select count(*) from meta where meta_key like \"$pattern\";" ) ]
}

meta_table_has_key_pattern_core_db() {
    pattern=$1
    core_db=$2
    [ "$core_db" ] && [ 0 -lt $($PARASITE_STAGING_MYSQL $core_db -Ne "select count(*) from meta where meta_key like \"$pattern\";" ) ]
}

core_db_has_updated_annotation() {
    core_db=$1
    species=$(echo $core_db | cut -d'_' -f1,2);
    full_bioproject=$(echo $core_db | cut -d'_' -f3);
    bioproject=$(${WORM_CODE}/parasite/scripts/production/get_bioproject.py $core_db);
    ftpbioproject=$(${WORM_CODE}/parasite/scripts/production/get_ftp_bioproject.py $core_db);
    fasta=$PARASITE_DATA/${species}_${full_bioproject}/${species}_${full_bioproject}.fa
    [[ -f $fasta ]]
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

if [[ -z ${ALL_CORE_DBS} ]]; then
   # SPECIES=$( ls $PARASITE_DATA )
   ALL_CORE_DBS=$(for f in $($PARASITE_STAGING_MYSQL -Ne "SHOW DATABASES LIKE \"%_core_${PARASITE_VERSION}_%\";"); do echo $f; done)
fi

for this_core_db in ${ALL_CORE_DBS} ; do
  if core_db_has_updated_annotation $this_core_db;
    then
      submit_busco $this_core_db; echo "Running BUSCO Assembly for ${this_core_db}";
      submit_busco_annotation $this_core_db; echo "Running BUSCO Annotation for ${this_core_db}";
    else
      meta_table_has_key_pattern_core_db %annotation.busco3% $this_core_db || { submit_busco_annotation $core_db; echo "Running BUSCO Annotation for $this_core_db"; }
      meta_table_has_key_pattern_core_db %assembly.busco3% $this_core_db || { submit_busco $this_core_db; echo "Running BUSCO Assembly for $this_core_db"; }
  fi
done