#!/usr/bin/bash
# Run Ensembl's ID mapping pipeline, and our custom upload.
# 
# The .ini file is for the default case of previous_staging->staging
# If ran in the terminal, the script will let you examine it and change it.
# If you change it to e.g. point to another server,
# upload_mapping_as_history.pl will die. It's ok to run it manually after it fails.
#
# Remove cache folder to restart after a failed cache dump.
set -euo pipefail

species=$1
sourcedbname=$2

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null && pwd )"

base_dir=$PARASITE_SCRATCH/id_mapping/WBPS${PARASITE_VERSION}/$species
mkdir -pv $base_dir
if [ ! -s $base_dir/$species.ini ] ; then
  echo "basedir=$base_dir" > $base_dir/$species.ini
  if [ -z "$sourcedbname" ]; then
    perl $DIR/mapping_conf.pl --species $species >> $base_dir/$species.ini
  else
    perl $DIR/mapping_conf.pl --species $species --sourcedbname $sourcedbname >> $base_dir/$species.ini
  fi
fi
[ -t 0 ] && ${EDITOR:-vi} $base_dir/$species.ini

pushd $ENSEMBL_CVS_ROOT_DIR/ensembl/misc-scripts/id_mapping
if [ ! -d $base_dir/cache ]; then
  perl dump_cache.pl --conf $base_dir/$species.ini
  rm -rf $base_dir/matrix $base_dir/mapping
fi
perl id_mapping.pl --conf $base_dir/$species.ini
popd
if [ ! -s $base_dir/debug/gene_mappings.txt ] || [ ! -s $base_dir/tables/stable_id_event_new.txt ] ; then 
   >&2 echo "No results in $base_dir - exiting"
   exit 1
fi

core_db=$(perl -MProductionMysql -ne 'print ProductionMysql->staging->core_db($_);'<<< $species )
join -t $'\t' -11 -21 \
  <( join -t $'\t' -11 -22 \
     <( $PARASITE_STAGING_MYSQL "$core_db" -Ne 'select gene_id, stable_id from gene' | sort -k1,1 ) \
     <( sort -k 2,2 $base_dir/debug/gene_mappings.txt ) \
     | perl -lanE 'say join "\t", @F[3,1]' | sort -k1,1  \
   ) \
  <(perl -lanE 'say join "\t", @F[0,-1] if grep {$_ eq "gene"} @F ' $base_dir/tables/stable_id_event_new.txt | sort -k1,1 ) \
  > $base_dir/previous_to_current_id.tsv

if [ ! -s $base_dir/previous_to_current_id.tsv ] ; then
   >&2 echo "No mapping in $base_dir/previous_to_current_id.tsv - exiting "
   exit 1
fi
echo "Obtained mapping: $base_dir/previous_to_current_id.tsv"
$DIR/upload_mapping_as_history.pl --species "$species" --mapping $base_dir/previous_to_current_id.tsv 
