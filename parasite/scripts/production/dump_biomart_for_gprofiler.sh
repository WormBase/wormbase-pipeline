#!/usr/bin/bash
# Dump for our g:Profiler friends
# Adapted from EnsEMBL, related:
# - Bio::EnsEMBL::Production::Pipeline::PipeConfig::MySQLDumping_conf
# - $ENSEMBL_CVS_ROOT_DIR/ensembl-production/modules/Bio/EnsEMBL/Production/Utils/MySQLDumping.sh
set -x
set -euo pipefail

STAGING_DIR=${PARASITE_SCRATCH}/dumps/WBPS${PARASITE_VERSION}/biomart/WBPS${PARASITE_VERSION}
RELEASE_DIR=/nfs/ftp/public/databases/wormbase/collaboration/UTE/WBPS${PARASITE_VERSION}


if [ ! -d "$STAGING_DIR" ] ; then
  TMP=${STAGING_DIR}.tmp
  mkdir -pv $TMP
  chmod 777 $TMP

  SRV_URL=$($PARASITE_STAGING_MYSQL details url)
  DBNAME=parasite_mart_${PARASITE_VERSION}
  OUTPUT_DIR=/homes/digri/test

  #Outputing txt files
  standaloneJob.pl Bio::EnsEMBL::Production::Pipeline::FileDump::MySQL_TXT \
    -db_url ${SRV_URL}${DBNAME} \
    -dbname ${DBNAME} \
    -output_dir ${TMP} \
    -dump_dir ${TMP} \
    -overwrite 1

  #Outputting sql files
  for table in $( $PARASITE_STAGING_MYSQL ${DBNAME} -Ne "SHOW TABLES" ); do
    mysqldump $( $PARASITE_STAGING_MYSQL-w details mysql ) \
    --quick --single-transaction \
    parasite_mart_${PARASITE_VERSION} $table > ${TMP}/$table.sql;
  done
    mv $TMP $STAGING_DIR
  fi

find $STAGING_DIR -type f -name '*.txt' | while read -r f; do
  gzip -v $f
done

sudo -u wormbase mkdir -pv $RELEASE_DIR/biomart
sudo -u wormbase rsync -av --delete $STAGING_DIR/ $RELEASE_DIR/biomart/

perl -MProductionMysql -E '
my %core_db_to_biomart_name;
my %bps;
for my $core_db (ProductionMysql->staging->core_databases) {
  #if ($core_db =~ /$m/) {
  my ($spe, $cies, $bp) = split "_", $core_db;
  push @{$bps{"${spe}_${cies}"}}, $bp;
  $core_db_to_biomart_name{$core_db} = ProductionMysql::core_db_to_biomart_name($core_db);
  #} else {
  #die "$core_db <no match>\n"
  #}
}
for my $core_db (sort keys %core_db_to_biomart_name){
  my ($species, $versions) = split "_core_", $core_db;
  my ($spe, $cies, $bp) = split "_", $species;
  my $pretty_name = sprintf("%s %s", ucfirst $spe, $cies);
  if ($bp){
    $pretty_name =~ s/sp(\d+)$/sp. $1/;
    $pretty_name =~ s/kr3021$/sp. KR3021/;
    $pretty_name =~ s/t(\d+)$/sp. T$1/;
    $pretty_name .= sprintf(" (%s)", uc $bp) if @{$bps{"${spe}_${cies}"}} > 1;
  } else { 
    my ($v) = $versions =~ /^(\d+)/;
    $pretty_name .= sprintf(" (E! comparator)", $v) if $v and $v > 70;
    $pretty_name .= sprintf(" (EG comparator)", $v) if $v and $v < 70;
  }
  say join "\t", $core_db, $pretty_name, $core_db_to_biomart_name{$core_db};

}

' | sudo -u wormbase tee $RELEASE_DIR/core_db_to_species_name_to_biomart_name.tsv

ls -latorhd $RELEASE_DIR/*






#old mysqldump was like:
#mysqldump $( $PARASITE_STAGING_MYSQL-w details mysql ) \
#  --quick --single-transaction \
#  -T $TMP \
#  parasite_mart_${PARASITE_VERSION}
#mv -v $TMP $STAGING_DIR