#!/usr/bin/env bash

dir=$PARASITE_SCRATCH/xrefs/WBPS${PARASITE_VERSION}/xref_update
mkdir -pv $dir

# Takes the E! xref sources template and adds wormbase-specific sources
$WORM_CODE/parasite/scripts/production/xref/our_ini.pl "$@" > $dir/xref_config.ini

# Converts the above to JSON
$WORM_CODE/parasite/scripts/production/xref/ini_to_json.pl $dir/xref_config.ini > $dir/xref_sources.json

# Sort out PERL5LIBs
# The pipeline looks for XrefMapper::$species on the module path
# see: Bio::EnsEMBL::Production::Pipeline::Xrefs::Base::get_xref_mapper
# We have overriden the xref modules- make one flavour of module for core wormbase species and another flavour for non-wormbase parasite species.

rm -rf $dir/lib/XrefMapper
mkdir -pv $dir/lib/XrefMapper
grep -q "ensembl/misc-scripts/xref_mapping" <<< "$PERL5LIB" || export PERL5LIB=$ENSEMBL_CVS_ROOT_DIR/ensembl/misc-scripts/xref_mapping:$PERL5LIB

PERL5LIB=$WORM_CODE/parasite/modules:$PERL5LIB

$WORM_CODE/parasite/scripts/production/xref/ini_to_species.pl $dir/xref_config.ini | while read species taxon; do
    cp $(perldoc -lm XrefMapper::$taxon.pm ) $dir/lib/XrefMapper/$species.pm
    sed -i "s/package XrefMapper::$taxon/package XrefMapper::$species/" $dir/lib/XrefMapper/$species.pm
done
  
grep -q "$dir/lib" <<< "$PERL5LIB" || export PERL5LIB=$dir/lib:$PERL5LIB

mkdir -pv $dir/sql
cp -v $ENSEMBL_CVS_ROOT_DIR/ensembl/misc-scripts/xref_mapping/sql/table.sql $dir/sql/table.sql
$ENSEMBL_CVS_ROOT_DIR/ensembl/misc-scripts/xref_mapping/xref_config2sql.pl $dir/xref_config.ini  > $dir/sql/populate_metadata.sql
   
# $dir/sql/populate_metadata.sql has to be newer than $dir/xref_config.ini, otherwise the script will try to repopulate it.
if ! [ -s "$dir/sql/populate_metadata.sql" ] ; then echo "Error, not populated: $dir/sql/populate_metadata.sql " ; return 1 ; fi

TMP=$dir/tmp perl -MProductionMysql -E 'map {say "$ENV{TMP}/$_"} (ProductionMysql->staging->species(@ARGV))' "$@" | xargs rm -rfv

#Create the source database we need
mysql-ps-prod-1-w -Ne "DROP DATABASE IF EXISTS WBPS${PARASITE_VERSION}_xref_source";
mysql-ps-prod-1-w -Ne "CREATE DATABASE WBPS${PARASITE_VERSION}_xref_source";
mysql-ens-meta-prod-1 mysqldump -d xref_source | mysql-ps-prod-1-w -D WBPS${PARASITE_VERSION}_xref_source  

XREF_SETUP_DIR=$dir
export XREF_SETUP_DIR
