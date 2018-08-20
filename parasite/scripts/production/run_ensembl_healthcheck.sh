#!/usr/bin/bash

# Options we like the most:
# -d "core db"
# --include_groups , --include_tests, --exclude_tests

cd $ENSEMBL_CVS_ROOT_DIR/ensj-healthcheck

PERL5LIB=`pwd`/perl:$PERL5LIB
#TODO
# --include_tests 'BUSCOResults CEGMAResults'

if [ $# -eq 0 ] ; then
  ARGS=( --include_groups EGCore --exclude_tests 'DisplayXref DisplayXrefId DuplicateObjectXref EGCompareCoreSchema EnaSeqRegionName ENASeqRegionSynonyms GeneSource GoTermCount MultiDbCompareNames PermittedEgMeta ProductionMasterTables ProteinTranslation SampleSetting SchemaPatchesApplied TranscriptSource UniParc_Coverage UniProtKB_Coverage UniProtKB_Coverage org.ensembl.healthcheck.testcase.eg_core.SeqRegionCoordSystem ' )
else
  ARGS=( "$@" )
fi
#PRINT_ANT=1 if you need
export NO_JAR=1
./run-configurable-testrunner.sh $(mysql-pan-prod details script) $($PARASITE_STAGING_MYSQL details suffix_1) $($PREVIOUS_PARASITE_STAGING_MYSQL details suffix_2) \
  --production.database ensembl_production_parasite \
  --output Info \
  "${ARGS[@]}"
