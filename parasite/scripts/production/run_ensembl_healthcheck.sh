
# Options we like the most:
# -d "core db"
# --include_groups , --include_tests, --exclude_tests

cd /nfs/production/panda/ensemblgenomes/wormbase/software/packages/ensembl/ensj-healthcheck


PERL5LIB=`pwd`/perl:$PERL5LIB
#TODO
# --include_tests 'BUSCOResults CEGMAResults'

[ $# -eq 0 ] && $@=" --include_groups EGCore --exclude_tests 'DisplayXref DisplayXrefId DuplicateObjectXref EGCompareCoreSchema EnaSeqRegionName ENASeqRegionSynonyms GeneSource GoTermCount MultiDbCompareNames PermittedEgMeta ProductionMasterTables ProteinTranslation SampleSetting SchemaPatchesApplied TranscriptSource UniParc_Coverage UniProtKB_Coverage UniProtKB_Coverage "


./run-configurable-testrunner.sh $(mysql-pan-prod details script) $($PARASITE_STAGING_MYSQL details suffix_1) $($PREVIOUS_PARASITE_STAGING_MYSQL details suffix_2) \
  --production.database ensembl_production_parasite \
  --output Info \
  "$@"
