

setenv IPR_PIPELINE_NAME $1
shift

setenv HIVEDB "wormbase_${IPR_PIPELINE_NAME}"
setenv PATH ${PATH}:$WORM_PACKAGES/ensembl/ensembl-hive/scripts
setenv EG_ENSEMBL_ROOT $WORM_PACKAGES/ensembl_genomes

setenv PIPELINE_DIR $PIPELINE/interpro_scratch

setenv PERL5LIB ${WORM_SW_ROOT}/packages/ensembl/branches/master/ensembl-production/modules:${PERL5LIB}

setenv HIVE_URL  mysql://wormadmin:worms@${WORM_DBHOST}:${WORM_DBPORT}/${HIVEDB}

setenv INTERPROSCAN_VERSION `interproscan.sh --version|perl -ne 'print $1 if /InterProScan version (\S+)/'`

#Delete the old hive database

mysql --host=${WORM_DBHOST} --port=${WORM_DBPORT} --user=wormadmin --password=worms -e "DROP DATABASE IF EXISTS $HIVEDB"

sleep 1
# setup the new pipelines

init_pipeline.pl Bio::EnsEMBL::Production::Pipeline::PipeConfig::ProteinFeatures_conf -registry ${ENSEMBL_REGISTRY} --host ${WORM_DBHOST} --port ${WORM_DBPORT} --user wormadmin --password worms -pipeline_dir ${PIPELINE_DIR} -pipeline_name ${IPR_PIPELINE_NAME} -production_lookup 0 -delete_existing=1 -interproscan_version $INTERPROSCAN_VERSION $*


