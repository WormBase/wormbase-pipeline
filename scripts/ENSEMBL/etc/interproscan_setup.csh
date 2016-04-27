

setenv IPR_PIPELINE_NAME $1
shift

setenv HIVEDB $IPR_PIPELINE_NAME
setenv PATH ${PATH}:$WORM_PACKAGES/ensembl/ensembl-hive/scripts
setenv EG_ENSEMBL_ROOT $WORM_PACKAGES/ensembl_genomes

setenv PIPELINE_DIR $PIPELINE/interpro_scratch

setenv PERL5LIB ${WORM_SW_ROOT}/packages/ensembl/branches/master/ensembl-production/modules:${PERL5LIB}

setenv HIVE_URL  mysql://wormadmin:worms@${WORM_DBHOST}:${WORM_DBPORT}/${HIVEDB}


#Delete the old hive database

mysql --host=${WORM_DBHOST} --port=${WORM_DBPORT} --user=wormadmin --password=worms -e "DROP DATABASE IF EXISTS ${HIVEDB}"

# setup the new pipelines

init_pipeline.pl Bio::EnsEMBL::Production::Pipeline::PipeConfig::InterProScan_conf -registry ${ENSEMBL_REGISTRY} -hive_host ${WORM_DBHOST} -hive_port ${WORM_DBPORT} -hive_user wormadmin -hive_password worms -hive_dbname ${HIVEDB} -pipeline_dir ${PIPELINE_DIR} -production_lookup 0 $*


