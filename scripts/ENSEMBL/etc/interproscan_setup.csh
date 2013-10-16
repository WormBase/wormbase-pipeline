

setenv IPR_PIPELINE_NAME $1
setenv HIVEDB $IPR_PIPELINE_NAME
setenv PATH ${PATH}:$WORM_PACKAGES/ensembl/ensembl-hive/scripts
setenv EG_ENSEMBL_ROOT $WORM_PACKAGES/ensembl_genomes

setenv PIPELINE_DIR $PIPELINE/interpro_scratch


foreach lib_path  (${EG_ENSEMBL_ROOT}/eg-*/lib)
    setenv PERL5LIB ${lib_path}:${PERL5LIB}
end

setenv HIVE_URL  mysql://wormadmin:worms@${WORM_DBHOST}:${WORM_DBPORT}/${HIVEDB}


#Delete the old hive database

mysql --host=${WORM_DBHOST} --port=${WORM_DBPORT} --user=wormadmin --password=worms -e "DROP DATABASE IF EXISTS ${HIVEDB}"

# setup the new pipelines

init_pipeline.pl Bio::EnsEMBL::Hive::PipeConfig::InterProScan_conf -registry ${ENSEMBL_REGISTRY} -species $SPECIES -hive_host ${WORM_DBHOST} -hive_port ${WORM_DBPORT} -hive_user wormadmin -hive_password worms -hive_dbname ${HIVEDB} -pipeline_dir ${PIPELINE_DIR} -ensembl_cvs_root_dir ${EG_ENSEMBL_ROOT}



