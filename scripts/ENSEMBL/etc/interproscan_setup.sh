

export IPR_PIPELINE_NAME="$1"
shift

export EHIVE_PROD2_DBPASS="$1"
shift

export HIVEDB="wormbase_${IPR_PIPELINE_NAME}"
export PATH=${PATH}:$WORM_PACKAGES/ensembl/ensembl-hive/scripts
export EG_ENSEMBL_ROOT=$WORM_PACKAGES/ensembl_genomes

export PIPELINE_DIR=$PIPELINE/interpro_scratch

export PERL5LIB=${WORM_SW_ROOT}/packages/ensembl/branches/master/ensembl-production/modules:${PERL5LIB}

export HIVE_URL=mysql://${EHIVE_PROD2_DBUSER}:${EHIVE_PROD2_DBPASS}@${EHIVE_PROD2_DBHOST}:${EHIVE_PROD2_DBPORT}/${HIVEDB}

export INTERPROSCAN_VERSION=$(interproscan.sh --version|perl -ne 'print $1 if /InterProScan version (\S+)/')

#Delete the old hive database

mysql --host=${EHIVE_PROD2_DBHOST} --port=${EHIVE_PROD2_DBPORT} --user=${EHIVE_PROD2_DBUSER} --password=${EHIVE_PROD2_DBPASS} -e "DROP DATABASE IF EXISTS $HIVEDB"

sleep 1
# setup the new pipelines

init_pipeline.pl Bio::EnsEMBL::Production::Pipeline::PipeConfig::ProteinFeatures_conf -registry ${ENSEMBL_REGISTRY} --host ${EHIVE_PROD2_DBHOST} --port ${EHIVE_PROD2_DBPORT} --user ${EHIVE_PROD2_DBUSER} --password ${EHIVE_PROD2_DBPASS} -pipeline_dir ${PIPELINE_DIR} -pipeline_name ${IPR_PIPELINE_NAME} -production_lookup 0 -delete_existing=1 -interproscan_version $INTERPROSCAN_VERSION $*


