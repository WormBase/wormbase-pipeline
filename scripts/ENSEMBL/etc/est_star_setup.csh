
# small csh script to help set up the EST STAR pipeline

# the MAX_INTRON_LENGTH environnment variable is expected to already be set
# this is the largest known intron in the genome

# the LOGIC_NAME environnment variable is expected to already be set
# this is the logic_name to be run (one of embl_star, washu_star, nembase_star)

# the SPECIES environnment variable is expected to already be set
# this is the normal WormBase short production name

# two or more different species can be run at the same time.
# you CAN NOT run two or more analyses in the same species at the same time.

if (! $?SPECIES) then
    echo "ERROR: SPECIES environment variable not set."
    exit 1
endif

if (! $?MAX_INTRON_LENGTH) then
    echo "ERROR: MAX_INTRON_LENGTH environment variable not set."
    exit 1
endif

if (! $?LOGIC_NAME) then
    echo "ERROR: You should set the LOGIC_NAME environment variable. (One of 'embl_star', 'washu_star', 'nembase_star')"
    exit 1
endif

if ($LOGIC_NAME == "embl_star") then
  setenv EST_FILE $wormbase/BUILD_DATA/cDNA/$SPECIES/EST
else if ($LOGIC_NAME == "washu_star") then
  setenv EST_FILE $wormbase/BUILD_DATA/cDNA/$SPECIES/Nematode.net
else if ($LOGIC_NAME == "nembase_star") then
  setenv EST_FILE $wormbase/BUILD_DATA/cDNA/$SPECIES/Nembase
else
    echo "ERROR: Unknown LOGIC_NAME: $LOGIC_NAME"
    exit 1
endif


setenv PIPELINE_CODE_DIR $WORM_PACKAGES/ensembl_genomes
setenv PIPELINE_DATA_DIR /nfs/nobackup/ensemblgenomes/wormbase/BUILD/EST-data/$SPECIES
set path = (/nfs/panda/ensemblgenomes/external/STAR $path)
setenv CORE_DB_NAME worm_ensembl_$SPECIES
setenv DB_HOST  mysql-wormbase-pipelines
setenv DB_PORT 4331
setenv DB_USER wormadmin
setenv DB_PASS worms
setenv HIVE_DB_HOST mysql-wormbase-pipelines
setenv HIVE_DB_PORT 4331
setenv HIVE_DB_USER wormadmin
setenv HIVE_DB_PASS worms
setenv FORCE 0
setenv PERL5LIB ${WORM_PACKAGES}/ensembl_genomes/eg-estalignfeature/lib/:${PERL5LIB}
setenv PERL5LIB ${WORM_PACKAGES}/ensembl_genomes/eg-pipelines/modules/:${PERL5LIB}
if (! -d $PIPELINE_DATA_DIR) then
    mkdir $PIPELINE_DATA_DIR
else
    # tidy up any existing files
    rm -rf $PIPELINE_DATA_DIR/*
endif

cd $PIPELINE_DATA_DIR

if (! -e $EST_FILE) then
    echo "ERROR: The EST sequence file ${EST_FILE} has not been set up yet."
    exit 1
endif

#if (-z $EST_FILE) then
#    echo "The EST sequence file ${EST_FILE} has no sequences in it. It is not worth running the pipeline on this."
#    exit 0
#endif

setenv HIVEDB ${USER}_EST_Pipeline_${SPECIES}
setenv HIVE_URL mysql://${HIVE_DB_USER}:${HIVE_DB_PASS}@${HIVE_DB_HOST}:${HIVE_DB_PORT}/${HIVEDB}

# Delete the old hive database and any old results
mysql --host=${DB_HOST} --port=${DB_PORT} --user=wormadmin --password=worms -e "DROP DATABASE IF EXISTS ${HIVEDB}"
mysql --host=${DB_HOST} --port=${DB_PORT} --user=wormadmin --password=worms -e "USE ${CORE_DB_NAME}; DELETE FROM dna_align_feature WHERE analysis_id IN (SELECT analysis_id FROM analysis WHERE logic_name = '${LOGIC_NAME}')"


perl ${PIPELINE_CODE_DIR}/ensembl-hive/scripts/init_pipeline.pl Bio::EnsEMBL::Hive::PipeConfig::EST_STAR_pipeline_conf \
   -registry ${ENSEMBL_REGISTRY} -species ${SPECIES} -ensembl_cvs_root_dir ${PIPELINE_CODE_DIR} \
   -hive_password ${HIVE_DB_PASS} -hive_port ${HIVE_DB_PORT} -hive_host ${HIVE_DB_HOST} -hive_user ${HIVE_DB_USER} -hive_dbname ${HIVEDB} \
   -pipeline_dir ${PIPELINE_DATA_DIR} -core_db_name ${CORE_DB_NAME} -est_file ${EST_FILE} -force 0 -logic_name ${LOGIC_NAME} -max_intron_length ${MAX_INTRON_LENGTH}


#${PIPELINE_CODE_DIR}/ensembl-hive/scripts/beekeeper.pl -reg_conf ${ENSEMBL_REGISTRY} -url ${HIVE_URL} -sync

#${PIPELINE_CODE_DIR}/ensembl-hive/scripts/beekeeper.pl -reg_conf ${ENSEMBL_REGISTRY} -url ${HIVE_URL} -submit_workers_max 100 -sleep 2 -life_span 2000 -loop

