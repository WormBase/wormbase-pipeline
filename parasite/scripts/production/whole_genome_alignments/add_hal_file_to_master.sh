# manual definition
export HAL_SUFFIX=$1
export HAL_FILE=${PARASITE_SCRATCH}/cactus_alignments/WBPS${PARASITE_VERSION}/multi/hal_files/WBPS${PARASITE_VERSION}_${HAL_SUFFIX}.hal

mv ${PARASITE_SCRATCH}/cactus_alignments/WBPS${PARASITE_VERSION}/WBPS${PARASITE_VERSION}_${HAL_SUFFIX}.hal ${HAL_FILE}

export HAL_COLLECTION_NAME=parasitehal-${HAL_SUFFIX}
export HAL_COMPARA_SERVER=$PARASITE_STAGING_MYSQL
export HAL_COMPARA_DB=ensembl_compara_parasite_${PARASITE_VERSION}_${ENSEMBL_VERSION}
export HAL_MASTER_SERVER=mysql-ps-prod-1
export HAL_MASTER_DB=ensembl_compara_master_parasite

export HAL_COMPARA_DB_URL=$($HAL_COMPARA_SERVER-w details url $HAL_COMPARA_DB)
export HAL_MASTER_DB_URL=$($HAL_MASTER_SERVER-w details url $HAL_MASTER_DB)

# automatic definition
export COMPARA_HAL_DIR=$(awk -F '/multi' '{print $1}' <<<$HAL_FILE)
export HAL_MLSS_URL="#base_dir#/multi$(awk -F '/multi' '{print $2}' <<<$HAL_FILE)"

# perl script extracting genomes from hal
PERL_SCRIPT=${WORM_CODE}/parasite/scripts/production/whole_genome_alignments/dump_hal_genome_ids.pl

# our version of the codebase
export ENSEMBL_ROOT_DIR=${ENSEMBL_CVS_ROOT_DIR}

# Check if the hal file exists
if [ -f "$HAL_FILE" ];
  then :;
  else echo "File $HAL_FILE does not exist. Exiting..."
  exit 1
fi

# hal genomes
HAL_GENOMES=$(perl $PERL_SCRIPT --hal_file $HAL_FILE)
HAL_GENOMES_COUNT=$(echo $HAL_GENOMES | tr " " "\n" | wc -l)

echo $HAL_GENOMES | tr " " "\n" > ${PARASITE_CONF}/compara.${HAL_SUFFIX}.species_list_for_hal


if [ "$HAL_GENOME_DB_IDS_COUNT" -ne "$HAL_GENOMES_COUNT" ];
    then
      echo "WARNING: Not every genome in your HAL file has a genome_db in the genome_db table of your master database.";
      echo "Updating master DB:";
      # add missing genomes to the genome_Db table and add assign them to collections
      perl ${WORM_CODE}/parasite/scripts/production/compara/parasite_update_compara_master_for_hal.pl \
        -reg_conf $PARASITE_CONF/compara.registry.pm \
        -collection $HAL_COLLECTION_NAME \
        -masterdbname $HAL_MASTER_DB \
        -taxdbname ncbi_taxonomy_parasite \
        -comparacode $ENSEMBL_CVS_ROOT_DIR/ensembl-compara \
        -sfile ${PARASITE_CONF}/compara.${HAL_SUFFIX}.species_list_for_hal \
        -locators
fi


# Datachecks
if [ "$HAL_GENOME_DB_IDS_COUNT" -ne "$HAL_GENOMES_COUNT" ];
    then
      echo "WARNING: Not every genome in your HAL file has a genome_db in the genome_db table of your master database despite we ran parasite_update_compara_master_for_hal.pl. Exiting"; exit 1;
fi

HAL_GENOME_DB_IDS=$(\
for halg in $HAL_GENOMES; do
    $HAL_MASTER_SERVER $HAL_MASTER_DB -Ne "SELECT genome_db_id FROM genome_db WHERE name LIKE '$halg';";
done)
HAL_GENOME_DB_IDS_COUNT=$(echo $HAL_GENOME_DB_IDS | tr " " "\n" | wc -l)
HAL_GENOME_DB_ID=$(echo $HAL_GENOME_DB_IDS | tr " " ",")

# Naming mapping:
NAMING_MAPPING="{
    $(for i in $HAL_GENOMES; do
        name=$(cut -d '.' -f1 <<< "$i" ); \
        assembly=$(cut -d '.' -f2- <<< "$i" ); \
        id="$($HAL_MASTER_SERVER $HAL_MASTER_DB -e "SELECT genome_db_id AS '' FROM genome_db WHERE name LIKE '$name'")"; \
        [ "$id" != "" ] && printf "$id  => '$i'"; \
    done | sed '/^[[:space:]]*$/d' | tr '\n' ',') \
}"

# Check:
NAMING_MAPPING_COUNT=$(echo $NAMING_MAPPING | tr "," "\n" | wc -l)
if [ "$NAMING_MAPPING_COUNT" -ne "$HAL_GENOMES_COUNT" ];
    then echo "WARNING: Not every genome in your HAL file has a genome_db in the genome_db table of your master database."; exit 1;
fi


perl ${ENSEMBL_ROOT_DIR}/ensembl-compara/scripts/pipeline/create_mlss.pl \
--compara $HAL_MASTER_DB_URL --reg_conf $PARASITE_CONF/compara.registry.pm \
--source wormbase --species_set_name collection-$HAL_COLLECTION_NAME \
--method_link_type CACTUS_HAL --genome_db_id ${HAL_GENOME_DB_ID} --url "'$HAL_MLSS_URL'" --release

$HAL_MASTER_SERVER-w $HAL_MASTER_DB -Ne "UPDATE method_link_species_set SET url=\"$HAL_MLSS_URL\" WHERE url=\"'$HAL_MLSS_URL'\";"

MLSS_HAL_ID=$($HAL_MASTER_SERVER $HAL_MASTER_DB -e "SELECT mlss.method_link_species_set_id AS '' FROM method_link_species_set AS mlss JOIN species_set_header AS ssh USING(species_set_id) WHERE ssh.name = 'collection-$HAL_COLLECTION_NAME'" | tail -1)

$HAL_MASTER_SERVER $HAL_MASTER_DB -e "SELECT * FROM method_link_species_set WHERE method_link_species_set_id=${MLSS_HAL_ID};"

export COMPARA_HAL_DIR="/hps/nobackup/flicek/wormbase/parasite/cactus_alignments/WBPS18"

init_pipeline.pl Bio::EnsEMBL::Compara::PipeConfig::RegisterHALFile_conf \
    $(mysql-ps-prod-1-w details hive) \
    -name hal_$HAL_SUFFIX \
    -mlss_id $MLSS_HAL_ID \
    -reg_conf $PARASITE_CONF/compara.registry.pm \
    -division $HAL_COLLECTION_NAME \
    -species_name_mapping "$NAMING_MAPPING" \
    -master_db $HAL_MASTER_DB_URL