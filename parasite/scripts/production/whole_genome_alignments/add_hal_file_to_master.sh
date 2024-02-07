export HAL_COMPARA_SERVER=$PARASITE_STAGING_MYSQL
export HAL_COMPARA_DB=ensembl_compara_parasite_${PARASITE_VERSION}_${ENSEMBL_VERSION}
export HAL_MASTER_SERVER=mysql-ps-prod-1
export HAL_MASTER_DB=ensembl_compara_master_parasite
export PROD_DB=mysql-ps-prod-1

export ENSEMBL_ROOT_DIR=${ENSEMBL_CVS_ROOT_DIR}
export HALXS_DIR=$ENSEMBL_ROOT_DIR/ensembl-compara/src/perl/modules/Bio/EnsEMBL/Compara/HAL/HALXS
export HALXS_TEST=${HALXS_DIR}/HALXS.c

export HAL_COMPARA_DB_URL=$($HAL_COMPARA_SERVER-w details url $HAL_COMPARA_DB)
export HAL_MASTER_DB_URL=$($HAL_MASTER_SERVER-w details url $HAL_MASTER_DB)

# perl script extracting genomes from hal
PERL_SCRIPT=${WORM_CODE}/parasite/scripts/production/whole_genome_alignments/dump_hal_genome_ids.pl

usage() {
    echo "Usage: $0 -e <HAL_suffix_file> -h <final_HAL_file>"
    echo "  -e  A line-delimited file containing the HAL suffixes given to the HAL names."
    echo "  -f  The final HAL file produced by the cactus pipeline."
    exit 1
}

while getopts ":e:f:" opt; do
    case $opt in
        e)
            HAL_SUFFIX_FILE=$OPTARG
            ;;
        f)
            FINAL_HAL_FILE=$OPTARG
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            usage
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            usage
            ;;
    esac
done

if [ -z "$HAL_SUFFIX_FILE" ] || [ -z "$FINAL_HAL_FILE" ]; then
    echo "Error: Both -e and -f options are required."
    usage
fi

if [ -f "$FINAL_HAL_FILE" ]; then
    echo "$(date +"%Y-%m-%d %H:%M:%S") Final HAL file: $FINAL_HAL_FILE"
    export COMPARA_HAL_DIR=$(dirname $FINAL_HAL_FILE)
    echo "$(date +"%Y-%m-%d %H:%M:%S") Whole genome alignments directory: $COMPARA_HAL_DIR"
else
    echo "${FINAL_HAL_FILE} file does not exist."
    usage
fi

OUTPUT_BEEKEEPER_COMMANDS_TXT=${COMPARA_HAL_DIR}/beekeeper_commands.txt

if [ -f "$HALXS_TEST" ]; then
    echo "$(date +"%Y-%m-%d %H:%M:%S") HALXS has been compiled."
else
    echo "$(date +"%Y-%m-%d %H:%M:%S") HALXS has not been compiled and we need to compile it now."
    cd $HALXS_DIR
    perl Makefile-Linuxbrew.PL && make
    cd $COMPARA_HAL_DIRs
fi


# Rest of the script logic using $HAL_SUFFIX_FILE and $FINAL_HAL_FILE variables

cat ${HAL_SUFFIX_FILE} | grep -v "#" |  while read HAL_SUFFIX; do 

    echo "$(date +"%Y-%m-%d %H:%M:%S") Working on ${HAL_SUFFIX}."
    
    export INITIAL_HAL_FILE=${PARASITE_SCRATCH}/cactus_alignments/WBPS${PARASITE_VERSION}/WBPS${PARASITE_VERSION}_${HAL_SUFFIX}.hal
    export HAL_FILE=${PARASITE_SCRATCH}/cactus_alignments/WBPS${PARASITE_VERSION}/multi/hal_files/WBPS${PARASITE_VERSION}_${HAL_SUFFIX}.hal

    echo "$(date +"%Y-%m-%d %H:%M:%S") Extracting $HAL_SUFFIX from $FINAL_HAL_FILE"
    halExtract --root $HAL_SUFFIX $FINAL_HAL_FILE ${INITIAL_HAL_FILE}; done | tee ${COMPARA_HAL_DIR}/split_hal_log_${HAL_SUFFIX}.txt

    if [ -f "$INITIAL_HAL_FILE" ]; then
        cp -i "$INITIAL_HAL_FILE" "$HAL_FILE"
    else
        echo "$(date +"%Y-%m-%d %H:%M:%S") ${INITIAL_HAL_FILE} file does not exist."
        exit 1;
    fi

    # Check if the hal file exists
    if [ -f "$HAL_FILE" ];
    then :;
    else echo "$(date +"%Y-%m-%d %H:%M:%S") File $HAL_FILE does not exist. Exiting."
    exit 1
    fi

    export HAL_COLLECTION_NAME=parasitehal-${HAL_SUFFIX}
    echo "$(date +"%Y-%m-%d %H:%M:%S") HAL collection name: parasitehal-${HAL_SUFFIX}"

    export HAL_MLSS_URL="#base_dir#/multi$(awk -F '/multi' '{print $2}' <<<$HAL_FILE)"


    # hal genomes
    echo "$(date +"%Y-%m-%d %H:%M:%S") Exporting genome names from $HAL_FILE to ${PARASITE_CONF}/compara.${HAL_SUFFIX}.species_list_for_hal"
    HAL_GENOMES=$(perl $PERL_SCRIPT --hal_file $HAL_FILE)
    HAL_GENOMES_COUNT=$(echo $HAL_GENOMES | tr " " "\n" | wc -l)
    echo $HAL_GENOMES | tr " " "\n" > ${PARASITE_CONF}/compara.${HAL_SUFFIX}.species_list_for_hal

    HAL_GENOME_DB_IDS=$(for halg in $HAL_GENOMES; do
        $HAL_MASTER_SERVER $HAL_MASTER_DB -Ne "SELECT genome_db_id FROM genome_db WHERE name LIKE '$halg';";
    done)
    HAL_GENOME_DB_IDS_COUNT=$(echo $HAL_GENOME_DB_IDS | tr " " "\n" | wc -l)
    HAL_GENOME_DB_ID=$(echo $HAL_GENOME_DB_IDS | tr " " ",")

    if [ "$HAL_GENOME_DB_IDS_COUNT" -ne "$HAL_GENOMES_COUNT" ];
        then
            echo "$(date +"%Y-%m-%d %H:%M:%S") WARNING: Not every genome in your HAL file has a genome_db in the genome_db table of your master database. This is expected as we have not ran compara for all the species in the alignment.";
            echo "$(date +"%Y-%m-%d %H:%M:%S") Updating master DB:";
            # add missing genomes to the genome_Db table and add assign them to collections
            perl ${WORM_CODE}/parasite/scripts/production/compara/parasite_update_compara_master_for_hal.pl \
                -reg_conf $PARASITE_CONF/compara.registry.pm \
                -collection $HAL_COLLECTION_NAME \
                -masterdbname $HAL_MASTER_DB \
                -taxdbname ncbi_taxonomy_parasite \
                -comparacode $ENSEMBL_CVS_ROOT_DIR/ensembl-compara \
                -sfile ${PARASITE_CONF}/compara.${HAL_SUFFIX}.species_list_for_hal \
                -locators
        else echo "$(date +"%Y-%m-%d %H:%M:%S") Every genome in your HAL file has a genome_db in the genome_db table. No need to update the master";
    fi


    # Datachecks
    HAL_GENOME_DB_IDS=$(\
    for halg in $HAL_GENOMES; do
        $HAL_MASTER_SERVER $HAL_MASTER_DB -Ne "SELECT genome_db_id FROM genome_db WHERE name LIKE '$halg';";
    done)
    HAL_GENOME_DB_IDS_COUNT=$(echo $HAL_GENOME_DB_IDS | tr " " "\n" | wc -l)
    HAL_GENOME_DB_ID=$(echo $HAL_GENOME_DB_IDS | tr " " ",")
    if [ "$HAL_GENOME_DB_IDS_COUNT" -ne "$HAL_GENOMES_COUNT" ];
        then
        echo "$(date +"%Y-%m-%d %H:%M:%S") ERROR: Not every genome in your HAL file has a genome_db in the genome_db table of your master database despite we ran parasite_update_compara_master_for_hal.pl. Exiting";
        exit 1;
    fi

    # Naming mapping:
    echo "$(date +"%Y-%m-%d %H:%M:%S") Obtaining naming mapping"
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
        then echo "$(date +"%Y-%m-%d %H:%M:%S") ERROR: Not every genome in your HAL file has a genome_db in the genome_db table of your master database."; exit 1;
    fi

    echo "$(date +"%Y-%m-%d %H:%M:%S") Creating MLSS for $HAL_FILE (collection-$HAL_COLLECTION_NAME)"
    perl ${ENSEMBL_ROOT_DIR}/ensembl-compara/scripts/pipeline/create_mlss.pl \
    --compara $HAL_MASTER_DB_URL --reg_conf $PARASITE_CONF/compara.registry.pm \
    --source wormbase --species_set_name collection-$HAL_COLLECTION_NAME \
    --method_link_type CACTUS_HAL --pw --genome_db_id ${HAL_GENOME_DB_ID} --url "'$HAL_MLSS_URL'" --release --force

    MLSS_HAL_ID=$($HAL_MASTER_SERVER $HAL_MASTER_DB -e "SELECT mlss.method_link_species_set_id AS '' FROM method_link_species_set AS mlss JOIN species_set_header AS ssh USING(species_set_id) WHERE ssh.name = 'collection-$HAL_COLLECTION_NAME'" | tail -1)
    re='^[0-9]+$'
    if ! [[ $MLSS_HAL_ID =~ $re ]]; then
        echo "$(date +"%Y-%m-%d %H:%M:%S") ERROR: MLSS_HAL_ID $MLSS_HAL_ID has not been set correctly. create_mlss did not run properly."
        exit 1
    fi
    echo "$(date +"%Y-%m-%d %H:%M:%S") collection-$HAL_COLLECTION_NAME) MLSS ID is ${MLSS_HAL_ID}"
    echo "$(date +"%Y-%m-%d %H:%M:%S") Testing database:"
    $HAL_MASTER_SERVER $HAL_MASTER_DB -e "SELECT * FROM method_link_species_set WHERE method_link_species_set_id=${MLSS_HAL_ID};"


    echo "$(date +"%Y-%m-%d %H:%M:%S") Fixing issue with quotes on method_link_species_set:"
    $HAL_MASTER_SERVER-w $HAL_MASTER_DB -Ne "UPDATE method_link_species_set SET url=\"$HAL_MLSS_URL\" WHERE url=\"'$HAL_MLSS_URL'\";"


    echo "$(date +"%Y-%m-%d %H:%M:%S") Creating RegisterHALFile pipeline for $HAL_FILE"
    init_pipeline.pl Bio::EnsEMBL::Compara::PipeConfig::RegisterHALFile_conf \
        $(${PROD_DB}-w details hive) \
        -name hal_$HAL_SUFFIX \
        -mlss_id $MLSS_HAL_ID \
        -reg_conf $PARASITE_CONF/compara.registry.pm \
        -division $HAL_COLLECTION_NAME \
        -species_name_mapping "$NAMING_MAPPING" \
        -master_db $HAL_MASTER_DB_URL
    
    
    hive_db=$(${PROD_DB} -Ne "SHOW DATABASES LIKE \"%${HAL_COLLECTION_NAME,,}%${ENSEMBL_VERSION}\";")

    # Count the number of lines in the output
    num_dbs=$(echo "$hive_db" | wc -l)

    if [ $num_dbs -eq 1 ]; then :;
    else
        # More than one database found, exit with an error
        echo "Error: Found zero or multiple databases matching the ${HAL_COLLECTION_NAME,,} pattern in ${PROD_DB}."
        exit 1
    fi

    hive_pipeline_url=$($PROD_DB-w details url $hive_db)
    echo "$(date +"%Y-%m-%d %H:%M:%S") RegisterHALFile pipeline hive db url: $hive_pipeline_url"

    echo "$(date +"%Y-%m-%d %H:%M:%S") Giving more memory to generate_pairwise_coverage_stats analysis:"
    tweak_pipeline.pl -url $hive_pipeline_url -tweak 'analysis[generate_pairwise_coverage_stats].resource_class=32Gb_job'

    echo "$(date +"%Y-%m-%d %H:%M:%S") Finished working on $HAL_SUFFIX"
    printf "#${HAL_COLLECTION_NAME} beekeeper command:\nexport COMPARA_HAL_DIR=${COMPARA_HAL_DIR}; beekeeper.pl -url $hive_pipeline_url -loop_until NO_WORK\n\n" >> ${OUTPUT_BEEKEEPER_COMMANDS_TXT}
done

printf "\n\n\n---------------------FINISHED WORKING--------------------\n\n\n"
echo "Beekeeper commands have been stored in $OUTPUT_BEEKEEPER_COMMANDS_TXT. Please run one pipeline at a time."
echo "Done"