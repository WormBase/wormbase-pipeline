#!/bin/bash 
set -e

# Set some defaults - we want to write the xref database into the prod db
# Source ENSEMBL_VERSION and PARASITE_VERSION from the environment - populated via `module load`

# The script also has the options to provide commands itself, you could integrate it differently. Maybetaging_host
# ${PARASITE_STAGING_MYSQL}-ensrw details suffix_END | sed 's/--[a-z]*END//g'
# Also gives host port user pass.
pass_file_for_prod_db=/nfs/panda/ensemblgenomes/external/mysql-cmds/ensrw/mysql-ps-prod
if ! [ -r $pass_file_for_prod_db ]; then echo "The db password is at $pass_file_for_prod_db. But you can not read it - such shame."; exit 1; fi
db_string=$( perl -ne 'print $1 if /run_mysql\(qw\((.*)\)\)/' < $pass_file_for_prod_db ) 
read host port user pass <<< $db_string
ensembl_vers=${ENSEMBL_VERSION}
parasite_vers=${PARASITE_VERSION}
name=$(whoami)
pass_file_for_staging_db=$( which ${PARASITE_STAGING_MYSQL}-ensrw )
if ! [ -r $pass_file_for_staging_db ]; then echo "You tried to get pass file from staging DB using an env var etc. and LOL it didn't work it told you $pass_file_for_staging_db" ; exit 1 ; fi

staging_db_string=$( perl -ne 'print $1 if /run_mysql\(qw\((.*)\)\)/' < $pass_file_for_staging_db )
read staging_host staging_port staging_user staging_pass <<< $staging_db_string
if ! [ "$staging_host" -a "$staging_port" -a "$staging_user" -a "$staging_pass" ]; then echo "Could not get staging DB credentials from string: $staging_db_string" ; exit 1 ; fi 

#1) config
while [[ $# > 0 ]]
do
key="$1"

case $key in
    -taxon)
    taxon="$2"
    shift # past argument
    ;;
    -alias)
    alias="$2"
    shift # past argument
    ;;
    -help)
    printf "Run the xref pipeline for 1 core database.\nMandatory arguments: -taxon, -alias (group_species_bioproject ex:steinernema_glaseri_prjna204943).\n"
    exit 1
    ;;
    *)
        printf "Invalid argument! \n"
	exit 1
    ;;
esac
shift # past argument or value
done

if ! [ ${alias} ]; then echo "alias is unset"; exit 1; fi;
if ! [ ${taxon} ]; then echo "taxon is unset"; exit 1; fi

XREF_TMP_DIR=${XREF_TMP_DIR:-/nfs/nobackup/ensemblgenomes/wormbase/parasite/xref/rel_$PARASITE_VERSION}

# we want to modify xref_config.ini in Ensembl's code so get the local copy
ENSEMBL_ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/ensembl"
rsync -a --delete "$ENSEMBL_CVS_ROOT_DIR/ensembl" "${ENSEMBL_ROOT_DIR}"

#2) parsing step
#a) append new species to xref_config.ini (if not already done)

printf '%s\n' '-----PARSING STEP-----'

species=`echo $alias | cut -d'_' -f 1,2`
dbname=${name}_${alias}_xref_${parasite_vers}_${ensembl_vers}_1

if ! grep -q $species ${ENSEMBL_ROOT_DIR}/misc-scripts/xref_mapping/xref_config.ini; then
echo "
[species $species]
taxonomy_id     = $taxon
aliases         = $alias
source          = EntrezGene::MULTI
source          = GO::MULTI
source          = RefSeq_dna::MULTI-invertebrate
source          = RefSeq_peptide::MULTI-invertebrate
source          = Uniprot/SPTREMBL::MULTI-invertebrate
source          = Uniprot/SWISSPROT::MULTI-invertebrate
source          = UniParc::MULTI" >> ${ENSEMBL_ROOT_DIR}/misc-scripts/xref_mapping/xref_config.ini
fi

printf "Config file ensembl/misc-scripts/xref_mapping/xref_config.ini updated.\n"

#b) convert the configuration into a database
cd ${ENSEMBL_ROOT_DIR}/misc-scripts/xref_mapping/
perl xref_config2sql.pl > sql/populate_metadata.sql

printf "Config table created.\n"

#c) parsing now

printf "We will now start parsing the sources \n"
printf "The output will be written to ${XREF_TMP_DIR}/${species}/${alias}_PARSER1.out\n"

mkdir -p "${XREF_TMP_DIR}/input_data" # Shared between species - all our worms need the same reference data anyway
mkdir -p "${XREF_TMP_DIR}/${species}/sql_dump"
#Not sure what will go inside these two - Wojtek
mkdir -p "${XREF_TMP_DIR}/${species}/mapping/xref"
mkdir -p "${XREF_TMP_DIR}/${species}/mapping/core"

cd ${ENSEMBL_ROOT_DIR}/misc-scripts/xref_mapping/
[ "$XREF_PARSER_SKIP" ] || perl xref_parser.pl \
  -user $user \
  -pass $pass \
  -host $host \
  -port $port \
  -species $species \
  -drop_db \
  -create \
  -dbname $dbname \
  -checkdownload \
  -stats \
 -download_dir ${XREF_TMP_DIR}/input_data \
  &> ${XREF_TMP_DIR}/${species}/${alias}_PARSER1.out

#d) xref database back up before mapping 1
printf "Backing up the xref database..\n"
printf "The xref database will be dumped in ${XREF_TMP_DIR}/${species}/sql_dump/${alias}_xref_after_parsing.sql\n"


mysqldump -P$port  -h$host  -u$user -p$pass $dbname  > ${XREF_TMP_DIR}/${species}/sql_dump/${alias}_xref_after_parsing.sql 

#3) mapping step 1
#a) create config file
printf '%s\n' '-----MAPPING STEP 1-----'

 
coredb=${alias}_core_${parasite_vers}_${ensembl_vers}_1

echo "xref
host=$host
port=$port
dbname=$dbname
user=$user
password=$pass
dir=${XREF_TMP_DIR}/${species}/mapping/xref

species=$species
taxon=parasite
host=$staging_host
port=$staging_port
dbname=$coredb
user=$staging_user
password=$staging_pass
dir=${XREF_TMP_DIR}/${species}/mapping/core

farm
queue=production-rh7
exonerate=/nfs/panda/ensemblgenomes/external/exonerate-2/bin/exonerate" > ${ENSEMBL_ROOT_DIR}/misc-scripts/xref_mapping/${alias}_xref_mapper.input

printf "Config file ensembl/misc-scripts/xref_mapping/${alias}_xref_mapper.input created.\n"

#b) mapping now

printf "We will now start mapping phase 1.\n"
printf "The output will be written to ${XREF_TMP_DIR}/${species}/${alias}_MAPPER1.out\n"

perl ${ENSEMBL_ROOT_DIR}/misc-scripts/xref_mapping/xref_mapper.pl -file ${ENSEMBL_ROOT_DIR}/misc-scripts/xref_mapping/${alias}_xref_mapper.input -dumpcheck >& ${XREF_TMP_DIR}/${species}/${alias}_MAPPER1.out

#d) xref database and core database back up before mapping 2
printf "Backing up the xref database..\n"
printf "The xref database will be dumped in ${XREF_TMP_DIR}/${species}/sql_dump/${alias}_xref_after_mapping1.sql\n"
mysqldump -P$port  -h$host -u$user -p$pass  $dbname  > ${XREF_TMP_DIR}/${species}/sql_dump/${alias}_xref_after_mapping1.sql

printf "Backing up the core database..\n"
printf "The core database will be dumped in ${XREF_TMP_DIR}/${species}/sql_dump/${alias}_core_after_mapping1.sql\n"
mysqldump  -P$staging_port  -h$staging_host -u$staging_user -p$staging_pass  $coredb  > ${XREF_TMP_DIR}/${species}/sql_dump/${alias}_core_after_mapping1.sql

#3) mapping step 2

printf "We will now start mapping phase 2.\n"
printf "The output will be written to ${XREF_TMP_DIR}/${species}/${alias}_MAPPER2.out\n"
perl ${ENSEMBL_ROOT_DIR}/misc-scripts/xref_mapping/xref_mapper.pl -file ${ENSEMBL_ROOT_DIR}/misc-scripts/xref_mapping/${alias}_xref_mapper.input -upload >& ${XREF_TMP_DIR}/${species}/${alias}_MAPPER2.out

