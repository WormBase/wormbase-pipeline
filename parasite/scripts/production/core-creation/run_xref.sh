#!/bin/bash 
set -e

# Set some defaults - we want to write the xref database into the prod db
# Source ENSEMBL_VERSION and PARASITE_VERSION from the environment - populated via `module load`
pass_file_for_prod_db=/nfs/panda/ensemblgenomes/external/mysql-cmds/ensrw/mysql-ps-prod
if ! [ -r $pass_file_for_prod_db ]; then echo "The db password is at $pass_file_for_prod_db. But you can not read it - such shame."; exit 1; fi
db_string=$( perl -ne 'print $1 if /run_mysql\(qw\((.*)\)\)/' < $pass_file_for_prod_db ) 
read host port user pass <<< $db_string
ensembl_vers=${ENSEMBL_VERSION}
parasite_vers=${PARASITE_VERSION}
name=$(whoami)

#1) config
while [[ $# > 0 ]]
do
key="$1"

case $key in
    -taxon)
    taxon="$2"
    shift # past argument
    ;;
    -name)
    name="$2"
    shift # past argument
    ;;
    -alias)
    alias="$2"
    shift # past argument
    ;;
    -ensembl_vers)
    ensembl_vers="$2"
    shift
    ;;
    -parasite_vers)
    parasite_vers="$2"
    shift
    ;;
    -host)
    host="$2"
    shift
    ;;
    -port)
    port="$2"
    shift
    ;;
    -pass)
    pass="$2"
    shift
    ;;
    -user)
    user="$2"
    shift
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

if ! [ ${name} ]; then echo "name is unset"; exit 1; fi
if ! [ ${alias} ]; then echo "alias is unset"; exit 1; fi;
if ! [ ${taxon} ]; then echo "taxon is unset"; exit 1; fi
if ! [ ${ensembl_vers} ]; then echo "ensembl_vers is unset"; exit 1; fi;
if ! [ ${parasite_vers} ]; then echo "parasite_vers is unset"; exit 1; fi;
if ! [ ${host} ]; then echo "host is unset"; exit 1; fi
if ! [ ${port} ]; then echo "port is unset"; exit 1; fi
if ! [ ${user} ]; then echo "user is unset"; exit 1; fi
if ! [ ${pass} ]; then echo "pass is unset"; exit 1; fi

XREF_TMP_DIR=${XREF_TMP_DIR:-/nfs/nobackup/ensemblgenomes/wormbase/parasite/xref}

# we want to modify xref_config.ini in Ensembl's code so get the local copy
ENSEMBL_ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/ensembl"
[ -d "${ENSEMBL_ROOT_DIR}" ] || git clone -q "$ENSEMBL_CVS_ROOT_DIR/ensembl" "${ENSEMBL_ROOT_DIR}"


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
printf "The output will be written to ${XREF_TMP_DIR}/${species}/${alias}_PARSER1.OUT\n"

mkdir -p "${XREF_TMP_DIR}/${species}/input_data"

cd ${ENSEMBL_ROOT_DIR}/misc-scripts/xref_mapping/
perl xref_parser.pl \
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
 -download_dir ${XREF_TMP_DIR}/${species}/input_data \
  &> ${XREF_TMP_DIR}/${species}/${alias}_PARSER1.OUT

#d) xref database back up before mapping 1
printf "Backing up the xref database..\n"
printf "The xref database will be dumped in ${XREF_TMP_DIR}/${species}/sql_dump/${alias}_xref_after_parsing.sql\n"

if [ ! -d "${XREF_TMP_DIR}/${species}/sql_dump" ]; then
cd ${XREF_TMP_DIR}/${species}
mkdir sql_dump
fi

mysqldump -P$port  -h$host  -u$user -p$pass $dbname  > ${XREF_TMP_DIR}/${species}/sql_dump/${alias}_xref_after_parsing.sql 

#3) mapping step 1
#a) create config file
printf '%s\n' '-----MAPPING STEP 1-----'

if [ ! -d "${XREF_TMP_DIR}/${species}/mapping" ]; then
cd ${XREF_TMP_DIR}/${species}
mkdir mapping
cd mapping
mkdir xref
mkdir core
fi
 
cd ${ENSEMBL_ROOT_DIR}/misc-scripts/xref_mapping/
coredb=${name}_${alias}_core_${parasite_vers}_${ensembl_vers}_1

if [ ! -f "${alias}_xref_mapper.input"  ];
then
echo "xref
host=$host
port=$port
dbname=$dbname
user=$user
password=$pass
dir=${XREF_TMP_DIR}/${species}/mapping/xref

species=$species
taxon=parasite
host=$host
port=$port
dbname=$coredb
user=$user
password=$pass
dir=${XREF_TMP_DIR}/${species}/mapping/core

farm
queue=production-rh7
exonerate=/nfs/panda/ensemblgenomes/external/exonerate-2/bin/exonerate" > ${alias}_xref_mapper.input
fi

printf "Config file ensembl/misc-scripts/xref_mapping/${alias}_xref_mapper.input created.\n"

#b) mapping now

printf "We will now start mapping phase 1.\n"
printf "The output will be written to ${XREF_TMP_DIR}/${species}/${alias}_MAPPER1.OUT\n"

perl xref_mapper.pl -file ${ENSEMBL_ROOT_DIR}/misc-scripts/xref_mapping/${alias}_xref_mapper.input -dumpcheck >& ${XREF_TMP_DIR}/${species}/${alias}_MAPPER1.out

#d) xref database and core database back up before mapping 2
printf "Backing up the xref database..\n"
printf "The xref database will be dumped in ${XREF_TMP_DIR}/${species}/sql_dump/${alias}_xref_after_mapping1.sql\n"
mysqldump -P$port  -h$host -u$user -p$pass  $dbname  > ${XREF_TMP_DIR}/${species}/sql_dump/${alias}_xref_after_mapping1.sql

printf "Backing up the core database..\n"
printf "The core database will be dumped in ${XREF_TMP_DIR}/${species}/sql_dump/${alias}_core_after_mapping1.sql\n"
mysqldump  -P$port  -h$host -u$user -p$pass  $coredb  > ${XREF_TMP_DIR}/${species}/sql_dump/${alias}_core_after_mapping1.sql

#3) mapping step 2

printf "We will now start mapping phase 2.\n"
printf "The output will be written to ${XREF_TMP_DIR}/${species}/${alias}_MAPPER2.OUT\n"
perl xref_mapper.pl -file ${ENSEMBL_ROOT_DIR}/misc-scripts/xref_mapping/${alias}_xref_mapper.input -upload >& ${XREF_TMP_DIR}/${species}/${alias}_MAPPER2.out

