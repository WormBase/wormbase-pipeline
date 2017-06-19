#!/bin/bash 
set -e

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
    printf "Run the xref pipeline for 1 core database.\nMandatory arguments: -taxon, -name (your user id ex: ms41), -alias (species_name_bioproject ex:steinernema_glaseri_prjna204943), -ensembl_vers, -parasite_vers, -host, -pass, -user, -port (connection details for the server where your core db is and where the xref dbs will be created).\n"
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

MYTMPDIR=/nfs/nobackup/ensemblgenomes/wormbase/parasite/xref
ENSEMBL_ROOT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ ! -d "${ENSEMBL_ROOT_DIR}/ensembl" ]; then
cd ${ENSEMBL_ROOT_DIR}
git clone https://github.com/Ensembl/ensembl.git
fi


#2) parsing step
#a) append new species to xref_config.ini (if not already done)

printf '%s\n' '-----PARSING STEP-----'

species=`echo $alias | cut -d'_' -f 1,2`
dbname=${name}_${alias}_xref_${parasite_vers}_${ensembl_vers}_1

if ! grep -q $species ${ENSEMBL_ROOT_DIR}/ensembl/misc-scripts/xref_mapping/xref_config.ini; then
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
source          = UniParc::MULTI" >> ${ENSEMBL_ROOT_DIR}/ensembl/misc-scripts/xref_mapping/xref_config.ini
fi

printf "Config file ensembl/misc-scripts/xref_mapping/xref_config.ini updated.\n"

#b) convert the configuration into a database
cd ${ENSEMBL_ROOT_DIR}/ensembl/misc-scripts/xref_mapping/
perl xref_config2sql.pl > sql/populate_metadata.sql

printf "Config table created.\n"

#c) parsing now

printf "We will now start parsing the sources \n"
printf "The output will be written to ${MYTMPDIR}/${species}/${alias}_PARSER1.OUT\n"

if [ ! -d "${MYTMPDIR}/${species}/input_data" ]; then
cd ${MYTMPDIR}
mkdir ${species}
cd ${species}
mkdir input_data
fi

cd ${ENSEMBL_ROOT_DIR}/ensembl/misc-scripts/xref_mapping/
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
 -download_dir ${MYTMPDIR}/${species}/input_data \
  &> ${MYTMPDIR}/${species}/${alias}_PARSER1.OUT

#d) xref database back up before mapping 1
printf "Backing up the xref database..\n"
printf "The xref database will be dumped in ${MYTMPDIR}/${species}/sql_dump/${alias}_xref_after_parsing.sql\n"

if [ ! -d "${MYTMPDIR}/${species}/sql_dump" ]; then
cd ${MYTMPDIR}/${species}
mkdir sql_dump
fi

mysqldump -P$port  -h$host  -u$user -p$pass $dbname  > ${MYTMPDIR}/${species}/sql_dump/${alias}_xref_after_parsing.sql 

#3) mapping step 1
#a) create config file
printf '%s\n' '-----MAPPING STEP 1-----'

if [ ! -d "${MYTMPDIR}/${species}/mapping" ]; then
cd ${MYTMPDIR}/${species}
mkdir mapping
cd mapping
mkdir xref
mkdir core
fi
 
cd ${ENSEMBL_ROOT_DIR}/ensembl/misc-scripts/xref_mapping/
coredb=${name}_${alias}_core_${parasite_vers}_${ensembl_vers}_1

if [ ! -f "${alias}_xref_mapper.input"  ];
then
echo "xref
host=$host
port=$port
dbname=$dbname
user=$user
password=$pass
dir=${MYTMPDIR}/${species}/mapping/xref

species=$species
taxon=parasite
host=$host
port=$port
dbname=$coredb
user=$user
password=$pass
dir=${MYTMPDIR}/${species}/mapping/core

farm
queue=production-rh7
exonerate=/nfs/panda/ensemblgenomes/external/exonerate-2/bin/exonerate" > ${alias}_xref_mapper.input
fi

printf "Config file ensembl/misc-scripts/xref_mapping/${alias}_xref_mapper.input created.\n"

#b) mapping now

printf "We will now start mapping phase 1.\n"
printf "The output will be written to ${MYTMPDIR}/${species}/${alias}_MAPPER1.OUT\n"

perl xref_mapper.pl -file ${ENSEMBL_ROOT_DIR}/ensembl/misc-scripts/xref_mapping/${alias}_xref_mapper.input -dumpcheck >& ${MYTMPDIR}/${species}/${alias}_MAPPER1.out

#d) xref database and core database back up before mapping 2
printf "Backing up the xref database..\n"
printf "The xref database will be dumped in ${MYTMPDIR}/${species}/sql_dump/${alias}_xref_after_mapping1.sql\n"
mysqldump -P$port  -h$host -u$user -p$pass  $dbname  > ${MYTMPDIR}/${species}/sql_dump/${alias}_xref_after_mapping1.sql

printf "Backing up the core database..\n"
printf "The core database will be dumped in ${MYTMPDIR}/${species}/sql_dump/${alias}_core_after_mapping1.sql\n"
mysqldump  -P$port  -h$host -u$user -p$pass  $coredb  > ${MYTMPDIR}/${species}/sql_dump/${alias}_core_after_mapping1.sql

#3) mapping step 2

printf "We will now start mapping phase 2.\n"
printf "The output will be written to ${MYTMPDIR}/${species}/${alias}_MAPPER2.OUT\n"
perl xref_mapper.pl -file ${ENSEMBL_ROOT_DIR}/ensembl/misc-scripts/xref_mapping/${alias}_xref_mapper.input -upload >& ${MYTMPDIR}/${species}/${alias}_MAPPER2.out

