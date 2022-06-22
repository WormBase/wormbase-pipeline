#!/bin/bash

# Everything below will go to the file 'log.out':
# Input parsing
INP1_FILE=$1
INP1_NAME=$2
INP2_FILE=$3
INP2_NAME=$4
PROJECT_NAME=$5
WORKING_DIR=${PARASITE_SCRATCH}/reciprocal_blast/${PROJECT_NAME}

echo "----------Info: Welcome to Dio's reciprocal best blast hit!------------"
echo '-----------------------------------------------------------------------'

printf "$(date '+%d/%m/%Y %H:%M:%S') - Info: Given Input:\n"
printf "PROJECT NAME: ${PROJECT_NAME}\n"
printf "INPUT FILE 1: ${INP1_FILE}\n"
printf "INPUT NAME 1: ${INP1_NAME}\n"
printf "INPUT FILE 2: ${INP2_FILE}\n"
printf "INPUT NAME 2: ${INP2_NAME}\n"
printf "WORKING DIRECTORY: ${WORKING_DIR}\n"

## To dump prot fa from a core db you can use:
# perl $WORM_CODE/scripts/ENSEMBL/scripts/dump_proteins.pl -canonical_only --host=$DBHOSTNAME --port=$DBPORT --user=$DBUSER --dbname=$core_db --outfile=$species.prot.fa

# Check blast is there
if ! command -v blastp &> /dev/null
then
    printf "$(date '+%d/%m/%Y %H:%M:%S') - Error: BLASTP software could not be detected or it's not executable:\n"
    exit
fi

if ! command -v makeblastdb &> /dev/null
then
    printf "$(date '+%d/%m/%Y %H:%M:%S') - Error: makeblastdb software could not be detected or it's not executable:\n"
    exit
fi
printf "$(date '+%d/%m/%Y %H:%M:%S') - Info: BLAST software has been detected and it's executable:\n"
command -v blastp
command -v makeblastdb

printf "$(date '+%d/%m/%Y %H:%M:%S') - Info: Creating Input/Output Directories...\n"
# Directories handling
INPUT_DIR=${WORKING_DIR}/input
OUTPUT_DIR=${WORKING_DIR}/output
INP1_INPUT=${INPUT_DIR}/${INP1_NAME}
INP2_INPUT=${INPUT_DIR}/${INP2_NAME}
mkdir -p ${WORKING_DIR}
mkdir -p ${INP1_INPUT}
mkdir -p ${INP2_INPUT}
mkdir -p ${OUTPUT_DIR}
printf "$(date '+%d/%m/%Y %H:%M:%S') - Info: Input/Output Directories have been created.\n"

# Copy input files
printf "$(date '+%d/%m/%Y %H:%M:%S') - Info: Copying and renaming input fasta files...\n"
INP1_PROT=${INP1_INPUT}/${INP1_NAME}.prot.fa
INP2_PROT=${INP2_INPUT}/${INP2_NAME}.prot.fa
rsync -avz ${INP1_FILE} ${INP1_PROT}
rsync -avz ${INP2_FILE} ${INP2_PROT}
printf "$(date '+%d/%m/%Y %H:%M:%S') - Info: Fasta files have been retrieved.\n"

# makeblastdb for both inputs
printf "$(date '+%d/%m/%Y %H:%M:%S') - Info: Running makeblastdb for both input fasta files...\n"
INP1_DB=${INP1_INPUT}/${INP1_NAME}
INP2_DB=${INP2_INPUT}/${INP2_NAME}
makeblastdb -in ${INP1_PROT} -dbtype prot -out ${INP1_DB}
makeblastdb -in ${INP2_PROT} -dbtype prot -out ${INP2_DB}
printf "$(date '+%d/%m/%Y %H:%M:%S') - Info: Blast DBs have been created for both files.\n"

printf "$(date '+%d/%m/%Y %H:%M:%S') - Info: Beginning reciprocal BLAST...\n"
printf "$(date '+%d/%m/%Y %H:%M:%S') - Info: Running BLASTP for ${INP1_PROT} to ${INP2_DB} DB...\n"
blastp -query "${INP1_PROT}" -db "${INP2_DB}" -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > ${OUTPUT_DIR}/${INP1_NAME}_to_${INP2_NAME}.blast.outfmt6
printf "$(date '+%d/%m/%Y %H:%M:%S') - Info: Running BLASTP for ${INP2_PROT} to ${INP1_DB} DB...\n"
blastp -query "${INP2_PROT}" -db "${INP1_DB}" -max_target_seqs 1 -outfmt 6 -evalue 1e-3 -num_threads 8 > ${OUTPUT_DIR}/${INP2_NAME}_to_${INP1_NAME}.blast.outfmt6
printf "$(date '+%d/%m/%Y %H:%M:%S') - Info: Done Blasting.\n"

printf "$(date '+%d/%m/%Y %H:%M:%S') - Info: Looking for reciprocal best blast hits,..\n"
OUT_RBH="${OUTPUT_DIR}/blast_RBH.txt"
OUT_SUMMARY="${OUTPUT_DIR}/blast_RBH_summary.txt"

rm $OUT_RBH
rm $OUT_SUMMARY

printf "queryHit\tdbHit\n" > $OUT_RBH
printf "queryHits\tdbHits\tbestHits\n" > $OUT_SUMMARY

while IFS=$'\t' read -r f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12;
  do
    #Determine RBH to DB blast results
    if grep -q "$f2"$'\t'"$f1"$'\t' ${OUTPUT_DIR}/${INP2_NAME}_to_${INP1_NAME}.blast.outfmt6; then
      printf "$f1\t$f2\n" >> $OUT_RBH
    fi;
  done < ${OUTPUT_DIR}/${INP1_NAME}_to_${INP2_NAME}.blast.outfmt6 | head

printf "$(date '+%d/%m/%Y %H:%M:%S') - Info: Gathering Statistics...\n"
queryHits=$(wc -l "${OUTPUT_DIR}/${INP1_NAME}_to_${INP2_NAME}.blast.outfmt6" | cut -d ' ' -f 1)
dbHits=$(wc -l "${OUTPUT_DIR}/${INP2_NAME}_to_${INP1_NAME}.blast.outfmt6" | cut -d ' ' -f 1)
bestHits=$(($(wc -l "$OUT_RBH" | cut -d ' ' -f 1)-1))
printf "${queryHits}\t${dbHits}\t${bestHits}" >> $OUT_SUMMARY

#Output end status message
printf "$(date '+%d/%m/%Y %H:%M:%S') - Info: Finished processing...\n"
printf "Statistics:\n"
cat $OUT_OUT_SUMMARY
printf "\n$(date '+%d/%m/%Y %H:%M:%S') - Info: Output Files...\n"
printf "OUTPUT RBH: ${OUT_RBH}\n"
printf "OUTPUT SUMMARY STATISTICS: ${OUT_SUMMARY}\n"



