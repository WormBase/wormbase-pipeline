FTP_CHECK_FILE=$1

#check if file exists
if [[ ! -f ${FTP_CHECK_FILE} ]];
then
    echo "${FTP_CHECK_FILE} does not exist. Exiting";
    exit 1;
fi

grep -v "^working" ${FTP_CHECK_FILE} | while read eline;
do
  if echo $eline | grep -q "genomesoftmaskedDNA_in_dump:";
  then
    echo "softmask error: $eline";
    db=$(echo $eline | cut -d' ' -f1);
    #echo $db;
    ftp_path=$(perl -MProductionMysql -e "print ProductionMysql->staging->core_db_to_local_ftp_path_n_filename(${db}).\"\n\";");
    echo "submitting job to fix it:";
    outfile=${PARASITE_SCRATCH}/dumps/WBPS${PARASITE_VERSION}/FTP/${ftp_path}.genomic_softmasked.fa
    MEM_MB=8000
    bsub -q production -M${MEM_MB} -R "select[mem > ${MEM_MB}] rusage[mem=${MEM_MB}]" \
    "perl ${WORM_CODE}/scripts/ENSEMBL/scripts/dump_genome.pl $($PARASITE_STAGING_MYSQL-w details script) --softmask --dbname ${db} --outfile /homes/digri/dio.fa; gzip -n -f /homes/digri/dio.fa;"
  fi;
    if echo $eline | grep -q "genomeDNA_in_dump:";
  then
    echo "genome fasta error: $eline";
    db=$(echo $eline | cut -d' ' -f1);
    #echo $db;
    ftp_path=$(perl -MProductionMysql -e "print ProductionMysql->staging->core_db_to_local_ftp_path_n_filename(${db}).\"\n\";");
    echo "submitting job to fix it:";
    outfile=${PARASITE_SCRATCH}/dumps/WBPS${PARASITE_VERSION}/FTP/${ftp_path}.genomic.fa
    MEM_MB=8000
    bsub -q production -M${MEM_MB} -R "select[mem > ${MEM_MB}] rusage[mem=${MEM_MB}]" \
    "perl ${WORM_CODE}/scripts/ENSEMBL/scripts/dump_genome.pl $($PARASITE_STAGING_MYSQL-w details script) --dbname ${db} --outfile /homes/digri/dio.fa; gzip -n -f /homes/digri/dio.fa;"
  fi;
done



