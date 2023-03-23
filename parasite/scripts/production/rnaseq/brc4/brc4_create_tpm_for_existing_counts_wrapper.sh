brc4_component_dir=$1
work_dir=$2
create_tpm_script=${WORM_CODE}/parasite/scripts/production/rnaseq/brc4/brc4_create_tpm_for_existing_counts.sh

ls -1 $brc4_component_dir | while read genome; do
  echo $genome
  genome_dir=${brc4_component_dir}/${genome}
  fl_file=${work_dir}/index/${genome}/${genome}.feature_length_per_gene.tsv
  log_dir=${work_dir}/recount_log;
  mkdir -p $log_dir
  if [ ! -f "$fl_file" ]; then echo "$fl_file does not exist. Exiting"; exit 1; fi;
  ls -1 $genome_dir | while read run_index; do
    study=$(echo $run_index | cut -d"_" -f1);
    run=$(echo $run_index | cut -d"_" -f2);
    run_dir=${genome_dir}/${run_index}/${run};
    log_prefix=${log_dir}/${run_index}_recount;
    ls -1 ${run_dir}/*.counts | grep -v "first" | grep -v "second" | while read sum_htseq_file; do
        echo "Creating TPM file for ${sum_htseq_file}";
        bash ${create_tpm_script} ${fl_file} ${sum_htseq_file};
        echo ""
      done
    printf "\n\n";
  done;
done



#genome=steinernema_carpocapsae_v1prjna202318
#
#fl_file=${work_dir}/index/${genome}/${genome}.feature_length_per_gene.tsv
#sum_htseq_file=genes.htseq-union.unstranded.counts