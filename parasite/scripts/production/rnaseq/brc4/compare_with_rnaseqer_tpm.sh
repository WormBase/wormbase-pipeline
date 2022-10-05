RESULTS_DIR=$1
RNASEQER_FTP="/nfs/ftp/public/databases/arrayexpress/data/atlas/rnaseq/"

for study in $(ls $RESULTS_DIR -1);
  do
  short_study=${study:0:6}
  study_dir=${RESULTS_DIR}/${study};
  # elength_file=$(ls ${RESULTS_DIR}/../../../temp/*/index/*/*.exon-lengths.tsv | head -n1)
  echo "study: "${study}
  # echo "   summing up metazoa pipeline output..."
  # fsnu=${study_dir}/genes.htseq-union.firststrand.nonunique.counts
  # ssnu=${study_dir}/genes.htseq-union.secondstrand.nonunique.counts
  # unuc=${study_dir}/genes.htseq-union.unstranded.counts
  #  cat ${study_dir}/genes.htseq-union.firststrand.nonunique.counts ${study_dir}/genes.htseq-union.secondstrand.nonunique.counts | awk '{a[$1]+=$2}END{for(i in a) print i,a[i]}' > ${study_dir}/genes.htseq-union.nonunique.counts;
  #  cat ${study_dir}/genes.htseq-union.firststrand.counts ${study_dir}/genes.htseq-union.secondstrand.counts | awk '{a[$1]+=$2}END{for(i in a) print i,a[i]}' > ${study_dir}/genes.htseq-union.counts;
  # echo "   create metazoa pipeline tpm file..."
  # join <(sort -k1 ${unuc}) <(sort -k1 ${elength_file}) | awk '{print $1"\t"$2/($3/1000)}' > ${unuc}.rpk
  # million_factor=$(awk '{Total=Total+$2} END{print Total/1000000}' ${unuc}.rpk)
  # awk -v million_factor="${million_factor}" '{print $1"\t"$2/million_factor}' ${unuc}.rpk > ${unuc}.tpm
  echo "   copying rnaseqer tpm file..."
  cp ${RNASEQER_FTP}/${short_study}/${study}/${study}.pe.genes.tpm.htseq2.irap.tsv ${study_dir}/;
  echo "   performing comparison and printing comparison plots..."
  export PATH=${PARASITE_PRODUCTION}/data/conda_envs/wbps-expression/bin:$PATH;
  Rscript ${WORM_CODE}/parasite/scripts/production/rnaseq/brc4/comparison_plots_tpm.R ${study_dir} ${study};
  echo "   Done"
  echo ""
  done


