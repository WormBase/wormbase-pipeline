RESULTS_DIR=$1
RNASEQER_FTP="/nfs/ftp/public/databases/arrayexpress/data/atlas/rnaseq/"


for study in $(ls $RESULTS_DIR -1);
  do
  short_study=${study:0:6}
  study_dir=${RESULTS_DIR}/${study};
  echo "study: "${study}
  echo "   summing up metazoa pipeline output..."
  fsnu=${study_dir}/genes.htseq-union.firststrand.nonunique.counts
  ssnu=${study_dir}/genes.htseq-union.secondstrand.nonunique.counts
  usnu=${study_dir}/genes.htseq-union.unstranded.nonunique.counts
  #  cat ${study_dir}/genes.htseq-union.firststrand.nonunique.counts ${study_dir}/genes.htseq-union.secondstrand.nonunique.counts | awk '{a[$1]+=$2}END{for(i in a) print i,a[i]}' > ${study_dir}/genes.htseq-union.nonunique.counts;
  #  cat ${study_dir}/genes.htseq-union.firststrand.counts ${study_dir}/genes.htseq-union.secondstrand.counts | awk '{a[$1]+=$2}END{for(i in a) print i,a[i]}' > ${study_dir}/genes.htseq-union.counts;
  #  echo "   copying rnaseqer file..."
  cp ${RNASEQER_FTP}/${short_study}/${study}/${study}.pe.genes.raw.htseq2.tsv ${study_dir}/;
  echo "   performing comparison and printing comparison plots..."
  export PATH=${PARASITE_PRODUCTION}/data/conda_envs/wbps-expression/bin:$PATH;
  Rscript ${WORM_CODE}/parasite/scripts/production/rnaseq/brc4/comparison_plots.R ${study_dir} ${study};
  echo "   Done"
  echo ""
  done


