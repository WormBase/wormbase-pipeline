RESULTS_DIR=$1
RNASEQER_FTP="/nfs/ftp/public/databases/arrayexpress/data/atlas/rnaseq/"


for study in $(ls $RESULTS_DIR -1);
  do
  short_study=${study:0:6}
  study_dir=${RESULTS_DIR}/${study};
  spe_cies_bp=$(basename ${RESULTS_DIR});

  echo "study: "${study}
  echo "   summing up metazoa pipeline output..."
  fsuc=${study_dir}/genes.htseq-union.firststrand.counts
  ssuc=${study_dir}/genes.htseq-union.secondstrand.counts
  usuc=${study_dir}/genes.htseq-union.unstranded.counts
  if [[ -f ${fsuc} && -f ${ssuc} ]]; then
    cat ${fsuc} ${ssuc} | awk '{a[$1]+=$2}END{for(i in a) print i,a[i]}' > ${study_dir}/${study}.genes.raw.htseq2.tsv;
  elif [[ -f ${usuc} ]]; then
    cat ${usuc} > ${study_dir}/{study}.genes.raw.htseq2.tsv;
  else
    echo "ERROR: No multiple (genes.htseq-union.firststrand.counts and genes.htseq-union.secondstrand.counts)
          or a single genes.htseq-union.unstranded.counts
          files found in ${study_dir}." 1>&2
    exit 1
  fi;



  cat ${study_dir}/genes.htseq-union.firststrand.counts ${study_dir}/genes.htseq-union.secondstrand.counts | awk '{a[$1]+=$2}END{for(i in a) print i,a[i]}' > ${study_dir}/${study}.counts;
  #  echo "   copying rnaseqer file..."
  cp ${RNASEQER_FTP}/${short_study}/${study}/${study}.pe.genes.raw.htseq2.tsv ${study_dir}/;
  echo "   performing comparison and printing comparison plots..."
  export PATH=${PARASITE_PRODUCTION}/data/conda_envs/wbps-expression/bin:$PATH;
  Rscript ${WORM_CODE}/parasite/scripts/production/rnaseq/brc4/comparison_plots.R ${study_dir} ${study};
  echo "   Done"
  echo ""
  done


