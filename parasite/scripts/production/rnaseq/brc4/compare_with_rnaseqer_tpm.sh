RESULTS_DIR=$1
RNASEQER_FTP="/nfs/ftp/public/databases/arrayexpress/data/atlas/rnaseq"
RFTP_SHORT_STUDIES=$(ls $RNASEQER_FTP -1)
export PATH=${PARASITE_PRODUCTION}/data/conda_envs/wbps-expression/bin:$PATH;

for study_run in $(ls $RESULTS_DIR -1);
  do
  study=$(echo $study_run | cut -d"_" -f1);
  run=$(echo $study_run | cut -d"_" -f2);
  short_run=${run:0:6};
  run_last="${run: -1}"
  run_dir=${RESULTS_DIR}/${study_run}/${run};
  brc4_tpm=$(ls ${run_dir}/*.tpm -1 | grep -v nonunique | head -n1);
  if [ -f "$brc4_tpm" ]; then :; else break; fi;
  if $(echo ${RFTP_SHORT_STUDIES} | grep -q ${short_run}); then :; else break; fi;
  if [ -d "${RNASEQER_FTP}/${short_run}/${run}" ];
    then RNDIR="${RNASEQER_FTP}/${short_run}/${run}";
  elif [ -d "${RNASEQER_FTP}/${short_run}/00${run_last}/${run}" ];
    then RNDIR="${RNASEQER_FTP}/${short_run}/00${run_last}/${run}";
  else break;
  fi;
  if [ -f "${RNDIR}/${run}.pe.genes.tpm.htseq2.irap.tsv" ];
   then irap_tpm="${RNDIR}/${run}.pe.genes.tpm.htseq2.irap.tsv";
  else break;
  fi;
  echo "run: "${run}
  echo "   copying rnaseqer tpm file..."
  cp ${irap_tpm} ${run_dir}/irap.tpm.tsv;
  echo "   performing comparison and printing comparison plots..."
  Rscript ${WORM_CODE}/parasite/scripts/production/rnaseq/brc4/comparison_plots_tpm.R ${brc4_tpm} ${run_dir}/irap.tpm.tsv ${run};
  echo "   Done"
  echo ""
done


