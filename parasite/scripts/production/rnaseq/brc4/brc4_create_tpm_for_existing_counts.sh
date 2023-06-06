fl_file=$1
sum_htseq_file=$2

rpk_htseq_file="${sum_htseq_file}.rpk"
tpm_htseq_file="${sum_htseq_file}.tpm2"

sum_htseq_file_lines=$(cat $sum_htseq_file | grep -v "^__" | wc -l)

echo "joining"
join <(sort -k1,1 $sum_htseq_file) <(sort -k1,1 $fl_file) | awk '{print $1"\t"$2/($3/1000)}' > ${rpk_htseq_file}

rpk_htseq_file_lines=$(cat $rpk_htseq_file | wc -l)

echo "qc1"
if [ "${sum_htseq_file_lines}" -ne "${rpk_htseq_file_lines}" ];
  then echo "Problems with line count ${sum_htseq_file}. Exiting";
  exit 1;
fi

echo "mf"
million_factor=$(awk '{ sum += $2 } END { print sum/1000000 }' $rpk_htseq_file)

echo "tpm"
awk -v million_factor="$million_factor" '{print $1"\t"$2/million_factor}' $rpk_htseq_file > $tpm_htseq_file

tpm_htseq_file_lines=$(cat $tpm_htseq_file | wc -l)

echo "qc2"
if [ "${sum_htseq_file_lines}" -ne "${tpm_htseq_file_lines}" ];
  then echo "Problems with line count ${sum_htseq_file}. Exiting";
  exit 1;
fi

rm $rpk_htseq_file

echo "Done: Created $tpm_htseq_file"
exit 0;