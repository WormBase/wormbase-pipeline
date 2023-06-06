splitfa(){
i=1;
cat $1 | while read line ; do
  if [ ${line:0:1} == ">" ] ; then
    header="$line"
    echo "$header" >> seq"${i}".fasta
    ((i++))
  else
    seq="$line"
    echo "$seq" >> seq"${i}".fasta
  fi
done
}