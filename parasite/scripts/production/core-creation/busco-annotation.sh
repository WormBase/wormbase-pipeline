# Run busco in "proteins" mode
set -e

# USAGE
if [ ! "$1" ] ; then echo "Usage: $0 core_db [roundworm/flatworm]" ; exit 1 ; fi

#Set up DB credentials
DBHOSTNAME=$($PARASITE_STAGING_MYSQL details script | grep -Po "^--host \K(\S+)")
DBUSER=$($PARASITE_STAGING_MYSQL details script | grep -Po "\s+--user \K(\S+)")
DBPORT=$($PARASITE_STAGING_MYSQL details script | grep -Po "\s+--port \K(\S+)")

#Input
core_db=$1
phylum=${2}

species=$(echo $core_db | cut -d'_' -f1,2,3)

#Check if core_db in current db
CORE_DBS=$($PARASITE_STAGING_MYSQL -Ne "SHOW DATABASES LIKE \"%_core_${PARASITE_VERSION}_${ENSEMBL_VERSION}%\"")
if [[ ! $CORE_DBS =~ (^|[[:space:]])$core_db($|[[:space:]]) ]]; then echo "core_db not in $PARASITE_STAGING_MYSQL. Check and Run again."; exit 1; fi

#Check busco software
module load busco
#if [ ! -e "$BUSCO4_CONTAINER" ] ; then echo "Importing BUSCO singularity image didn't work!" ; exit 1 ; fi
#if [ ! -e "$BUSCO5_CONTAINER" ] ; then echo "Importing BUSCO singularity image didn't work!" ; exit 1 ; fi
if [ ! -e "$BUSCO3_CONTAINER" ] ; then echo "Importing BUSCO singularity image didn't work!" ; exit 1 ; fi
if [ ! -d "$BUSCO_LINEAGES" ] ; then echo "Importing BUSCO lineages didn't work!" ; exit 1 ; fi

#Check if phylum is valid
if [ "$phylum" == roundworm ] || [ 1 -eq "$($PARASITE_STAGING_MYSQL --column-names=FALSE $core_db -e 'select count(*) from meta where meta_value="Nematoda";' )" ] ; then
  species_parameter_for_augustus=caenorhabditis
  busco_library=$BUSCO_LINEAGES/nematoda_odb9
elif [ "$phylum" == flatworm ] || [ 1 -eq "$($PARASITE_STAGING_MYSQL --column-names=FALSE $core_db -e 'select count(*) from meta where meta_value="Platyhelminthes";')" ] ; then
  species_parameter_for_augustus=schistosoma
  busco_library=$BUSCO_LINEAGES/metazoa_odb9
else
   >&2 echo "$species does not look like a nematode or a platyhelminth- which is bad because we don't know what Augustus parameters and BUSCOs to use"
   exit 1
fi


#Dump proteins Run
BUSCO_TMP=$PARASITE_SCRATCH/busco-annotation/WBPS${PARASITE_VERSION}/$species
mkdir -p $BUSCO_TMP
cd $BUSCO_TMP


run_log_dpl=$BUSCO_TMP/run-dump-proteins.$(date "+%Y-%m-%d").out
perl $WORM_CODE/scripts/ENSEMBL/scripts/dump_proteins.pl -canonical_only --host=$DBHOSTNAME --port=$DBPORT --user=$DBUSER --dbname=$core_db --outfile=$species.prot.fa \
  | tee $run_log_dpl

## Check if we got the expected number of protein sequences we thought we would get
EXPECTED_PROTEIN_SEQ_NUMBER=$($PARASITE_STAGING_MYSQL $core_db --column-names=FALSE -e 'select count(*) from transcript where biotype="protein_coding" AND canonical_translation_id IS NOT NULL;')
DUMPED_PROTEIN_SEQ_NUMBER=$(grep -c ">" $species.prot.fa)
if [[ ! "$EXPECTED_PROTEIN_SEQ_NUMBER" -eq "$DUMPED_PROTEIN_SEQ_NUMBER" ]]; then
  printf "\n\nERROR: Did not get the expected amount of sequences dumped:\nExpected protein sequences: %s\nDumped protein sequences: %s.\n" "$EXPECTED_PROTEIN_SEQ_NUMBER" "$DUMPED_PROTEIN_SEQ_NUMBER" \
  | tee -a $run_log_dpl; fi
printf "Expected protein sequences: %s\nDumped protein sequences: %s.\n" "$EXPECTED_PROTEIN_SEQ_NUMBER" "$DUMPED_PROTEIN_SEQ_NUMBER" >> $run_log_dpl



# BUSCO Protein Run
run_log=$BUSCO_TMP/run-busco.$(date "+%Y-%m-%d").out

singularity exec $BUSCO3_CONTAINER run_BUSCO.py -sp $species_parameter_for_augustus -l $busco_library -o $species -i $species.prot.fa -c 8 -m proteins -f \
  | tee $run_log


# Output Handling
## BUSCO doesn't reliably use status codes
## https://gitlab.com/ezlab/busco/issues/84
if grep '^CRITICAL\|^ERROR' $run_log ; then
  echo "Run log said a worrying thing. Bailing out!"
  exit 1
fi

result=${BUSCO_TMP}/run_${species}/short_summary_${species}.txt

if [ ! -f "$result" ] ; then >&2 echo "Could not find the result file $result - did BUSCO succeed? " ; exit 1 ; fi
if [ ! "$core_db" ] ; then  >&2 echo "No core db - go read $result. Finishing " ; exit; fi

#Update the DB
${PARASITE_STAGING_MYSQL}-w $core_db -e 'delete from meta where meta_key like "annotation.busco%"'
${PARASITE_STAGING_MYSQL}-w $core_db -e 'delete from meta where meta_key like "annotation.busco3%"'
echo "Parsing the result file: $result"
perl -ne 'print "annotation.busco3_complete\t$1\nannotation.busco3_duplicated\t$2\nannotation.busco3_fragmented\t$3\nannotation.busco3_missing\t$4\nannotation.busco3_number\t$5\n" if /C:(-?[0-9.]+).*D:(-?[0-9.]+).*F:(-?[0-9.]+).*M:(-?[0-9.]+).*n:(\d+)/' $result \
  | while read meta_key meta_value ; do
  ${PARASITE_STAGING_MYSQL}-w $core_db -e "insert into meta (meta_key, meta_value) values (\"$meta_key\", \"$meta_value\");"
done

#Delete the protein fasta
rm ${BUSCO_TMP}/${species}.prot.fa

echo "test done"
#echo "Parsed the results and inserted into the meta table. DONE"

sleep 2
exit