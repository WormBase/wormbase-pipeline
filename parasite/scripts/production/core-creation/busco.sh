# Run busco

set -euo pipefail

if [ ! -f "$1" ] ; then echo "Usage: $0 <species/species.fa>. NEMATODES ONLY" ; exit 1 ; fi 
 
fasta=$1
species=$(basename $(dirname $fasta) )
core_db=$($PARASITE_STAGING_MYSQL -e 'show databases' | grep "${species}_core_${PARASITE_VERSION}_${ENSEMBL_VERSION}" | head -n 1)
if [ ! "core_db" ] ; then echo "Could not find core db for species $species " ; exit 1 ; fi

if [ 1 -ne $($PARASITE_STAGING_MYSQL --column-names=FALSE $core_db -e 'select count(*) from meta where meta_value="Nematoda";' ) ] ; then
   >&2 echo "$species does not look like a nematode - which is bad because we don't know what Augustus parameters and BUSCOs to use"
   exit 1
fi

BUSCO_TMP=/nfs/nobackup/ensemblgenomes/wormbase/parasite/busco/WBPS${PARASITE_VERSION}/$species

mkdir -p $BUSCO_TMP

cd $BUSCO_TMP

module load busco

python3 $BUSCO/BUSCO.py -sp caenorhabditis -l $BUSCO/nematoda_odb9 -o $species -i $fasta -c 8 -m genome -f -r 

result=$BUSCO_TMP/run_$species/short_summary_${species}.txt

if [ ! -f "$result" ] ; then >&2 echo "Could not find the result file $result - did BUSCO succeed? " ; exit 1 ; fi

perl -ne 'print "assembly.busco_complete\t$1\nassembly.busco_duplicated\t$2\nassembly.busco_fragmented\t$3\nassembly.busco_missing\t$4\nassembly.busco_number\t$5\n" if /C:([0-9.]+).*D:([0-9.]+).*F:([0-9.]+).*M:([0-9.]+).*n:(\d+)/' $result \
  | while read meta_key meta_value ; do
  ${PARASITE_STAGING_MYSQL}-ensrw $core_db -e "insert into meta (meta_key, meta_value) values (\"$meta_key\", \"$meta_value\");"
done 
