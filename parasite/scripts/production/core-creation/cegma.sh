# Run cegma

set -e

if [ ! -f "$1" ] ; then echo "Usage: $0 <species/species.fa>" ; exit 1 ; fi 
 
fasta=$1
species=$(basename $(dirname $fasta) )
core_db=$($PARASITE_STAGING_MYSQL -e 'show databases' | grep "${species}_core_${PARASITE_VERSION}_${ENSEMBL_VERSION}" | head -n 1)
if [ ! "core_db" ] ; then echo "Could not find core db for species $species " ; exit 1 ; fi

tmp=$(dirname $fasta)/.tmp
mkdir -pv $tmp
cd $tmp
echo "Running CEGMA"
/nfs/panda/ensemblgenomes/wormbase/software/packages/cegma/runcegma_kog.sh -g $fasta

cegma_complete=$(perl -ne 'print $1 if /Complete\s*[0-9.]+\s*([0-9.]+).*/' ./output.completeness_report)
cegma_partial=$(perl -ne 'print $1 if /Partial\s*[0-9.]+\s*([0-9.]+).*/' ./output.completeness_report)
 
if ! [ "$cegma_complete" ] || ! [ "$cegma_partial" ] ; then
    echo "Could not parse: ./output.completeness_report"
    exit 1
fi
 
${PARASITE_STAGING_MYSQL}-ensrw $core_db -e "insert into meta (meta_key, meta_value) values (\"assembly.cegma_complete\", \"$cegma_complete\");"
${PARASITE_STAGING_MYSQL}-ensrw $core_db -e "insert into meta (meta_key, meta_value) values (\"assembly.cegma_partial\", \"$cegma_partial\");"
echo "Complete!"
rm -v $tmp
