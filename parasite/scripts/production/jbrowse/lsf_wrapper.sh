
if [ $# -ne 3 ] ; then
  >&2 echo "Usage: $0 <ftp dir> <jbrowse_path> <out_dir> "
  exit 1
fi
log=$PARASITE_SCRATCH/log/jbrowse/$(date -u +"%Y-%m-%dT%H:%M" )
mkdir -p $log

#Conveniently Ensembl has all the right modules installed for jbrowse
# You could install them yourself - see Makefile.pl in the repo
PERL5LIB=/nfs/public/rw/ensembl/libs/perl/share/perl5/:$PERL5LIB
export PERL5LIB

#species=($(ls -d $1/species/*))

for s in /nfs/nobackup/ensemblgenomes/wormbase/parasite/production/log/jbrowse/2018-03-05T16:52/*out.failed ; do
  s="$(basename $s)"
  s="$(cut -f 1 -d . <<< $s)"
  bsub -M4000 -o $log/${s}.out -e $log/${s}.err  "perl $(dirname $0)/create_species_config.pl --ftp_dir $1 --jbrowse_path $2 --out_dir $3 --species $s"
done

