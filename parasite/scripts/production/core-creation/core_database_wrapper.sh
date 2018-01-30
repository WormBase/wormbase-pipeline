#!/bin/bash 
set -e
DIR=`dirname $0`

#this wrapper creates one or several core dbs from a YAML config file then performs some basic healthchecks on them

#Creating core database with worm_lite.pl

if ! [ -f $1 ] ; then echo "Usage: $0 <species conf>" ; exit 1 ; fi

out_file=$( basename "$1" | sed 's/\..*$/.wormlite.out/' ) 
SPECIES_CONF=$1.auto

printf "Expanding the config \n"
perl -MCoreCreation::Conf -e "print CoreCreation::Conf->new(\"$1\")-> dump()" > $SPECIES_CONF

printf "Launching worm_lite.pl \n"
printf "The output of this script will be written in $out_file \n\n"

if perl -w $WORM_CODE/scripts/ENSEMBL/scripts/worm_lite.pl -yfile $SPECIES_CONF -allspecies -setup -load_genes -load_dna &> $out_file; then
    printf "worm_lite.pl exited successfully.\n\n"
else
    printf "worm_lite has errored. Please check the output file\n"
    printf "Exiting now.\n\n"
    exit 1
fi

#performing basic healthchecks
git=`grep -o 'cvsdir: .*' $SPECIES_CONF | awk '{print $2}'`

# TODO this used to be for a list, now it's for one species. Does no longer work for multiple species.
species_list=`grep '^[a-z0-9]*_[a-z0-9]*_[a-z0-9]*:' $SPECIES_CONF`
fasta_list=`grep 'fasta: .*' $SPECIES_CONF | awk '{print $2}'`
gff3_list=`grep 'gff3: .*' $SPECIES_CONF | awk '{print $2}'`
coredb=$(perl -MCoreCreation::Conf -e "print CoreCreation::Conf->new(\"/nfs/production/panda/ensemblgenomes/wormbase/parasite/data/release10/ascaris_suum_prjna62057/ascaris_suum_prjna62057.conf\")->{core_database}->{dbname} ;")
host=$(perl -MCoreCreation::Conf -e "print CoreCreation::Conf->new(\"/nfs/production/panda/ensemblgenomes/wormbase/parasite/data/release10/ascaris_suum_prjna62057/ascaris_suum_prjna62057.conf\")->{core_database}->{host};")
port=$(perl -MCoreCreation::Conf -e "print CoreCreation::Conf->new(\"/nfs/production/panda/ensemblgenomes/wormbase/parasite/data/release10/ascaris_suum_prjna62057/ascaris_suum_prjna62057.conf\")->{core_database}->{port};")

END=`printf "%s\n" $species_list | wc -l`

for i in $(seq 1 $END);
do
species=`sed -n ${i}p <<< "$species_list" | sed s'/.$//'`
fasta=`sed -n ${i}p <<< "$fasta_list"`
gff3=`sed -n ${i}p <<< "$gff3_list"`

health_out="/nfs/panda/ensemblgenomes/wormbase/parasite/core-creation/$coredb.healthcheck.out"
printf "Submitting healthcheck job for $species \n"
printf "The output of this job will be written in $health_out \n\n"
set -x
$DIR/healthcheck.sh -d $coredb -g $fasta -a $gff3 -u ensro -p $port -s $host > $health_out
set +x
done

printf "All healthchecks done!\n"


