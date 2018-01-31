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

$DIR/healthcheck.pl $1
