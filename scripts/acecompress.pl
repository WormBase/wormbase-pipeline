#!/usr/local/bin/perl5.6.1 -w
#
# acecompress.pl
#
# dl
#
# Usage : acecompress.pl [-options]
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2003-01-28 17:14:17 $

use strict;
use Getopt::Long;

my $homol;         # homol data   BLAT files
my $feature;       # feature data Good_intron files    
my %objects;       # hash for number of homol lines per object
my $ace_object;  
my %acedata;


GetOptions (
	    "homol"       => \$homol,
	    "feature"     => \$feature
	    );

my $file = shift;

open (FILE, "<$file");
while (<FILE>) {

    # BLAT homology data
    if ($homol) {
	if (/Homol_data : \"(\S+)\"/) {
	    $ace_object = $1;
	    $objects{$ace_object}++;
	    next;
	}
	if (/DNA_homol/) {
	    push (@{$acedata{$ace_object}{DNA_homol}},$_);
	}
    }
    # Confirmed_intron feature data
    elsif ($feature) {
	if (/Feature_data : \"(\S+)\"/) {
	    $ace_object = $1;
	    $objects{$ace_object}++;
	    next;
	}
	if (/Confirmed_intron/) {
	    push (@{$acedata{$ace_object}{Confirmed_intron}},$_);
	}
    }
}
close FILE;


# print output

foreach my $obj (keys %objects) {
    print "\n// $obj\n\n";

    # homol (BLAT homology)
    if ($homol) {
	print "Homol_data : \"$obj\"\n";
	foreach my $line (@{$acedata{$obj}{DNA_homol}}) {
	    print "$line";
	}
    }
    # feature (confirmed_intron)
    elsif ($feature) {
	print "Feature_data : \"$obj\"\n";
	foreach my $line (@{$acedata{$obj}{Confirmed_intron}}) {
	    print "$line";
	}
    }  
}

exit(0);
