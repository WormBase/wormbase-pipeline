#!/usr/local/bin/perl
#
# usage: repeatmasker2ace.pl directory > outfile

use strict;
my ($score,$query,$qstart,$qend,$rep,$repfam,$try1,$try2,$try3);
my ($repstart,$repend);
my $dir   = shift;
my $type  = shift;
my @files = glob("$dir/*.dna.out");

foreach my $file (@files) {

	open(F,$file);
	my $count = 1;
	while (<F>) {
		/(\d+)\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+\(\d+\)\s+\S\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ and do {
			($score,$query,$qstart,$qend,$rep,$repfam,$try1,$try2,$try3) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);
			if ($try1 =~ /\(/) {
				$repstart = $try2;
				$repend   = $try3; 
			}
			else {
				$repstart = $try1;
				$repend   = $try2; 
			}
		};
		if ($query =~ /\S+/) {
			print "Sequence \"$query\"\n";
			print "Motif_homol \"$rep\" \"$type\" $score $qstart $qend $repstart $repend\n\n";	
		}
		else {
			print STDERR "Problem with $file line $count\n";
		}
		$count++;
	}
	close F;
}
