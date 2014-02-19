#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  fix_treefam_from_ace.pl
#
#        USAGE:  ./fix_treefam_from_ace.pl INFILE
#
#  DESCRIPTION: that is a crude script to fix thr treefam dump file
#               needs a working ACeDB server on wormsrv2 
#
#      VERSION:  1.0
#      CREATED:  01/29/08 16:58:10 GMT
#===============================================================================
use Ace;
use strict;

# hardcoded settings
#my $db=Ace->connect(-host => 'wormsrv2',-port => 23100);
my $db=Ace->connect(-path => shift);

while (<>){
        if (/\>/){
	 	my @a=split;
		print STDERR join(@a,"\n") unless $a[2];
		my $cds=$db->fetch(CDS => $a[2]);
		my $gene=$cds->Gene;
		my $public_name=$gene->Public_name;
		die("can't get cds for $cds / gene for $gene / public_name for $public_name") unless ($cds && $gene && $public_name);
		s/>\S+/>$gene/;
		s/\S+_/${public_name}_/;
		print;
	}else {print}
}


