#!/usr/local/bin/perl

# converts ensembl.pep file with human peptide entries into ace file

use strict;

while (<>) {
	if (/\>(ENS\S+)/) {
		print "\nProtein : \"ENSEMBL:$1\"\n";
		print "Peptide \"ENSEMBL:$1\"\n";
		print "Species \"Homo sapiens\"\n";
		print "Database \"Ensembl\" \"$1\" \"$1\"\n\n";
		print "Peptide : \"ENSEMBL:$1\"\n"; 	
	}
	else { 
		print $_;
	}
}

__END__
Protein : "ENSEMBL:ENSP00000266127"
Peptide "ENSEMBL:ENSP00000266127"
Species "Homo sapiens"
Database "Ensembl" "ENSP00000266127" "ENSP00000266127"

Peptide : "ENSEMBL:ENSP00000266127"
AFYPPRLFAPPPPWFFPPISPP
