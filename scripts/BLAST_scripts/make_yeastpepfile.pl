#!/usr/local/bin/perl5.6.1 -w

# krb 020829

# converts yeastX.pep file with yeast peptide entries into ace file
# Puts SGDID as Accssion in Database field.

use strict;

while (<>) {
	if (/\>ORFP:([\w-]+)\s+([\w-',\(\)]+)\s+SGDID:(\w+)/) {
		print "\nProtein : \"SGD:$1\"\n";
		print "Peptide \"SGD:$1\"\n";
		print "Species \"Saccharomyces cerevisiae\"\n";
		print "DB_remark \"SGD gene name is $2\"\n";
		print "Database \"SGD\" \"$1\" \"$3\"\n\n";
		print "Peptide : \"SGD:$1\"\n"; 	
	}
	else { 
		print $_;
	}
}

__END__

