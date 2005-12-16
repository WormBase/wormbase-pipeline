#!/usr/local/bin/perl
#
# dbfetch
#
# Usage: dbfetch <database> <name>
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2005-12-16 11:18:55 $


my $database = shift;
my $file     = shift;
my $seq      = "";

$/=">";

open (SHOTGUN, "<$database");
while (<SHOTGUN>) {
    chop $_;
    /^(\S+)\s+(\S+)/;
    if (($file eq $1) || ($file eq $2)) {
	$seq = $_;
	last;
    }
}
close SHOTGUN;

$/="\n";

print ">$seq\n";

exit(0);
