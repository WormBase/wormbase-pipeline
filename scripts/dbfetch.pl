#!/usr/local/bin/perl
#
# dbfetch
#
# Usage: dbfetch <database> <name>
#
# Last edited by: $Author: dl1 $
# Last edited on: $Date: 2005-07-11 12:14:50 $


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
