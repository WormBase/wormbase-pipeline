#!/usr/local/bin/perl5.6.0

$|=1;
use lib "/wormsrv2/scripts/";   
use strict;
use Wormbase;

my $maintainers = "All";
my $log         = "/tmp/whoswho.log";
my $line_count  = 0;

# move to dl1 cvs checkout directory
chdir ("/nfs/griffin2/dl1/wormbase/wormbase/scripts");

open (LOG, ">$log");
print LOG "Current checkouts from the CVS repository\n";
print LOG "==============================================================\n\n";

open (CVS, "/usr/local/bin/cvs editors |");
while (<CVS>) {
    print LOG;
    $line_count++;
}
close CVS;

print LOG "\n============================================================\n\n";

unless ($line_count == 0) {
    &mail_maintainer("CVS repository report:",$maintainers,$log);
}

exit(0);
