#!/usr/local/bin/perl5.6.1 -w
#
# EMBLDump.pl :  makes EMBL dumps from camace.
# 
#  Last updated on: $Date: 2003-04-14 07:39:25 $
#  Last updated by: $Author: krb $


$0 =~ s/^\/.+\///;
system ("touch /wormsrv2/logs/history/$0.`date +%y%m%d`");

use strict;
use lib "/wormsrv2/scripts";
use Wormbase;

# variables

# Need new Transcript dumping code for this part
my $exec        = &giface;
#my $exec = "/nfs/team71/acedb/edgrif/TEST/KEITH/giface ";
my $dbdir       = "/wormsrv2/camace";
my $tace        = &tace;
my $giface      = "$exec $dbdir";
my $outfilename = "/nfs/disk100/wormpub/tmp/EMBLdump.$$";
our $rundate    = `date +%y%m%d`; chomp $rundate;

# Dump the EMBL file from camace

my $query = "Query Find Genome_Sequence From_Laboratory = HX AND Finished AND DNA\ngif EMBL $outfilename\n";
print "Query = $query\n";
print "Exec = $giface\n";
open (READ, "echo '$query' | $giface |") or die ("Could not open $giface\n"); 
while (<READ>) {
 if ($_ =~ /\/\//) {next};
 if ($_ =~ /acedb/) {next};
}                   
close READ;


print "Outfile is $outfilename";
exit;








