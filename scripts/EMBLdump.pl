#!/usr/local/bin/perl5.6.1 -w
#
# EMBLDump.pl :  makes EMBL dumps from camace.
# 
#  Last updated on: $Date: 2002-12-13 12:52:04 $
#  Last updated by: $Author: dl1 $


$0 =~ s/^\/.+\///;
system ("touch /wormsrv2/logs/history/$0.`date +%y%m%d`");

use strict;
use lib "/wormsrv2/scripts";
use Wormbase;

# variables
my $exec        = &giface;
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

######################################
# do some  modifications of the file #
######################################

# make CDS2locus hash
my %CDS2locus;
my ($locus,$cds,$processing);    
    
$ENV{'ACEDB'} = "/wormsrv2/current_DB";
my $command = "Table-maker -p \"/wormsrv2/autoace/wquery/CDS2locus.def\"\nqui
t\n";


open (TACE, "echo '$command' | $tace | ");
while (<TACE>) {
    chomp;
    s/acedb\> //g;      # only need this is using 4_9i code, bug fixed in 4_9k onward (should be redundant)
    next if ($_ eq "");
    next if (/\/\//);
    s/\"//g;
    (/^(\S+)\s/);

    ($locus,$cds) = split /\t/;
    $CDS2locus{$cds} = $locus;
    print "Assigning $locus to $cds\n";
}
close TACE;

# cycle through the EMBL dump
    
open (EMBL2, ">/nfs/disk100/wormpub/tmp/EMBLdump.mod") or  die "Can't process new EMBL dump file\n";
open (EMBL,  "<$outfilename.embl") or die "Can't process EMBL dump file\n";
while (<EMBL>) {
    
    # gene line
    if (/\/gene=\"(\S+)\"/) {
	$processing = $1;
	# is this a defined locus?
	if ($CDS2locus{$processing} ne "") {
	    print EMBL2 "FT                   /gene=\"$CDS2locus{$processing}\"\n";
	    print EMBL2 "FT                   /standard_name=\"$processing\"\n";
	    print EMBL2 "FT                   /product=\"C. elegans $CDS2locus{$processing} protein\n";
	    print EMBL2 "FT                   (corresponding sequence $processing)\"\n";}
	else {
	    print EMBL2;
	    print EMBL2 "FT                   /standard_name=\"$processing\"\n";
	    print EMBL2 "FT                   /product=\"Hypothetical protein $processing\"\n";
	}
    }
    else {
	print EMBL2;
    }
}
close EMBL;
close EMBL2;

# copy modified copy back onto output file
system ("mv -f /nfs/disk100/wormpub/tmp/EMBLdump.mod $outfilename.embl");

print "Outfile is $outfilename";
exit;








