#!/usr/local/bin/perl5.6.1 -w
#
# EMBLDump.pl :  makes EMBL dumps from camace.
# 
#  Last updated on: $Date: 2003-05-19 15:15:59 $
#  Last updated by: $Author: dl1 $


$0 =~ s/^\/.+\///;
system ("touch /wormsrv2/logs/history/$0.`date +%y%m%d`");

use strict;
use lib "/wormsrv2/scripts";
use Wormbase;

# variables

# Need new Transcript dumping code for this part
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

# make clone2name hash
                                                          
my %clone2type;
my ($clone,$type,$processing);

$ENV{'ACEDB'} = "/nfs/disk100/wormpub/DATABASES/current_DB";
my $command = "Table-maker -p \"/wormsrv2/autoace/wquery/clone2type.def\"\nquit\n";

open (TACE, "echo '$command' | $tace | ");
while (<TACE>) {
    chomp;
    s/acedb\> //g;      # only need this is using 4_9i code, bug fixed in 4_9k onward (should be redundant)
    next if ($_ eq "");
    next if (/\/\//);
    s/\"//g;
    (/^(\S+)\s/);
     
    ($clone,$type) = split /\t/;
    $type =~ tr/[A-Z]/[a-z]/;
    $clone2type{$clone} = $type;
    ($clone2type{$clone} = "YAC") if ($type eq "yac");

#    print "Assigning $type to clone $clone\n";
}
close TACE;
                                                          
# cycle through the EMBL dump
                                                          
open (EMBL2, ">/nfs/disk100/wormpub/tmp/EMBLdump.mod") or  die "Can't process new EMBL dump file\n";
open (EMBL,  "<$outfilename.embl") or die "Can't process EMBL dump file\n";
while (<EMBL>) {

    # DE   Caenorhabditis elegans cosmid C05G5
    
    if (/^DE   Caenorhabditis elegans cosmid (\S+)/) {
	if ($clone2type{$1} eq "") {
	    print EMBL2 "DE   Caenorhabditis elegans clone $1\n";
	    }
	elsif ($clone2type{$1} eq "other") {
	    print EMBL2 "DE   Caenorhabditis elegans clone $1\n";
	}
	else {
	    print EMBL2 "DE   Caenorhabditis elegans $clone2type{$1} $1\n";
	}
	    next;
    }

    # species line
    if (/\/organism/) {
	print EMBL2 "FT                   /db_xref=\"taxon:6239\"\n";
	print EMBL2 "$_";
	print EMBL2 "FT                   /strain=\"Bristol N2\"\n";
	print EMBL2 "FT                   /mol_type=\"genomic DNA\"\n";
    }
    
    print EMBL2;

}
close EMBL;
close EMBL2;
                                                          
# copy modified copy back onto output file
system ("mv -f /nfs/disk100/wormpub/tmp/EMBLdump.mod $outfilename.embl");

print "Outfile is $outfilename";
exit;








