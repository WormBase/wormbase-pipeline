#!/usr/local/bin/perl5.8.0 -w
#
# EMBLDump.pl :  makes EMBL dumps from camace.
# 
#  Last updated on: $Date: 2003-12-03 10:06:15 $
#  Last updated by: $Author: krb $

use strict;
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use File::Copy qw(mv);



##############################
# command-line options       #
##############################
my $test;   # use test environment in ~wormpub/TEST_BUILD/

GetOptions ("test"         => \$test);

my $basedir     = "/wormsrv2";
$basedir        = glob("~wormpub")."/TEST_BUILD" if ($test); 



###############################
# misc. variables             #
###############################

my $giface      = &giface;
my $dbdir       = "$basedir/camace";
my $tace        = &tace;
my $outfilename = "/nfs/disk100/wormpub/tmp/EMBLdump.$$";
my $current_DB = "/nfs/disk100/wormpub/DATABASES/current_DB";




#########################
# make history log item
#########################
$0 =~ s/^\/.+\///;
system ("touch $basedir/logs/history/$0.`date +%y%m%d`");



#############################################
# Use giface to dump EMBL files from camace #
#############################################

my $query = "Query Find Genome_Sequence From_laboratory = HX AND Finished AND DNA\ngif EMBL $outfilename\n";

open(READ, "echo '$query' | $giface $dbdir |") or die ("Could not open $giface $dbdir\n"); 
while (<READ>) {
 next if ($_ =~ /\/\//);
 next if ($_ =~ /acedb/);
}                   
close(READ);


########################################################
# make clone2sv hash from info in camace
########################################################

# this is needed to fix missing sequence version and accession info in 
# dumped EMBL files

my %clone2sv;

my $command = "Table-maker -p \"$basedir/autoace/wquery/clone2sv.def\"\nquit\n";

open (TACE, "echo '$command' | $tace $dbdir | ");
while (<TACE>) {
  chomp;
  next if ($_ eq "");
  next if (/\/\//);
  s/acedb\> //g;      # only need this is using 4_9i code, bug fixed in 4_9k onward (should be redundant)
  s/\"//g;

  if($_ =~ m/^(\S+)\s/){
    my ($clone,$sv) = split /\t/;
    $clone2sv{$clone} = $sv;
  }
}
close TACE;


 
###############################################
# make clone2name hash from info in current_DB
###############################################
my %clone2type;
my ($clone,$type,$processing);

$command = "Table-maker -p \"$basedir/autoace/wquery/clone2type.def\"\nquit\n";

open (TACE, "echo '$command' | $tace $current_DB | ");
while (<TACE>) {
  chomp;
  s/acedb\> //g;      # only need this is using 4_9i code, bug fixed in 4_9k onward (should be redundant)
  next if ($_ eq "");
  next if (/\/\//);
  s/\"//g;

  if ($_ =~ m/^(\S+)\s+(\S+)/){    
    ($clone,$type) = split /\t/;
    $type =~ tr/[A-Z]/[a-z]/;
    $clone2type{$clone} = $type;
    ($clone2type{$clone} = "YAC") if ($type eq "yac");
  }
  elsif ($_ =~ m/^(\S+)\s/){    
    print "WARNING: Missing type information?  $_\n";
  }
}
close TACE;


######################################################################                     
# cycle through the EMBL dump file, replacing info where appropriate #
######################################################################

open (OUT, ">/nfs/disk100/wormpub/tmp/EMBLdump.mod") or  die "Can't process new EMBL dump file\n";
open (EMBL, "<$outfilename.embl") or die "Can't process EMBL dump file\n";

my $id = "";

while (<EMBL>) {

  # store ID line to reformat later
  if(/^ID   NDB_ID/){
    $id = $_;
    chomp($id);
  }

  # DE   Caenorhabditis elegans cosmid C05G5    
  if (/^DE   Caenorhabditis elegans cosmid (\S+)/) {
    my $clone = $1;

    # grab accession using sequence version
    my $accession = $clone2sv{$clone};
    substr($accession,-2) = "";

    # correct ID line
    $id =~ s/NDB_ID/CE${clone}/;

    # now print out other lines which have been skipped above
    print OUT "$id\n";
    print OUT "XX\n";
    print OUT "AC   $accession;\n";
    print OUT "XX\n";
    print OUT "SV   $clone2sv{$clone};\n";
    print OUT "XX\n";

    # can now reset $id
    $id = "";

    if (!defined($clone2type{$clone})){
      print OUT "DE   Caenorhabditis elegans clone $clone\n";
      print "WARNING: no clone type for $_" 
    }
    elsif ($clone2type{$clone} eq "other") {
      print OUT "DE   Caenorhabditis elegans clone $clone\n";
    }
    else {
      print OUT "DE   Caenorhabditis elegans $clone2type{$clone} $clone\n";
    }
    next;
  }
  
  # species line
  if (/\/organism/) {
    print OUT "FT                   /db_xref=\"taxon:6239\"\n";
    print OUT "$_";
    print OUT "FT                   /strain=\"Bristol N2\"\n";
    print OUT "FT                   /mol_type=\"genomic DNA\"\n";
    next;
  }

  # don't print out first few lines until they have been converted 
  # can only do this when it gets to DE line
  print OUT if ($id eq "");
  
}
close EMBL;
close OUT;
                                                          
# copy modified copy back onto output file
mv("/nfs/disk100/wormpub/tmp/EMBLdump.mod","$outfilename.embl");

print "\nOutfile is $outfilename\n\n";

exit(0);








