#!/usr/local/bin/perl5.8.0 -w
#
# EMBLDump.pl :  makes EMBL dumps from camace.
# 
#  Last updated on: $Date: 2005-03-22 11:24:09 $
#  Last updated by: $Author: dl1 $

use strict;
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use File::Copy;

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

#my $giface      = &giface;

my $giface =  glob("~edgrif/TEST/DAN/giface");
my $dbdir       = "$basedir/camace";
my $tace        = &tace;
my $outfilename = "/nfs/disk100/wormpub/tmp/EMBLdump.$$";
my $current_DB  = "/nfs/disk100/wormpub/DATABASES/current_DB";
my $mod_file    = "/nfs/disk100/wormpub/tmp/EMBLdump.mod";

if( $test ) {
  $giface    = glob("~edgrif/TEST/DAN/giface");
  $outfilename   = "/nfs/disk100/wormpub/test/EMBLdump.$$";
  $mod_file  = "/nfs/disk100/wormpub/test/EMBLdump.mod";
}

#########################
# make history log item
#########################
$0 =~ s/^\/.+\///;
system ("touch $basedir/logs/history/$0.`date +%y%m%d`");

#############################################
# Use giface to dump EMBL files from camace #
#############################################


my $query = "nosave\nquery find CDS where Method = \"Genefinder\"\nkill\nquery find CDS where Method = \"twinscan\"\nkill\nquery find Genome_sequence From_laboratory = HX AND Finished AND DNA\ngif EMBL $outfilename\nquit\nn\n";
$query = "query find Genome_sequence AH6\ngif EMBL $outfilename\n" if $test;

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
my ($clone,$type);

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

#############################################
# make CDS2CGC hash from info in current_DB #
#############################################

my %cds2cgc;
my %cds2gene;
my ($cds,$cgc,$wbgene);

$command = "Table-maker -p \"$basedir/autoace/wquery/cgc2cds.def\"\nquit\n";

open (TACE, "echo '$command' | $tace $current_DB | ");
while (<TACE>) {
    chomp;
    s/acedb\> //g;        # only need this if using 4_9i code
    next if ($_ eq "");
    next if (/\/\//);
    s/\"//g;

    if (/^(\S+)\s+(\S+)\s+(\S+)/) {
	($wbgene,$cgc,$cds) = ($1,$2,$3);
	$cds2cgc{$cds} = $cgc;
	$cds2gene{$cds} = $wbgene;
    }
}
close TACE;

######################################################################                     
# cycle through the EMBL dump file, replacing info where appropriate #
######################################################################

open (OUT, ">$mod_file") or  die "Can't process new EMBL dump file\n";
open (EMBL, "<$outfilename.embl") or die "Can't process EMBL dump file\n";

my $id = "";
our $reference_remove = 0;
our $author_change;

while (<EMBL>) {

  # print ID line and next XX. store id
  if(/^ID\s+CE(\S+)/){
    $id = $1;
    print OUT "${_}XX\n";
    next;
  }

  # print ID line and next XX
  if( /^AC/ ) {
    print OUT "AC   ",substr($clone2sv{$id},0,-2),";\nXX\n";
    next;
  }

  # DE   Caenorhabditis elegans cosmid C05G5    
  if (/^DE   Caenorhabditis elegans cosmid (\S+)/) {
    my $clone = $1;

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

  # print OC line and next XX, tag for WormBase inclusion
  if (/^RP\s+(\S+)/) {
      $author_change = $1;
      print OUT "$_";
      print OUT "RX   MEDLINE; 99069613.\n";
      print OUT "RX   PUBMED; 9851916.\n";
      print OUT "RG   WormBase Consortium\n";
      print OUT "RA   ;\n";
      print OUT "RT   \"Genome sequence of the nematode C. elegans: a platform for investigating\n";
      print OUT "RT   biology\";\n";
      print OUT "RL   Science 282(5396):2012-2018(1998).\n";
      print OUT "XX\n";
      print OUT "RN   [2]\n";
      print OUT "$_";
      next;
  }

  if (/^RN   \[2\]/) {
      $reference_remove = 6;
      next;
  }

  if ($reference_remove > 0) {
      $reference_remove--;
      next;
  }


  if (/^RL   E-mail: jes/) {
      next;
  }


  # locus_tag name.....
  if (/\/gene=\"(\S+)\"/) {
      $cds = $1;
      if ($cds2cgc{$cds}) {
	  print OUT "FT                   /gene=\"" . $cds2cgc{$cds}  ."\"\n";
	  print OUT "FT                   /locus_tag=\"$cds\"\n";
	  next;
      }
      else {
	  print OUT "FT                   /gene=\"$cds\"\n";
	  print OUT "FT                   /locus_tag=\"$cds\"\n";
	  next;
      }
  }	  

  next if (/^CC   For a graphical/);
  next if (/^CC   see:-/);

  if (/^CC   name=/) {
      print OUT "CC   For a graphical representation of this sequence and its analysis\n";
      print OUT "CC   see:- http://www.wormbase.org/perl/ace/elegans/seq/sequence?\n";
      print OUT "CC   name=$clone;class=Sequence\n";
      $reference_remove = 1;
      next;
  }

  # standard_name......

  # don't print out first few lines until they have been converted 
  # can only do this when it gets to DE line
  print OUT if ($id eq "");
  
}
close EMBL;
close OUT;
                                                          
# copy modified copy back onto output file

my $status = move("$mod_file","$outfilename.embl");
print "ERROR: Couldn't move file: $!\n" if ($status == 0);

print "\nOutfile is $outfilename\n\n";

exit(0);








