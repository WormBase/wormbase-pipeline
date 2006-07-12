#!/usr/local/bin/perl5.8.0 -w
#
# EMBLdump.pl :  makes modified EMBL dumps from camace.
# 
#  Last updated on: $Date: 2006-07-12 18:44:50 $
#  Last updated by: $Author: pad $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use File::Copy;
use Storable;

##############################
# command-line options       #
##############################

my $test;
my $single;
my $debug;
my $store;
my $wormbase;
my $quicktest;
my $version;

GetOptions (
	    "test"         => \$test,
	    "debug=s"      => \$debug,
	    "store:s"      => \$store,
	    "quicktest"    => \$quicktest,
	    "single=s"     => \$single,
	    "version=s"    => \$version,
	    );


if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
                            -test    => $test,
                           );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

# who will receive log file?
my $maintainers = "All";

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}


###############################
# misc. variables             #
###############################

my $basedir        = $wormbase->wormpub;
my $giface         = $wormbase->giface;
my $dbdir          = $wormbase->database('camace');
my $tace           = $wormbase->tace;
my $outfilename    = "$basedir/tmp/EMBLdump.$$";
#my $current_DB     = $wormbase->database('current');
my $mod_file       = "$basedir/tmp/EMBLdump.mod";
#my (%clone2accession, %clone2type, %cds2cgc,);
my $dir;
$log->write_to("You are embl dumping from $dbdir\n\n");


#############################################
# Use giface to dump EMBL files from camace #
#############################################

my $command;
if (!$single) {
$command  = "nosave\n"; # Don't really want to do this
$command .= "query find CDS where Method = \"Genefinder\"\nkill\n";# remove Genefinder predictions
$command .= "query find CDS where Method = \"twinscan\"\nkill\n";# remove twinscan predictions
$command .= "query find Genome_sequence From_laboratory = HX AND Finished AND DNA\ngif EMBL $outfilename\n";# find Genome_sequences and EMBL dump
$command .= "quit\nn\n";# say you don't want to save and exit
}

# quicktest mode only works on B0250
if ($quicktest) {
    $command    = "query find Genome_sequence B0250\ngif EMBL $outfilename\nquit\n";
}

if ($single) {
  $command  = "nosave\n"; # Don't really want to do this
  $command .= "query find CDS where Method = \"Genefinder\"\nkill\n";# remove Genefinder predictions
  $command .= "query find CDS where Method = \"twinscan\"\nkill\n";# remove twinscan predictions
  $command .= "query find Genome_sequence $single From_laboratory = HX AND Finished AND DNA\ngif EMBL $outfilename\n";# find Genome_sequences and EMBL dump
  $command .= "quit\nn\n";# say you don't want to save and exit
}

open (READ, "echo '$command' | $giface $dbdir |") or die ("Could not open $giface $dbdir\n"); 
while (<READ>) {
    next if ($_ =~ /\/\//);
    next if ($_ =~ /acedb/);
}
close (READ);


#########################
# Get COMMONDATA hashes #
#########################

# This is needed to fix missing sequence version and accession info in dumped EMBL files
# Where to get data from? If the build has been released COMMON_DATA has been removed, 
# Therefore need to be able to specify the latest/last build COMMON_DATA.

if ($version) {
$dir = $wormbase->database("WS$version")."/COMMON_DATA";
print "***$dir\n"
}

else {
  $dir = $wormbase->common_data;
}

#CommonData hash Key: Genome sequence Value: Sequence version integer
my %clone2accession = $wormbase->FetchData('clone2accession', undef, "$dir");
#CommonData hash Key: Clone/Sequence name Value: Type information (cosmid|fosmid|yac|Other);
my %clone2type = $wormbase->FetchData('clone2type', undef, "$dir");
my %cds2cgc  = $wormbase->FetchData('cds2cgc', undef, "$dir");
#my %cds2gene = $wormbase->FetchData('cds2wbgene_id');


######################################################################
# cycle through the EMBL dump file, replacing info where appropriate #
######################################################################

open (OUT, ">$mod_file") or  die "Can't process new EMBL dump file\n";
open (EMBL, "<$outfilename.embl") or die "Can't process EMBL dump file\n";

my $id = "";
my $cds;
my $clone;
my $ID2;

our $reference_remove = 0;
our $author_change;

while (<EMBL>) {
  # Store the necessary default ID line elements ready for use in the new style EMBL ID lines.
  if(/^ID\s+CE(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/){
    $ID2 = "XXX; linear; genomic $3 STD; "."$4 "."$5 "."$6";
    $id = $1;
    next;
  }
  # print new format ID line and AC lines with XX lines once the accession lookup has been done.
  if( /^AC/ ) {
    print "$_\n";
    print OUT "ID   $clone2accession{$id}; $ID2\nXX\n";
    print OUT "AC   $clone2accession{$id};\nXX\n";
    next;
  }
  # DE   Caenorhabditis elegans cosmid C05G5    
  if (/^DE   Caenorhabditis elegans cosmid (\S+)/) {
    $clone = $1;
    # can now reset $id
    $id = "";

    if (!defined($clone2type{$clone})){
      print OUT "DE   Caenorhabditis elegans clone $clone\n";
      print "WARNING: no clone type for $_";
      $log->write_to("WARNING: no clone type for $_");
    }
    elsif ($clone2type{$clone} eq "other") {
      print OUT "DE   Caenorhabditis elegans clone $clone\n";
    }
    elsif ($clone2type{$clone} eq "yac") {
      print OUT "DE   Caenorhabditis elegans YAC $clone\n";
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
    #      print OUT "RX   MEDLINE; 99069613.\n"; # Stripped by EMBL
    print OUT "RX   PUBMED; 9851916.\n";
    print OUT "RG   WormBase Consortium\n";
    print OUT "RA   The C. elegans Sequencing Consortium;\n"; may be required?
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

print "\nOutfile is $outfilename.*\n\n";

$log->mail("$maintainers");

exit(0);

__END__
