#!/usr/local/bin/perl5.8.0 -w
#
# EMBLdump.pl :  makes modified EMBL dumps from camace.
# 
#  Last updated on: $Date: 2011-03-04 14:24:55 $
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
my $database;
GetOptions (
	    "test"         => \$test,
	    "debug=s"      => \$debug,     # Only emails specified recipient and turns on extra printing.
	    "store:s"      => \$store,
	    "quicktest"    => \$quicktest, # only dumps B0250
	    "single=s"     => \$single,    # only dumps specified clone
	    "version=s"    => \$version,   # Specifies what WS version COMMON_DATA to use.
	    "database=s"   => \$database,  # Requires full path to database eg.c ~wormpub/DATABASES/BACKUPS/camace_backup.061116
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
  $log->write_to("DEBUG = \"$debug\"\n\n");
  ($maintainers = $debug . '\@sanger.ac.uk');
}


###############################
# misc. variables             #
###############################

my $basedir = $wormbase->wormpub;
my $giface = $wormbase->giface;
my $dbdir;
if ($database) {$dbdir = $database;}
else {$dbdir = $wormbase->database('camace');}
my $refdb = $wormbase->database('current');
my $tace = $wormbase->tace;

my $outfilename = "$basedir/tmp/EMBLdump.$$";
my $mod_file = "$basedir/tmp/EMBLdump.mod";
my $dir;
my $flag;
my %cds2status;
my %cds2cgc;
my %clone2type;
my %clone2accession;
my %specialclones;
my %rnagenes;
#These are clones that were swapped with St Louis
$specialclones{'CU457737'} = "AC006622";
$specialclones{'CU457738'} = "AC006690";
$specialclones{'CU457739'} = "AF043691";
$specialclones{'CU457740'} = "AF098988";
$specialclones{'CU457741'} = "AF043695";
$specialclones{'CU457742'} = "AF045637";
$specialclones{'CU457743'} = "AF003149";
$specialclones{'CU457744'} = "AF099926";

$log->write_to("You are embl dumping from $dbdir\n\n");

#Fetch additional info#
&fetch_database_info;

#############################################
# Use giface to dump EMBL files from camace #
#############################################

my $command;
if ($quicktest) {  # quicktest mode only works on B0250
  $single = "B0250";
}

if ($single) {
  $command  = "nosave\n"; # Don't really want to do this
  $command .= "query find CDS where Method = \"Genefinder\"\nkill\ny\n";# remove Genefinder predictions
  $command .= "query find CDS where Method = \"twinscan\"\nkill\ny\n";# remove twinscan predictions
  $command .= "query find CDS where Method = \"jigsaw\"\nkill\ny\n";# remove jigsaw predictions
  $command .= "query find Genome_sequence $single From_laboratory = HX AND Finished AND DNA\ngif EMBL $outfilename\n";# find sequence and dump
  $command .= "quit\nn\n";# say you don't want to save and exit
}

else {
  $command  = "nosave\n"; # Don't really want to do this
  $command .= "query find CDS where Method = \"Genefinder\"\nkill\ny\n";# remove Genefinder predictions
  $command .= "query find CDS where Method = \"twinscan\"\nkill\ny\n";# remove twinscan predictions
  $command .= "query find CDS where Method = \"jigsaw\"\nkill\ny\n";# remove jigsaw predictions
  $command .= "query find Genome_sequence From_laboratory = HX AND Finished AND DNA\ngif EMBL $outfilename\n";# find sequence and dump
  $command .= "quit\nn\n";# say you don't want to save and exit
}

$log->write_to("$command\n");
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
  $log->write_to("***$dir\n");
}

else {
  $dir = $wormbase->database("current")."/COMMON_DATA";
}

#CommonData hash Key: Genome sequence Value: Sequence version integer
if (!(-e $dir."/clone2accession.dat")) {
  $dir = $wormbase->common_data;
  $log->write_to("Switching to the build copy of COMMON_DATA $dir!!!\n");
}

%clone2accession = $wormbase->FetchData('clone2accession', undef, "$dir");
#CommonData hash Key: Clone/Sequence name Value: Type information (cosmid|fosmid|yac|Other);
%clone2type = $wormbase->FetchData('clone2type', undef, "$dir");
%cds2cgc  = $wormbase->FetchData('cds2cgc', undef, "$dir");
#my %cds2gene = $wormbase->FetchData('cds2wbgene_id', undef "$dir");
%rnagenes  = $wormbase->FetchData('rna2cgc', undef, "$dir");



######################################################################
# cycle through the EMBL dump file, replacing info where appropriate #
######################################################################

open (OUT, ">$mod_file") or  die "Can't process new EMBL dump file\n";
open (EMBL, "<$outfilename.embl") or die "Can't process EMBL dump file\n";

my $id = "";
my $cds;
my $locus;
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
    my $acc = $clone2accession{$id};
    print OUT "ID   $clone2accession{$id}; $ID2\nXX\n";
#    next;
#AC * _AC006622
#AC   CU457737; AC006622;

    if (defined$specialclones{$acc}) {
      print OUT "AC * _$specialclones{$acc}\n";
      print OUT "AC   $clone2accession{$id}; $specialclones{$acc};\nXX\n";
    }
    else { 
      print OUT "AC   $clone2accession{$id};\nXX\n";
    }
    next;
  }
  # DE   Caenorhabditis elegans cosmid C05G5    
  if (/^DE   Caenorhabditis elegans cosmid (\S+)/) {
    $clone = $1;
    # can now reset $id
    $id = "";
    
    if (!defined($clone2type{$clone})){
      print OUT "DE   Caenorhabditis elegans clone $clone\n";
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
    print OUT "RA   Caenorhabditis elegans Sequencing Consortium;\n";
    print OUT "RT   \"Genome sequence of the nematode C. elegans: a platform for investigating\n";
    print OUT "RT   biology\";\n";
    print OUT "RL   Science 282(5396):2012-2018(1998).\n";
    print OUT "XX\n";
    print OUT "RN   [2]\n";
    print OUT "$_";
    next;
  }
  if (/^RN   \[2\]/) {
    $reference_remove = 5;
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
    $locus = $1;
    #strip the isoform letter
    if ($locus =~ /(\S+\.\d+)\l/){
#    if ($locus =~ /(\S+\.\d+)/){
      $cds = $1;
    }
    else {
      $cds = $locus;
    }
    
    if ((defined$cds2cgc{$cds}) || (defined$cds2cgc{$locus}) ) {
      print OUT "FT                   /gene=\"" . $cds2cgc{$locus}  ."\"\n" if defined($cds2cgc{$locus});
      print OUT "FT                   /gene=\"" . $cds2cgc{$cds}  ."\"\n" if (($cds ne $locus) && defined($cds2cgc{$cds}));
      print OUT "FT                   /locus_tag=\"$locus\"\n";
      next unless defined($flag);#allows non-coding_transcript isoform addition
    }
    elsif ((defined$rnagenes{$cds}) || (defined$rnagenes{$locus})) {
      print OUT "FT                   /gene=\"" . $rnagenes{$cds}  ."\"\n" if defined($rnagenes{$cds});
      print OUT "FT                   /gene=\"" . $rnagenes{$locus}  ."\"\n" if ($cds ne $locus);
      print OUT "FT                   /locus_tag=\"$locus\"\n";
      next unless defined($flag);#allows non-coding_transcript isoform addition
    }
    else {
      print OUT "FT                   /locus_tag=\"$locus\"\n";
      next unless defined($flag);#allows non-coding_transcript isoform addition
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

##################################################################################################################
   #  Modify the /product lines for genes that have EST evidence
   # /product="Hypothetical protein XXX.X"
   # /product="C. elegans protein XXX.X, partially confirmed by transcript evidence"
   # /product="C. elegans protein XXX.X, confirmed by transcript evidence"
  if (/\/product=\"Hypothetical protein (\S+.\S+)\"/) {
    $cds = $1;
    if (!defined$cds2status{$cds}) {
      $log->write_to("Warning data missing for $cds\n$_\n");
      print OUT "$_";
      next;
    }
    if ($cds2status{$cds} eq ("Confirmed")) {
      print OUT "FT                   \/product=\"C. elegans protein $cds,\nFT                   confirmed by transcript evidence\"\n";
      next;
    }
    elsif ($cds2status{$cds} eq ("Partially_confirmed")) {
      print OUT "FT                   \/product=\"C. elegans protein $cds,\nFT                   partially confirmed by transcript evidence\"\n";
      next;
    }
    else {
      print OUT "FT                   \/product=\"C. elegans predicted protein $cds\"\n";
      next;
    }
  }

 # Modify Product lines for RNA genes
  if (/\/product=\"Hypothetical RNA transcript (\S+.\S+)\"/) {
    my $name = $1;
    if (defined$rnagenes{$name}) {
      print OUT "FT                   \/product=\"C. elegans RNA transcript $name ($rnagenes{$name})\"\n";
      next;
    }
    else {
      print OUT "FT                   \/product=\"C. elegans RNA transcript $name\"\n";
      next;
    }
  }
#############################################################################################################
  ###########################################################
  ## Feature Table edits that are necessary for submission.##
  ###########################################################


  # RNA line edits for ncRNA type data.
  # All RNAs should be converted if they aren't tRNA or rRNA
  #
  # 1st FT line should be one of 3
  # FT    ncRNA
  # FT    rRNA
  # FT    tRNA
  # Supported bio types for ncRNA
  #  /ncRNA_class="miRNA"
  #  /ncRNA_class="siRNA"
  #  /ncRNA_class="scRNA"              
  #  /ncRNA_class="Other"
  #  /ncRNA_class="snoRNA"
  #  /ncRNA_class="snRNA"
  # Nothing else counts
  
  if ((/^FT\s+(\w{1}RNA)\s+/) || (/^FT\s+(\w{2}RNA)\s+/) || (/^FT\s+(\w{3}RNA)\s+/) || (/^FT\s+(misc_RNA)\s+/)) {
    my $mol = $1;

    # Supported bio-types
    if (($mol eq "tRNA") || ($mol eq "rRNA")) {
      print OUT "$_";
      next;
    }
    elsif (($mol eq "snoRNA") || ($mol eq "miRNA") || ($mol eq "siRNA") || ($mol eq "scRNA") || ($mol eq "snRNA")) {
      if ($mol =~ (/^\w{2}RNA/)) {
	s/$mol/ncRNA/g;
      }
      elsif ($mol =~ (/^\w{3}RNA/)) {
	s/$mol/ncRNA /g;
      }
      print OUT "$_";
      print OUT "FT                   /ncRNA_class=\"$mol\"\n";
      next;
    }
    
    # Un-supported bio-types
    elsif (($1 eq "ncRNA") || ($1 eq "misc_RNA")) {
      if ($1 eq "misc_RNA") {
	s/$1/ncRNA   /g;
      }
      else {
	s/$1/ncRNA/g;
      }
      print OUT "$_";
      print OUT "FT                   /ncRNA_class=\"Other\"\n" unless (/,/);
      $flag = "$_" if (/,/);
      next;
    }
    elsif ($1 eq "snlRNA") {
      s/$1/ncRNA /g;
      print OUT "$_";
      print OUT "FT                   /ncRNA_class=\"Other\"\n";
      print OUT "FT                   /note=\"$mol\"\n";
      next;
    }
    else {
      print OUT "RNA LINE ERROR:$_\n";
      next;
    }
  }

  #non_coding_transcript Isoform notes section
  #FT                   /gene="
  # this is a bit messy as it relies on a skip of the next when converting /gene=" up above.
  # This might be a bit inclusive as it works on all Transcripts with complexed structure.
  if ((/\/gene/) && defined($flag)){
    print OUT "FT                   /ncRNA_class=\"Other\"\n";
    print OUT "FT                   /note=\"non-functional Isoform from coding locus\"\n";
    undef $flag;
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
# It's useful to keep a copy when debugging
if ($debug) {
  $log->write_to("debug selected therefore a copy of the original embl dump will be made for backup \n");
  system ("mv $outfilename.embl $outfilename.bk");
}

my $status = move("$mod_file","$outfilename.embl");
$log->write_to("ERROR: Couldn't move file: $!\n") if ($status == 0);

$log->write_to("\nOutfile is $outfilename.*\n\n");

$log->mail();

exit(0); 

###################################################
#                 SUBROUTINES                     #
###################################################


#######################
# fetch_database_info #
#######################

sub fetch_database_info {
  my $query;
  my $def_dir = $wormbase->database('current')."/wquery";
  my @queries = ("${def_dir}/SCRIPT:CDS_status.def");
  foreach $query(@queries) {
    my $command = "Table-maker -p $query\nquit\n";
    
    open (TACE, "echo '$command' | $tace $refdb |");
    while (<TACE>) {
      #my $status;
      chomp;
      s/\"//g;
      next unless (/^([A-Z,0-9,.]+?\w)\s+(\w+)/) ;
      my $cds = $1;
      my $status = $2;
      $cds2status{$cds} = "$status";
    }
    close TACE;
  }
}

__END__
