#!/usr/local/bin/perl5.8.0 -w
#
# EMBLdump.pl :  makes modified EMBL dumps from camace.
# 
#  Last updated on: $Date: 2008-05-09 14:38:19 $
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
  print "DEBUG = \"$debug\"\n\n";
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
$specialclones{'CU457737'} = "AC006622";
$specialclones{'CU457738'} = "AC006690";
$specialclones{'CU457739'} = "AF043691";
$specialclones{'CU457740'} = "AF098988";
$specialclones{'CU457741'} = "AF043695";
$specialclones{'CU457742'} = "AF045637";
$specialclones{'CU457743'} = "AF003149";
$specialclones{'CU457744'} = "AF099926";

my %mirgenes; #this is bad as all are hard coded.
$mirgenes{'F09A5.7'}  = "miRNA (mir-56) secondary product";
$mirgenes{'T16G12.11'}  = "mir-2120";
$mirgenes{'Y51H4A.402'}  = "mir-2021";
$mirgenes{'Y41E3.214'}  = "mir-1833";
$mirgenes{'F39B1.2'}  = "mir-1829c";
$mirgenes{'F20D1.11'}  = "mir-1829b";
$mirgenes{'K09A9.7'}  = "mir-1829a";
$mirgenes{'T22A3.9'}  = "mir-1828";
$mirgenes{'Y87G2A.21'}  = "mir-1824";
$mirgenes{'T07D10.9'}  = "mir-1823";
$mirgenes{'R11.5'}  = "mir-1823";
$mirgenes{'Y71A12B.20'}  = "mir-1818";
$mirgenes{'M04C9.8'}  = "mir-1019";
$mirgenes{'F33E2.9'}  = "mir-795";
$mirgenes{'T07D10.7'}  = "mir-794";
$mirgenes{'M03B6.6'}  = "mir-793";
$mirgenes{'T28F3.10'}  = "mir-789.2";
$mirgenes{'Y51H4A.34'}  = "mir-789.1";
$mirgenes{'F08G12.12'}  = "mir-787";
$mirgenes{'F54B11.12'}  = "mir-392";
$mirgenes{'T27D12.5'}  = "mir-355";
$mirgenes{'Y105E8A.31'}  = "mir-354";
$mirgenes{'E01F3.2'}  = "mir-273";
$mirgenes{'Y41C4A.20'}  = "mir-272";
$mirgenes{'C06H2.8'}  = "mir-268";
$mirgenes{'Y38E10A.27'}  = "mir-267";
$mirgenes{'T23F6.6'}  = "mir-265";
$mirgenes{'ZK384.5'}  = "mir-262";
$mirgenes{'F25D1.6'}  = "mir-259";
$mirgenes{'Y102A5D.4'}  = "mir-257";
$mirgenes{'ZK455.9'}  = "mir-254";
$mirgenes{'W02B12.14'}  = "mir-252";
$mirgenes{'F59F3.7'}  = "mir-251";
$mirgenes{'F55A11.12'}  = "mir-250";
$mirgenes{'ZK593.10'}  = "mir-246";
$mirgenes{'F55D12.7'}  = "mir-245";
$mirgenes{'F56A12.4'}  = "mir-241";
$mirgenes{'C34E11.6'}  = "mir-239.2";
$mirgenes{'C34E11.5'}  = "mir-239.1";
$mirgenes{'K01F9.5'}  = "mir-238";
$mirgenes{'C13B4.3'}  = "mir-234";
$mirgenes{'W03G11.5'}  = "mir-233";
$mirgenes{'F13H10.7'}  = "mir-232";
$mirgenes{'C29E6.7'}  = "mir-124";
$mirgenes{'K01F9.1'}  = "mir-90";
$mirgenes{'F10C2.8'}  = "mir-87";
$mirgenes{'Y56A3A.34'}  = "mir-86";
$mirgenes{'F49E12.11'}  = "mir-85";
$mirgenes{'B0395.4'}  = "mir-84";
$mirgenes{'K01F9.3b'}  = "mir-80";
$mirgenes{'K01F9.3a'}  = "mir-80";
$mirgenes{'C12C8.4'}  = "mir-79";
$mirgenes{'Y40H7A.12'}  = "mir-78";
$mirgenes{'T21B4.13'}  = "mir-77";
$mirgenes{'F16A11.4'}  = "mir-71";
$mirgenes{'T07C5.6'}  = "mir-62";
$mirgenes{'F55A11.9'}  = "mir-61";
$mirgenes{'B0035.17'}  = "mir-59";
$mirgenes{'T09A5.13'}  = "mir-57";
$mirgenes{'F09A5.8'}  = "mir-56";
$mirgenes{'F09A5.6'}  = "mir-55";
$mirgenes{'F09A5.5'}  = "mir-54";
$mirgenes{'F36H1.8'}  = "mir-53";
$mirgenes{'Y37A1B.16'}  = "mir-52";
$mirgenes{'F36H1.7'}  = "mir-51";
$mirgenes{'F19C6.6'}  = "mir-49";
$mirgenes{'K02B9.5'}  = "mir-47";
$mirgenes{'ZK525.3'}  = "mir-46";
$mirgenes{'ZK930.11'}  = "mir-45";
$mirgenes{'ZK930.10'}  = "mir-44";
$mirgenes{'ZK930.9'}  = "mir-43";
$mirgenes{'ZK930.8'}  = "mir-42";
$mirgenes{'Y62F5A.8'}  = "mir-41";
$mirgenes{'Y62F5A.7'}  = "mir-40";
$mirgenes{'Y62F5A.6'}  = "mir-39";
$mirgenes{'Y62F5A.5'}  = "mir-38";
$mirgenes{'Y62F5A.4'}  = "mir-37";
$mirgenes{'Y62F5A.3'}  = "mir-36";
$mirgenes{'Y62F5A.2'}  = "mir-35";
$mirgenes{'M04C9.7'}  = "mir-2";
$mirgenes{'C32C4.6'}  = "lsy-6";
$mirgenes{'F56A12.3'}  = "lin-58";
$mirgenes{'C05G5.6'}  = "let-7";

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
  $command .= "query find Genome_sequence $single From_laboratory = HX AND Finished AND DNA\ngif EMBL $outfilename\n";# find sequence and dump
  $command .= "quit\nn\n";# say you don't want to save and exit
}

else {
  $command  = "nosave\n"; # Don't really want to do this
  $command .= "query find CDS where Method = \"Genefinder\"\nkill\ny\n";# remove Genefinder predictions
  $command .= "query find CDS where Method = \"twinscan\"\nkill\ny\n";# remove twinscan predictions
  $command .= "query find Genome_sequence From_laboratory = HX AND Finished AND DNA\ngif EMBL $outfilename\n";# find sequence and dump
  $command .= "quit\nn\n";# say you don't want to save and exit
}

print "$command\n";
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
if (!(-e "/nfs/disk100/wormpub/BUILD/autoace/COMMON_DATA/clone2accession.dat")) {
print "Copying COMMON_DATA around as the build has finished!!!\n";
system ("scp ~/DATABASES/current_DB/COMMON_DATA/clone2accession.dat  /nfs/disk100/wormpub/BUILD/autoace/COMMON_DATA/");
system ("scp ~/DATABASES/current_DB/COMMON_DATA/clone2type.dat /nfs/disk100/wormpub/BUILD/autoace/COMMON_DATA/");
system ("scp ~/DATABASES/current_DB/COMMON_DATA/cds2cgc.dat /nfs/disk100/wormpub/BUILD/autoace/COMMON_DATA/");
}
%clone2accession = $wormbase->FetchData('clone2accession', undef, "$dir");
#CommonData hash Key: Clone/Sequence name Value: Type information (cosmid|fosmid|yac|Other);
%clone2type = $wormbase->FetchData('clone2type', undef, "$dir");
%cds2cgc  = $wormbase->FetchData('cds2cgc', undef, "$dir");
#my %cds2gene = $wormbase->FetchData('cds2wbgene_id', undef "$dir");


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
    my $acc = $clone2accession{$id};
    print OUT "ID   $clone2accession{$id}; $ID2\nXX\n";

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
    $cds = $1;
    if (defined$cds2cgc{$cds}) {
      print OUT "FT                   /gene=\"" . $cds2cgc{$cds}  ."\"\n";
      print OUT "FT                   /locus_tag=\"$cds\"\n";
#      print OUT "FT                   /note=\"non-functional Isoform from coding locus\"\n" if ($flag eq 1);
      next;
    }
    elsif (defined$mirgenes{$cds}) {
      print OUT "FT                   /gene=\"" . $mirgenes{$cds}  ."\"\n";
      print OUT "FT                   /ncRNA_class=\"Other\"\n";
      print OUT "FT                   /note=\"miRNA\"\n";
      print OUT "FT                   /locus_tag=\"$cds\"\n";
      next;
    }
    else {
      print OUT "FT                   /ncRNA_class=\"Other\"\n" if ($flag eq 1);
#      print OUT "FT                   /note=\"non-functional Isoform from coding locus\"\n" if ($flag eq 1);
      print OUT "FT                   /locus_tag=\"$cds\"\n";
      $flag = "0";
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

   #  Modify the /product lines for genes that have EST evidence
   # /product="Hypothetical protein XXX.X"
   # /product = "C. elegans protein XXX.X, partially confirmed by transcript evidence"
   # /product = "C. elegans protein XXX.X, confirmed by transcript evidence"
    if (/\/product\=\"Hypothetical protein (\S+.\S+)\"/) {
      $cds = $1;
      if ($cds2status{$cds} eq ("Confirmed")) {
	print OUT "FT                   \/product = \"C. elegans protein $cds,\nFT                   confirmed by transcript evidence\"\n";
	next;
      }
      elsif ($cds2status{$cds} eq ("Partially_confirmed")) {
	print OUT "FT                   \/product = \"C. elegans protein $cds,\nFT                   partially confirmed by transcript evidence\"\n";
	next;
      }
      else {
	print OUT "$_";
	next;
      }
    }

  ###########################################################
  ## Feature Table edits that are necessary for submission.##
  ###########################################################

  #FT   ncRNA(5 char -> 8 characters misc_RNA)
  #if (/^FT\s+(\w{2}RNA)\s+join/) {
  #  print OUT "$_";
  #  next;
  #}
  #if ((/^FT\s+(\w{2}RNA)/) || (/^FT\s+(\w{3}RNA)/) || (/^FT\s+(misc_RNA)/)) {
  if ((/^FT\s+(\w{2}RNA)\s+/) || (/^FT\s+(\w{3}RNA)\s+/) || (/^FT\s+(misc_RNA)\s+/)) {
    chomp;
    my $mol = $1;
    if ($1 eq "snlRNA") {
      s/$1/ncRNA /g;
      print OUT "$_\n";
      print OUT "FT                   /ncRNA_class=\"Other\"\n";
      print OUT "FT                   /note=\"$mol\"\n";
      next;
    }
    elsif ($1 eq "snoRNA") {
      s/$1/ncRNA /g;
      print OUT "$_\n";
      print OUT "FT                   /ncRNA_class=\"$mol\"\n";
      next;
    }
    elsif ($1 eq "ncRNA") {
      print OUT "$_\n";
      print OUT "FT                   /ncRNA_class=\"Other\"\n" unless (/,/);
      if (/,/) {$flag = "1";}
      next;
    }
    elsif (/^FT\s+(\w{2}RNA)/) {
      s/$1/ncRNA/g;
      print OUT "$_\n";
      print OUT "FT                   /ncRNA_class=\"$mol\"\n";
      next;
    }
    elsif (/^FT\s+(misc_RNA)/) {
      s/$1/ncRNA   /g;
      print OUT "$_\n";
      next;
    }
    else {
      print OUT "11$_\n";
      next;
    }
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
