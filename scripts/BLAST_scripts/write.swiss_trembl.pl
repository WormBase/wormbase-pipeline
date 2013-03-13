#!/usr/local/ensembl/bin/perl

# Marc Sohrmann (ms2@sanger.ac.uk)

# takes as input a swissprot or trembl .fasta file,
# and deletes all worm, fly, human and yesast entries

use lib $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'}."/BLAST_scripts";
use strict;
use Getopt::Long;
use DB_File;
use GSI;
use Wormbase;
use File::Copy;

my ($swiss, $trembl, $debug, $database, $list, $old);
my ($test, $store);
# $old is for switch back to protein model
GetOptions (
	    "swiss"      => \$swiss,
	    "trembl"     => \$trembl,
	    "database:s" => \$database,
	    "list"       => \$list,
	    "old"        => \$old,
	    "debug:s"    => \$debug,
	    "test"       => \$test,
	    "store:s"    => \$store
	  );

my $wormbase;
if ( $store ) {
  $wormbase = Storable::retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

my $log = Log_files->make_build_log($wormbase);

my $acefiles  = $wormbase->acefiles;
my $wormpipe_dump     = $wormbase->farm_dump;
my $output_swiss      = "$wormpipe_dump/swissproteins.ace";
my $output_trembl     = "$wormpipe_dump/tremblproteins.ace";
my $swiss_list_txt    = "$wormpipe_dump/swisslist.txt";
my $trembl_list_txt   = "$wormpipe_dump/trembllist.txt";

my $db_files        = "$ENV{'PIPELINE'}/swall_data";

my $blast_files = "$wormpipe_dump/*blastp.ace $wormpipe_dump/*X*.ace ";

# for species not being built, the blastp files will live in the build
# directory
my %acc = $wormbase->species_accessors;
foreach my $sp (keys %acc) {
  my $blastp_file = $acc{$sp}->acefiles . "/${sp}_blastp.ace";

  if (-e $blastp_file) {
    $blast_files .= " $blastp_file";
  }
  
  my $blastx_file = $acc{$sp}->acefiles . "/${sp}_blastx.ace";
  if (-e $blastx_file) {
    $blast_files .= " $blastx_file";
  }
}

# extract and write lists of which proteins have matches
unless ( $list ){
  print "Grabbing proteins from these files: $blast_files\n";

  open (DATA,"cat $blast_files |");
  my (%swisslist, %trembllist);
  while (<DATA>) {
    if (/Pep_homol\s+\"/) {
      if( /SW:(\S+)\"/ ) {
	$swisslist{$1} = 1;
      }
      elsif( /Pep_homol\s+\"TR:(\S+)\"/ ) {
	$trembllist{$1} = 1;
      }
    }
  }
  close DATA;
  
  open (SWISS,">$swiss_list_txt");
  open (TREMBL,">$trembl_list_txt");
  foreach (keys %swisslist) { print SWISS "$_\n"; }
  foreach (keys %trembllist) { print TREMBL "$_\n"; }
  
  close SWISS;
  close TREMBL;
}
# now extract info from dbm files and write ace files

$db_files = "$database" if defined $database;  # can use other database files if desired

my %input2output;
$input2output{"$swiss_list_txt"}  = [ ("$output_swiss", "SwissProt", "SW", "$db_files/slimswissprot" ) ];
$input2output{"$trembl_list_txt"} = [ ("$output_trembl", "TrEMBL", "TR", "$db_files/slimtrembl_f" ) ];

my @lists_to_dump;
$db_files = "$database" if defined $database;  # can use other database files if desired

my %ORG;
my %DES;

unless ($swiss || $trembl) {
    die "usage -swiss for swissprot, -trembl for trembl, -database directory where .gsi database file is\n";
}

my %file_mapping = ( 
	"$db_files/swissprot2org" => '/tmp/swissprot2org',
	"$db_files/trembl2org"    => '/tmp/trembl2org',
	"$db_files/swissprot2des" => '/tmp/swissprot2des',
	"$db_files/trembl2des"    => '/tmp/trembl2des',
);
while (my($from,$to)=each %file_mapping ){
	unlink $to if -e $to;
	copy $from,$to or die "Copy of $from to $to failed: $!";
}

if ($swiss) {
  unless (-s "/tmp/swissprot2org") {
    die "swissprot2org not found or empty";
  }
  tie (%ORG, 'DB_File', "/tmp/swissprot2org", O_RDWR|O_CREAT, 0666, $DB_HASH) or $log->log_and_die("cannot open /tmp/swissprot2org DBM file\n");
  unless (-s "/tmp/swissprot2des") {
    die "swissprot2des not found or empty";
  }
  tie (%DES, 'DB_File', "/tmp/swissprot2des", O_RDWR|O_CREAT, 0666, $DB_HASH) or $log->log_and_die("cannot open /tmp/swissprot2des DBM file\n");
  &output_list($swiss_list_txt,$input2output{"$swiss_list_txt"});
  untie %ORG;
  untie %DES;
}

if ($trembl) {
  unless (-s "/tmp/trembl2org") {
    die "trembl2org not found or empty";
  }
  tie (%ORG, 'DB_File', "/tmp/trembl2org", O_RDWR|O_CREAT, 0666, $DB_HASH) or $log->log_and_die("cannot open /tmp/trembl2org DBM file\n");
  unless (-s "/tmp/trembl2des") {
    die "trembl2des not found or empty";
  }
  tie (%DES, 'DB_File', "/tmp/trembl2des", O_RDWR|O_CREAT, 0666, $DB_HASH) or $log->log_and_die("cannot open /tmp/trembl2des DBM file\n");
  push( @lists_to_dump,$trembl_list_txt);
  &output_list($trembl_list_txt,$input2output{"$trembl_list_txt"});
  untie %ORG;
  untie %DES;
}

$log->mail;
exit(0);
    

sub output_list
  {
    my $list = shift;
    my $parameters = shift;
    my $db = $$parameters[1];
    my $fetch_db = $$parameters[3].".gsi";
    my $prefix = $$parameters[2];
    my $outfile = $$parameters[0];
    
    #used to be in fetch.pl
    ########################
    # open the GSI file
    GSI::openGSI ("$fetch_db");
    ########################


    open (LIST,"<$list") or die "cant open input file $list\n";
    open (ACE, ">$outfile") or die "cant write to $outfile\n";
    while (<LIST>) {
      #	print;  # ID line 
      chomp;
      /^(\S+)/;

      my $id = $1;
      
      # access gsi database to get info about protein
      my ($file , $fmt , $offset) = GSI::getOffset ($id);
      unless( "$file" eq "-1" ) {
	open (DB , "$file") || die "cannot open db $file\n";
	seek (DB , $offset , 0);
	my $header = <DB>;
	chomp $header;
	my $seq = "";
	while ((my $line = <DB>) =~ /^[^\>]/) {
	  $seq .= "$line"."\n";
	}
	close DB;
	$seq =~ s/\n//g;
	
	$header =~ /^>(\S+)/;
	my $accession = $1;
	
	if ($seq) {
	  # that is a cleanup for already existing dodgy description lines
	  my $description = $DES{$accession}=~/Full=([^;]*)/?$1:$DES{$accession};

	  print ACE "Protein : \"$prefix:$accession\"\n";
	  print ACE "Peptide \"$prefix:$accession\"\n";
	  print ACE "Species \"$ORG{$accession}\"\n";
	  print ACE "Description \"$description\"\n";
	  print ACE "Database UniProt UniProtAcc $accession\n";
	  print ACE "Database UniProt UniProtID $id\n" if $id;
	  print ACE  "\n";
	  print ACE  "Peptide : \"$prefix:$accession\"\n";
	  print ACE "$seq\n";
	  print ACE "\n";
	}
	else {
	  print "// Couldn't fetch sequence for $accession\n\n";
	}
      }  
    }
  }


__END__

=pod

=head2   NAME - write.swiss_trembl.pl

=head1 USAGE

=over 4

=item write.swiss_trembl.pl [-options]

-swiss  get SwissProt entries

-trembl get TrEMBL  entries

-debug  read/write to different directory

-database specify a non-default .gsi directory to use

=back

This script creates ace files containing the details of any proteins that have blast matches.

The input list are simply a list of ID's of matching proteins.

The .gsi databases are written by fasta2gsi.pl -f /scratch/ensembl/wormpipe/swall_data/slimswissprot whenever the swissprot/trembl are updated.
