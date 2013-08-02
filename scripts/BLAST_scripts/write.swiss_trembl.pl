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

my ($swiss, $debug, $database, $list, $old, $wormpipe_dump);
my ($test, $store);
# $old is for switch back to protein model
GetOptions (
	    "swiss"      => \$swiss,
	    "database:s" => \$database,
	    "list"       => \$list,
	    "old"        => \$old,
	    "dumpdir:s"  => \$wormpipe_dump,
	    "debug:s"    => \$debug,
	    "test"       => \$test,
	    "store:s"    => \$store,
	  );

$wormpipe_dump ||= $ENV{'PIPELINE'}.'/dumps';

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
my $output_swiss      = "$wormpipe_dump/swissproteins.ace";
my $output_trembl     = "$wormpipe_dump/tremblproteins.ace";
my $swiss_list_txt    = "$wormpipe_dump/swisslist.txt";

my $db_files        = "$ENV{'PIPELINE'}/swall_data";

my $blast_files;
if (-e "/software/worm") {
  $blast_files = "$wormpipe_dump/*blastp.ace $wormpipe_dump/*X*.ace ";
} else {
  $blast_files = "$wormpipe_dump/*blastp.ace $wormpipe_dump/*x*.ace ";
}


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
  my %swisslist;
  while (<DATA>) {
    if (/Pep_homol\s+\"/) {
      if( /SW:(\S+)\"/ ) {
	$swisslist{$1} = 1;
      }
    }
  }
  close DATA;
  
  open (SWISS,">$swiss_list_txt");
  foreach (keys %swisslist) { print SWISS "$_\n"; }
  
  close SWISS;
}
# now extract info from dbm files and write ace files

$db_files = "$database" if defined $database;  # can use other database files if desired

my %input2output;
$input2output{"$swiss_list_txt"}  = [ ("$output_swiss", "SwissProt", "SW", "$db_files/slimswissprot" ) ];

my @lists_to_dump;
$db_files = "$database" if defined $database;  # can use other database files if desired

my %ORG;
my %DES;

unless ($swiss ) {
    die "usage -swiss for swissprot, -database directory where .gsi database file is\n";
}

my %file_mapping = ( 
	"$db_files/swissprot2org" => '/tmp/swissprot2org',
	"$db_files/swissprot2des" => '/tmp/swissprot2des',
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

$log->mail;
exit(0);

sub output_list
  {
    my $list = shift;
    my $parameters = shift;
    my $db = $$parameters[1];
    my $fasta_file = $$parameters[3];
    my $fetch_db = "${fasta_file}.gsi";
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
	open (DB , "$fasta_file") || die "cannot open db $fasta_file\n";
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

-debug  read/write to different directory

-database specify a non-default .gsi directory to use

=back

This script creates ace files containing the details of any proteins that have blast matches.

The input list are simply a list of ID's of matching proteins.

The .gsi databases are written by fasta2gsi.pl -f /scratch/ensembl/wormpipe/swall_data/slimswissprot whenever the swissprot are updated.
