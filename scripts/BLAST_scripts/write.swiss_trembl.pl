#!/usr/local/ensembl/bin/perl

# Marc Sohrmann (ms2@sanger.ac.uk)

# takes as input a swissprot or trembl .fasta file,
# and deletes all worm, fly, human and yesast entries

use lib $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'}."/BLAST_scripts";
use strict;
use Getopt::Long;
use GDBM_File;
use GSI;
use Wormbase;

my ($swiss, $trembl, $debug, $database, $list, $old, $species);
my ($test, $store);
# $old is for switch back to protein model
GetOptions (
	    "swiss"      => \$swiss,
	    "trembl"     => \$trembl,
	    "database:s" => \$database,
	    "list"       => \$list,
	    "old"        => \$old,
	    "debug:s"    => \$debug,
	    "species=s"  => \$species,
	    "test"       => \$test,
	    "store:s"    => \$store
	  );

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
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
my $blastx_file       = "$wormpipe_dump/blastx_ensembl.ace";
my $blastp_file       = "$wormpipe_dump/blastp_ensembl.ace";
my $ensembl_info_file = "$wormpipe_dump/ensembl_protein_info.ace";
my $swiss_list_txt    = "$wormpipe_dump/swisslist.txt";
my $trembl_list_txt   = "$wormpipe_dump/trembllist.txt";

my $db_files        = "/acari/work2a/wormpipe/swall_data";

my @blastp_databases = qw(worm_pep worm_brigpep);
my $blast_files = "$wormpipe_dump/*blastp.ace $wormpipe_dump/*blastx.ace ";

# extract and write lists of which proteins have matches
unless ( $list ){
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

my @lists_to_dump;
$db_files = "$database" if defined $database;  # can use other database files if desired

my %input2output;
$input2output{"$swiss_list_txt"}  = [ ("$output_swiss", "SwissProt", "SW", "$db_files/slimswissprot" ) ];
$input2output{"$trembl_list_txt"} = [ ("$output_trembl", "TrEMBL", "TR", "$db_files/slimtrembl" ) ];

my @lists_to_dump;
$db_files = "$database" if defined $database;  # can use other database files if desired

my %ORG;
my %DES;

unless ($swiss || $trembl) {
    die "usage -swiss for swissprot, -trembl for trembl, -database directory where .gsi database file is\n";
}

if ($swiss) {
  unless (-s "$db_files/swissprot2org") {
    die "swissprot2org not found or empty";
  }
  tie %ORG,'GDBM_File', "$db_files/swissprot2org",&GDBM_WRCREAT, 0666 or die "cannot open swissprot2org DBM file $db_files/swissprot2org";
  unless (-s "$db_files/swissprot2des") {
    die "swissprot2des not found or empty";
  }
  tie %DES,'GDBM_File', "$db_files/swissprot2des",&GDBM_WRCREAT, 0666 or die "cannot open swissprot2des DBM file $db_files/swissprot2des";
  &output_list($swiss_list_txt,$input2output{"$swiss_list_txt"});
  untie %ORG;
  untie %DES;
}

if ($trembl) {
  unless (-s "$db_files/trembl2org") {
    die "trembl2org not found or empty";
  }
  tie %ORG,'GDBM_File', "$db_files/trembl2org",&GDBM_WRCREAT, 0666 or die "cannot open trembl2org DBM file";
  unless (-s "$db_files/trembl2des") {
    die "trembl2des not found or empty";
  }
  tie %DES,'GDBM_File', "$db_files/trembl2des",&GDBM_WRCREAT, 0666 or die "cannot open trembl2des DBM file";
  push( @lists_to_dump,$trembl_list_txt);
  &output_list($trembl_list_txt,$input2output{"$trembl_list_txt"});
  untie %ORG;
  untie %DES;
}

system("cat $output_swiss $output_trembl > $ensembl_info_file");

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
	  print ACE "Protein : \"$prefix:$accession\"\n";
	  print ACE "Peptide \"$prefix:$accession\"\n";
	  print ACE "Species \"$ORG{$accession}\"\n";
	  print ACE "Description \"$DES{$accession}\"\n";
	  if ("$prefix" eq "SW" ) {
	    #print ACE "Database SwissProt SwissProt_ID $id\n";
	    print ACE "Database SwissProt SwissProt_AC $accession\n";
	  }
	  else {
	    print ACE "Database TREMBL TrEMBL_AC $id\n";
	  }
	
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

The .gsi databases are written by fasta2gsi.pl -f /lustre/work1/ensembl/wormpipe/swall_data/slimswissprot whenever the swissprot/trembl are updated.
