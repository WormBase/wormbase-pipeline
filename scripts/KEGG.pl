#!/software/bin/perl -w

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Storable; 
use Log_files;
use DBI;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species, $no_load);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "no_load"    => \$no_load
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species
			     );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

# MYSQL CONFIG VARIABLES
my $dbhost = "cbi3";
my $dbuser = "genero";
my $dbname = &get_latest_uniprot;

$log->write_to("connecting to $dbhost:$dbname as $dbuser . . \n");
my $dbh = DBI -> connect("DBI:mysql:$dbname:$dbhost", $dbuser, "",{RaiseError => 1});
my $tax_id = $wormbase->ncbi_tax_id;

my $query =  "select dbxref.entry_id, dbxref.database_id, dbxref.primary_id, dbxref.tertiary_id from taxonomy,dbxref where taxonomy.ncbi_tax_id = \"".$tax_id."\" and taxonomy.entry_id = dbxref.entry_id and (database_id = \"Enzyme\" or database_id = \"WormBase\");";

$log->write_to("\tquerying . . \n");
my $sth = $dbh->prepare($query);
$sth->execute;
my %enzyme;

while (my ($entry, $db, $enzid, $geneid) = $sth->fetchrow()) {
    $enzyme{$entry}->{'gene'} = $geneid if($db eq 'WormBase');
    $enzyme{$entry}->{'kegg'} = $enzid if($db eq 'Enzyme');
}
$sth->finish;

my $out = $wormbase->acefiles."/KEGG.ace";
open (KEGG,">$out") or $log->log_and_die("cant write $out : $!\n");
$log->write_to("writing to $out\n");
foreach my $entry (keys %enzyme) {
    if( $enzyme{$entry}->{'gene'} and $enzyme{$entry}->{'kegg'}) {
	print KEGG "\nGene : ",$enzyme{$entry}->{'gene'},"\n","Database NemaPath KEGG_id \"",$enzyme{$entry}->{'kegg'}."\"\nDatabase KEGG KEGG_id \"",$enzyme{$entry}->{'kegg'},"\"\n";
    }
}
close KEGG;
$log->write_to("finished writing acefile\n");

unless ($no_load) {
    $log->write_to("\tloading to ".$wormbase->orgdb."\n");
    $wormbase->load_to_database($wormbase->orgdb, $out, 'KEGG', $log);
}

$log->write_to("DONE\n");

$log->mail;
exit;


sub get_latest_uniprot {
    $log->write_to("Getting latest uniprot table from M&M database . .\n");
    my $dbh = DBI -> connect("DBI:mysql:information_schema:$dbhost", $dbuser, "",{RaiseError => 1});
    my $query = "show databases";
    my $sth = $dbh->prepare($query);
    $sth->execute;
    my $up;
    my $year=0;
    my $month=0;
    while (my $db = $sth->fetchrow()) {
      if ($db =~ /uniprot_(\d+)_(\d+)/) {
	if ($1 > $year) {
	  $year = $1;
	  $month = $2;
	  $up = $db;
	} elsif ($1 == $year && $2 > $month) {
	  $month = $2;
	  $up = $db;
	}

      }
    }
    $log->write_to("\tusing UniProt version $up\n");
    return $up;
}
