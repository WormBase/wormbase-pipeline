#!/usr/bin/env perl 

use strict;
use lib $ENV{CVS_DIR};

use Getopt::Long;
use Storable; 

use Wormbase;
use Log_files;
use Modules::EBI_datafetch;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species, $noload, $acefile);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "noload"     => \$noload,
            "acefile:s"  => \$acefile,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species
			     );
}
my %tran2gene = $wormbase->FetchData('worm_gene2geneID_name');

$acefile = $wormbase->acefiles . "/KEGG.ace" if not defined $acefile;
my $log = Log_files->make_build_log($wormbase);
my $species_full = $wormbase->full_name;
my $fetcher = EBI_datafetch->new();

$log->write_to("\tquerying EBI web service . . \n");
my $data;
eval { 
  $data = $fetcher->get_uniprot_info_for_species($species_full);
};
$@ and $log->log_and_die("Could not successfully fetch data from EBI web-service: $@\n");
  
open(my $out_fh,">$acefile") or 
    $log->log_and_die("cant write $acefile : $!\n");
foreach my $entry_acc (keys %$data) {
  my $entry = $data->{$entry_acc};

  next if not $entry->{ec_number};
  my @wb_names = @{$entry->{wormbase_names}};
  next if not @wb_names;

  print "GOT HERE for $entry_acc\n";
  # if there are multiple names, we are assuming they are all attached to the same gene!
  my $gene;
  foreach my $cds (@wb_names) {
    if (exists $tran2gene{$cds}) {
      $gene = $tran2gene{$cds};
      last;
    }
  }
  if (not defined $gene) {
    $log->write_to("WARNING: for entry $entry_acc, could not find WB gene id for @wb_names\n");
    next;
  }

  printf $out_fh "\nGene : %s\n", $gene;
  printf $out_fh "Database NemaPath KEGG_id \"%s\"\n", $entry->{ec_number};
  printf $out_fh "Database KEGG KEGG_id \"%s\"\n", $entry->{ec_number};
}
close($out_fh) or $log->log_and_die("Could not successfully close output file $acefile\n");

unless ($noload) {
  $log->write_to("\tloading to ".$wormbase->orgdb."\n");
  $wormbase->load_to_database($wormbase->orgdb, $acefile, 'KEGG', $log);
}


$log->mail;
exit(0);
