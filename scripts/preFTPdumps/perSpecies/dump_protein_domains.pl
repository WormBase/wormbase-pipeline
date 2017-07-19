#!/usr/bin/env perl
#
# script to dump the interpro domains for each species into a flatfile in REPORTs

use strict;
use Getopt::Long;
use IO::File;
use Storable;

use lib $ENV{CVS_DIR};

use Ace;
use Wormbase;
use Log_files;

my($test,$debug, $species, $store, $wb, $database, $outfile);

GetOptions(
  'test'       => \$test,
  'debug=s'    => \$debug,
  'species=s'  => \$species,
  'store=s'    => \$store,
  'outfile=s'  => \$outfile,
  'database=s' => \$database,
)||die();

my $wormbase;
if ($store) { 
  $wb = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")
}else {
  $wb = Wormbase->new( -debug => $debug, 
                       -test => $test,
                       -organism => $species);
}

$species = $wb->species;
my $full = $wb->long_name;
my $log = Log_files->make_build_log($wb);

$database = $wb->autoace if not defined $database;

$outfile = $wb->reports . "/protein_domains.csv" if not defined $outfile;
my $out = IO::File->new($outfile,'w');

$log->log_and_die("cannot open $outfile\n") unless defined $out;
$log->write_to("dumping domains for $species to $outfile\n");

my $db = Ace->connect(-path => $database) or $log->log_and_die("Could not connect to $database\n");

my $cdses = $db->fetch_many(-query => "find CDS Corresponding_protein; Species=\"$full\"") 
    ||$log->log_and_die(Ace->error);

while (my $cds = $cdses->next){
  my $protein = $cds->Corresponding_protein;
  my @motifs = grep {/INTERPRO:/} $protein->Motif_homol;
  print $out join("\t",($cds->Gene,$cds->Gene->Public_name,$protein,
                      map {"$_". ($_->fetch->Title?' "'.$_->fetch->Title.'"':'')} @motifs)), "\n";
}
$db->close();

$log->mail();
exit(0);
