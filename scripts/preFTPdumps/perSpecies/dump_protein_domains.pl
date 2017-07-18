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

my($test,$debug, $species, $store, $wb, $outfile);

GetOptions(
  'test'      => \$test,
  'debug=s'   => \$debug,
  'species=s' => \$species,
  'store=s'   => \$store,
  'outfile=s' => \$outfile,
)||die();

my $wormbase;
if ($store) { 
  $wb = Storable::retrieve($store) or croak("Can't restore wormbase from $store\n")
}else {
  $wb = Wormbase->new( -debug => $debug, -test => $test,-organism => $species)
}

if (not defined $outfile) {
  my $fname = "protein_domains.csv";
  $outfile = $wb->reports."/$fname";
}

$species = $wb->species;

my $full = $wb->long_name;
my $log = Log_files->make_build_log($wb);
my $out = IO::File->new($outfile,'w');

$log->log_and_die("cannot open $outfile\n") unless defined $out;
$log->write_to("dumping domains for $species to $outfile\n");

my $db = Ace->connect(-path => $wb->autoace)||$log->log_and_die(Ace->error);

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
