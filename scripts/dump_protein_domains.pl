#!/usr/bin/env perl
# dump the protein domains into a flatfile

use feature qw(say);

use Ace;
use Wormbase;
use Getopt::Long;
use IO::File;
use strict;


my($species,$store,$test,$debug);
GetOptions(
  'species=s' => \$species,
  'store=s'   => \$store,
  'test'      => \$test,
  'debug=s'   => \$debug
)||die();

my $wb;

if ($store) {
  $wb = Storable::retrieve($store)
      or croak("cannot restore wormbase from $store");
}
else { 
  $wb = Wormbase->new( -debug => $debug, 
                       -test => $test, 
                       -organism => $species,
		     )
}

my $log = Log_files->make_build_log($wb);

my $db = Ace->connect(-path => $wb->autoace)||$log->log_and_die(Ace->error);
my $outfile =$wb->autoace.'/MISC_OUTPUT/protein_domains.tsv';
my $out = IO::File->new($outfile,'w');
$log->log_and_die("cannot open $outfile\n") unless defined $out;

my $cdses = $db->fetch_many(-query => "find CDS Corresponding_protein; Species=\"${\$wb->long_name}\"") 
   ||$log->log_and_die(Ace->error);

while (my $cds = $cdses->next){
  my $protein = $cds->Corresponding_protein;
  my @motifs = grep {/INTERPRO:/} $protein->Motif_homol;
  say $out join("\t",($cds->Gene,$cds->Gene->Public_name,$protein,
                map {"$_". ($_->fetch->Title?' "'.$_->fetch->Title.'"':'')} @motifs));
}

$log->mail();
