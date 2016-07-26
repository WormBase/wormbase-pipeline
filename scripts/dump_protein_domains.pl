#!/usr/bin/env perl
# script to dump the interpro domains for each species into a flatfile in MISC_OUTPUT
# should be run post-merge

use feature qw(say);

use Ace;
use Wormbase;
use Getopt::Long;
use IO::File;
use strict;

my($test,$debug);
GetOptions(
  'test'      => \$test,
  'debug=s'   => \$debug
)||die();

my $wb = Wormbase->new(-debug => $debug, 
                       -test => $test, 
		      );

my $log = Log_files->make_build_log($wb);
my $db = Ace->connect(-path => $wb->autoace)||$log->log_and_die(Ace->error);

my %species_accessors=$wb->species_accessors;
$species_accessors{ref($wb)} = $wb;

while (my($k,$v)=each %species_accessors){

  my $outfile =$v->autoace.'/MISC_OUTPUT/protein_domains.tsv';
  my $out = IO::File->new($outfile,'w');
  $log->log_and_die("cannot open $outfile\n") unless defined $out;
  $log->write_to("dumping domains for $k to $outfile\n");

  my $cdses = $db->fetch_many(-query => "find CDS Corresponding_protein; Species=\"${\$v->long_name}\"") 
   ||$log->log_and_die(Ace->error);

  while (my $cds = $cdses->next){
    my $protein = $cds->Corresponding_protein;
    my @motifs = grep {/INTERPRO:/} $protein->Motif_homol;
    say $out join("\t",($cds->Gene,$cds->Gene->Public_name,$protein,
                map {"$_". ($_->fetch->Title?' "'.$_->fetch->Title.'"':'')} @motifs));
  }
}

$log->mail();
