#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Storable;
use Ace;

use lib $ENV{CVS_DIR};
use Wormbase;
use Log_files;

use lib "$ENV{CVS_DIR}/ONTOLOGY";
use GAF;


my ( $help, $debug, $test, $store, $wormbase );
my ( $outfile, $acedbpath, $verbose );

GetOptions(
  "help"       => \$help,
  "debug=s"    => \$debug,
  "test"       => \$test,
  "store:s"    => \$store,
  "database:s" => \$acedbpath,
  "output:s"   => \$outfile,
  "verbose"    => \$verbose,
    );

if ($store) {
  $wormbase = retrieve($store)
      or croak("Can't restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new(
    -debug => $debug,
    -test  => $test,
      );
}

my $log = Log_files->make_build_log($wormbase);
my $date = &get_GAF_date();
my $taxid = $wormbase->ncbi_tax_id;
my $full_name = $wormbase->full_name;

$acedbpath ||= $wormbase->autoace;
my $tace = $wormbase->tace;

$log->write_to("connecting to database... $acedbpath\n");

my $db = Ace->connect( -path => $acedbpath, -program => $tace )
  or $log->log_and_die( "Connection failure: " . Ace->error );

my ($gene_info, $it, $count);

$gene_info = &get_gene_info( $acedbpath, $wormbase, $full_name );
$log->write_to( scalar(keys %$gene_info) . " genes read\n" ) if $verbose;

$outfile ||= $wormbase->ontology . "/anatomy_association." . $wormbase->get_wormbase_version_name . ".wb";
open(my $outfh, ">$outfile" ) or $log->log_and_die("cannot open $outfile : $!\n");

$it = $db->fetch_many( -query => 'find Expr_pattern Anatomy_term' );
while ( my $obj = $it->next ) {
  $count++;
  $log->write_to("$count objects processed\n") if $verbose and ( $count % 1000 == 0 );
  
  my ( %genes, %at, %auth);
 
  foreach my $g ($obj->Gene) {
    next if not exists $gene_info->{$g} or $gene_info->{status} eq 'Dead';
    $genes{$g->name}++;
  }
  foreach my $at ($obj->Anatomy_term) {
    $at{$at} = "";
    foreach my $tag ($at->col) {
      next if $tag->name =~ /Life_stage/;
      $at{$at} = $tag->name;
    }
  }

  # take first reference only
  my ($ref) = $obj->Reference;
  
  foreach my $g ( keys %genes ) {
    foreach my $a ( keys %at ) {
      my $qual = $at{$a};
      &print_wormbase_GAF_line($outfh, 
                               $g, 
                               $gene_info->{$g}->{public_name}, 
                               $qual, 
                               $a, 
                               "WB_REF:$ref",
                               "IEP", 
                               "WB:" . $obj->name, 
                               "A",  
                               $gene_info->{$g}->{sequence_name},
                               $taxid, 
                               $date);
    }
  }
}

close($outfh);
$db->close;
$log->mail;

exit(0);
