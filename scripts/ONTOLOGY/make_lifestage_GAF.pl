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

$outfile ||= $wormbase->ontology . "/lifestage_association." . $wormbase->get_wormbase_version_name . ".wb";
open(my $outfh, ">$outfile" ) or $log->log_and_die("cannot open $outfile : $!\n");

&print_wormbase_GAF_header($outfh);

$it = $db->fetch_many( -query => 'find Expr_pattern Life_stage' );
while ( my $obj = $it->next ) {
  $count++;
  $log->write_to("$count objects processed\n") if $verbose and ( $count % 1000 == 0 );
  
  my ( %genes, %ls, %auth);
 
  foreach my $g ($obj->Gene) {
    next if not exists $gene_info->{$g} or $gene_info->{status} eq 'Dead';
    $genes{$g->name}++;
  }
  foreach my $ls ($obj->Life_stage) {
    my $qual = $ls->right;
    $ls{$ls->name} = ($qual) ? $qual : "";
  }

  foreach my $g ( keys %genes ) {
    foreach my $ls ( keys %ls ) {
      my $qual = $ls{$ls};
      &print_wormbase_GAF_line($outfh, 
                               $g, 
                               $gene_info->{$g}->{public_name}, 
                               $qual, 
                               $ls, 
                               "WB_EP:" . $obj->name,
                               "IDA", 
                               "",
                               "L",  
                               $gene_info->{$g}->{sequence_name},
                               $taxid, 
                               $date);
    }
  }
}

$it = $db->fetch_many(-query => "Find Gene RNAseq_FPKM AND Species = \"$full_name\" AND NOT Dead");
while (my $g = $it->next) {
  foreach my $ls ($g->RNASeq_FPKM) {
    my ($sra_acc, $max_fpkm);
    foreach my $data ($ls->col) {
      if ($data->name >= 10.0) {  # semi-arbitrary cut-off for confirmed expression
        if ($data->right and $data->right->name eq 'From_analysis') {
          my $ana = $data->right->right->name;
          if ($ana =~ /\.(SRX\S+)/) {
            if (not defined $sra_acc or $data->name > $max_fpkm) {
              $sra_acc = $1;
              $max_fpkm = $data->name;
            }
          }
        }
      }
    }
    if (defined $sra_acc) {
      &print_wormbase_GAF_line($outfh, 
                               $g, 
                               $gene_info->{$g}->{public_name}, 
                               "",
                               $ls, 
                               "ENA:$sra_acc",
                               "IEA", 
                               "",
                               "L",  
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
