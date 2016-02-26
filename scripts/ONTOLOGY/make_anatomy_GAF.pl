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
my ( %output_hash, $outfile, $acedbpath, $verbose );

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
      or croak("Cannot restore wormbase from $store\n");
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
my $tace = $wormbase->tace;

$acedbpath ||= $wormbase->autoace;
$outfile ||= $wormbase->ontology . "/anatomy_association." . $wormbase->get_wormbase_version_name . ".wb";

$log->write_to("connecting to database... $acedbpath\n");

my $db = Ace->connect( -path => $acedbpath, -program => $tace )
  or $log->log_and_die( "Connection failure: " . Ace->error );

my ($gene_info, $it, $count);

$gene_info = &get_gene_info( $acedbpath, $wormbase, $full_name );
$log->write_to( scalar(keys %$gene_info) . " genes read\n" ) if $verbose;

$log->write_to("Querying Expr_pattern objects...\n") if $verbose;

$it = $db->fetch_many( -query => 'find Expr_pattern Anatomy_term' );
while ( my $obj = $it->next ) {
  $count++;

  
  my ( %genes, %at );
 
  foreach my $g ($obj->Gene) {
    next if not exists $gene_info->{$g} or $gene_info->{status} eq 'Dead';
    $genes{$g->name}++;
  }
  foreach my $at ($obj->Anatomy_term) {
    $at{$at} = {};
    foreach my $tag ($at->col) {
      # Note from Caltech Oct 2015:

      # "We would like to remove 'Life_stage' from column 4 of the
      # anatomy_association GAF -anatomy_association.WSxxx.wb- as it is not
      # informative.
      # That column contains qualifiers from the Anatomy_term Qualifier hash, such
      # as Certain, Uncertain, Partial

      # Anatomy_term ?Anatomy_term XREF Expr_pattern #Qualifier

      # When we started to capture life_stage specific expression we added
      # Life_stage in the qualifier hash to maintain the relationship
      # anatomy-developmental stage. But in the GAF the tag Life stage per se can
      # be confusing."
      next if $tag->name =~ /Life_stage/;
      $at{$at}->{$tag->name} = 1;
    }
    # if there were no qualifiers, add in an empty one
    $at{$at}->{""} = 1 if not keys %{$at{$at}};
  }

  # take first reference only
  my ($ref) = $obj->Reference;
  
  foreach my $g ( keys %genes ) {
    foreach my $a ( keys %at ) {
      foreach my $qual (keys %{$at{$a}}) {
        $output_hash{$g}->{$a}->{$ref}->{$qual}->{$obj->name} = 1;
      }
    }
  }
}

$log->write_to("Querying Expression_cluster objects...\n") if $verbose;

$it = $db->fetch_many( -query => 'find Expression_cluster Anatomy_term' );
while ( my $obj = $it->next ) {
  my ( %genes, %at );
 
  foreach my $g ($obj->Gene) {
    next if not exists $gene_info->{$g} or $gene_info->{status} eq 'Dead';
    $genes{$g->name}++;
  }
  foreach my $at ($obj->Anatomy_term) {
    foreach my $tag ($at->col) {
      if ($tag->name eq 'Enriched') {
        $at{$at} = $tag->name;
      }
    }
  }

  # again, take first reference only
  my ($ref) = $obj->Reference;
  
  foreach my $g ( keys %genes ) {
    foreach my $a ( keys %at ) {
      my $qual = $at{$a};
      $output_hash{$g}->{$a}->{$ref}->{$qual}->{$obj->name} = 1;
    }
  }
}

open(my $outfh, ">$outfile" ) or $log->log_and_die("cannot open $outfile : $!\n");
&print_wormbase_GAF_header($outfh);

foreach my $g (sort keys %output_hash) {
  foreach my $at (sort keys %{$output_hash{$g}}) {
    foreach my $ref (sort keys %{$output_hash{$g}->{$at}}) {
      foreach my $qual (sort keys %{$output_hash{$g}->{$at}->{$ref}}) {
        my @ec_objs = sort keys  %{$output_hash{$g}->{$at}->{$ref}->{$qual}};
        
        &print_wormbase_GAF_line($outfh, 
                                 $g, 
                                 $gene_info->{$g}->{public_name}, 
                                 $qual, 
                                 $at, 
                                 "WB_REF:$ref",
                                 "IDA", 
                                 join("|", map { "WB:$_" } @ec_objs), 
                                 "A",  
                                 $gene_info->{$g}->{sequence_name},
                                 $taxid, 
                                 $date);
      }
    }
  }
}


close($outfh);
$db->close;
$log->mail;

exit(0);
