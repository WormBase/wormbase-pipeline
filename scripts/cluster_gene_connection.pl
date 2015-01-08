#!/usr/bin/env perl

=pod

=head2 NAME - cluster_gene_connection.pl

This script makes connections between microarray expression clusters/patterns and Genes.
Requires Expression_cluster and Expr_pattern data is in database


=cut 

use lib $ENV{'CVS_DIR'};
use Wormbase;
use strict;
use Ace;
use Getopt::Long;
use Storable;
use Log_files;

my ($debug, $store, $test, $no_load, $database, $acefile);

GetOptions ( "debug:s"    => \$debug,
	     "test"       => \$test,
	     "store:s"    => \$store,
             "noload"     => \$no_load,
             "database=s" => \$database,
             "acefile=s"  => \$acefile,
	   );

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}
my $log = Log_files->make_build_log($wormbase);
$database = $wormbase->autoace if not defined $database;

$log->write_to("Connecting to database...$database\n");

my $db = Ace->connect(-path => $database) or  $log->log_and_die("Connection failure: ". Ace->error);

$log->write_to("done\n");

$acefile = $wormbase->acefiles."/cluster_gene.ace"  if not defined $acefile;
open(my $outfh, ">$acefile") or $log->log_and_die("cannot open $acefile : $!\n");

foreach my $class ("Expression_cluster", "Expr_pattern") {

  my $query="find $class WHERE Microarray_results";

  my @objs = $db->find($query);

  $log->write_to("Processing " . scalar(@objs) . " $class objects with Microarray_results...\n");

  if(scalar @objs < 10 ) {
    $log->log_and_die("There are less than 10 $class objects with Microarray_results in the database which is probably not right!");
  }

  foreach my $obj (@objs) {
    
    print $outfh "\n$class : \"$obj\"\n";
    print $outfh "-D Gene\n";
    print $outfh "\n$class : \"$obj\"\n";
    
    my @mrs = $obj->Microarray_results(1);
    foreach my $mr (@mrs) {
      my @genes = $mr->Gene;
      foreach my $gene (@genes) {
        print $outfh "Gene\t\"$gene\"\n";
      }
    }
  }
  
  #
  # For Expression_clusters without Microarray_results, many will have
  # direct Gene connections from citace. We want to mostly retain these, 
  # but remove the Dead genes
  #
  $query = "FIND $class WHERE Gene AND NOT Microarray_results";
  @objs = $db->find($query);

  $log->write_to("Processing " . scalar(@objs) . " $class objects without Microarray_results but with direct gene connections...\n");
  foreach my $obj (@objs) {
    my (@reject_genes);

    my @genes = $obj->Gene;
    foreach my $g (@genes) {
      if ($g->Status ne 'Live') {
        push @reject_genes, $g;
      } 
    }

    foreach my $rg (@reject_genes) {
      print $outfh "\n$class : \"$obj\"\n";
      print $outfh "-D Gene \"$rg\"\n";
    }

    if (@reject_genes) {
      $log->write_to("Removed from $obj because no longer live: @reject_genes\n");
    }
  }

  $log->write_to("Processed objects $class objects without Microarray_results but direct Gene connections\n");
}

close($outfh) or $log->log_and_die("Could not cleanly close output file\n");
$db->close;

$wormbase->load_to_database($database, $acefile, 'cluster_genes', $log) unless $no_load;

$log->mail;
exit(0);
