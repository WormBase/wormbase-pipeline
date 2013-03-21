#!/usr/bin/env perl

=pod

=head2 NAME - cluster_gene_connection.pl

This script makes connections between microarray expression clusters and Genes (original by Igor). Requires Expression_cluster data is in database


=cut 

use lib $ENV{'CVS_DIR'};
use Wormbase;
use strict;
use Ace;
use Getopt::Long;
use Storable;
use Log_files;

my ($debug, $store, $test,$no_load);
GetOptions ( "debug:s" => \$debug,
	     "test"    => \$test,
	     "store:s" => \$store,
             "noload"  => \$no_load,
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
my $database = $wormbase->autoace;

$log->write_to("Connecting to database...$database\n");

my $db = Ace->connect(-path => $database) or  $log->log_and_die("Connection failure: ". Ace->error);

$log->write_to("done\n");

my $acefile = $wormbase->acefiles."/cluster_gene.ace";
open OUT, ">$acefile" or $log->log_and_die("cannot open $acefile : $!\n");

my $query="find expression_cluster WHERE Microarray_results";

my @clusters=$db->find($query);

$log->write_to("Processing " .scalar @clusters. " objects with Microarray_results...\n");
if(scalar @clusters < 10 ) {
  $log->log_and_die("There are less than 10 clusters in the database which is probably not right!");
}

my $i=0;
foreach my $cluster (@clusters) {
  print OUT "\nExpression_cluster : \"$cluster\"\n";
  print OUT "-D Gene\n";
  print OUT "\nExpression_cluster : \"$cluster\"\n";

  my @mrs=$cluster->Microarray_results(1);
  foreach my $mr (@mrs) {
    my @genes=$mr->Gene;
    foreach my $gene (@genes) {
      print OUT "Gene\t\"$gene\"\n";
    }
  }
}
$log->write_to("Processed objects with Microarray_results\n");

#
# For Expression_clusters without Microarray_results, many will have
# direct Gene connections from citace. We want to mostly retain these, 
# but remove the Dead genes
#
$query = "FIND Expression_cluster WHERE Gene AND NOT Microarray_results";
@clusters = $db->find($query);

$log->write_to("Processing " .scalar @clusters. " objects without Microarray_results but with direct gene connections...\n");
foreach my $cluster (@clusters) {
  my (@keep_genes, @reject_genes);

  my @genes = $cluster->Gene;
  foreach my $g (@genes) {
    if ($g->Status eq 'Live') {
      push @keep_genes, $g;
    } else {
      push @reject_genes,$g;
    }
  }

  print OUT "\nExpression_cluster : \"$cluster\"\n";
  print OUT "-D Gene\n";
  print OUT "\nExpression_cluster : \"$cluster\"\n";
  foreach my $g (@keep_genes) {
    print OUT "Gene \"$g\"\n";
  }

  if (@reject_genes) {
    $log->write_to("Removed from $cluster because no longer live: @reject_genes\n");
  }
}
close(OUT);
$db->close;

$log->write_to("Processed objects without Microarray_results but direct Gene connections\n");

$wormbase->load_to_database($database, $acefile, 'cluster_genes', $log) unless $no_load;



$log->mail;
exit(0);
