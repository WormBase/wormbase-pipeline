#!/usr/bin/perl -w

=pod

=head2 NAME - cluster_gene_connection.pl

This script makes connections between microarray expression clusters and Genes (original by Igor). Requires Expression_cluster data is in database

script_template.pl MANDATORY arguments:

=over 4

=item None

=back

Optional args 

=over 4

=item debug , test, store

=back

=cut 

use lib $ENV{'CVS_DIR'};
use Wormbase;
use strict;
use Ace;
use Getopt::Long;
use Storable;
use Log_files;

my ($debug, $store, $test);
GetOptions ( "debug:s" => \$debug,
	     "test"    => \$test,
	     "store:s" => \$store
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

my $query="find expression_cluster";

my @clusters=$db->find($query);

$log->write_to(scalar @clusters. " objects retrieved.\n");
if(scalar @clusters < 10 ) {
  $log->log_and_die("There are less than 10 clusters in the database which is probably not right!");
}

my $i=0;
foreach my $cluster (@clusters) {
    print OUT "Expression_cluster : \"$cluster\"\n";
    print OUT "-D Gene\n";
    print OUT "\n";
    $i++;
    if ($i % 10 == 0) {
      $log->write_to("$i objects processed.\n");
    }
    print OUT "Expression_cluster : \"$cluster\"\n";
    my @mrs=$cluster->Microarray_results(1);
    foreach my $mr (@mrs) {
	my @genes=$mr->Gene;
	foreach my $gene (@genes) {
	    print OUT "Gene\t\"$gene\"\n";
	}
    }
    print OUT "\n";
}

$log->write_to("$i objects processed.\n");

$log->mail;
exit(0);
