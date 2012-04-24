#!/bin/env perl
#===============================================================================
#
#         FILE:  get_all_elegans_orthologues.pl
#
#      CREATED:  03/08/06 13:26:19 BST (mh6@sanger.ac.uk)
#===============================================================================

use strict;
use IO::File;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

use Getopt::Long;

my $comparadb;

GetOptions(
	'database=s' => \$comparadb,
)||die();


$comparadb||='worm_compara';

my $compara_db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -host   => 'farmdb1',        # change that
    -user   => 'wormro',       # and that
    -dbname => $comparadb
);


my $proteinTreeAdaptor = $compara_db->get_ProteinTreeAdaptor();

my @proteinTrees = @{$proteinTreeAdaptor->fetch_all()};

while( my $proteinTree = shift @proteinTrees){
	my @leaves = @{$proteinTree->get_all_leaves()};

	my %spec;
	map {$spec{$_->taxon_id}=1} @leaves;
	next unless (scalar(keys %spec) == scalar @leaves);
        next unless scalar @leaves >= 15;

	printf "====== %d =========\n",scalar @leaves;
	foreach my $leaf ( @leaves ) {
	      printf ">%s %s\n%s\n",$leaf->stable_id,$leaf->taxon->binomial,$leaf->translation->seq;
        }
}   


