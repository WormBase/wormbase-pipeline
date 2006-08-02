#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  check_genes.pl
#
#        USAGE:  ./check_genes.pl 
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:   (), <>
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  11/07/06 17:49:26 BST
#     REVISION:  ---
#===============================================================================

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Long;
use YAML;
use Test::More qw(no_plan);

my ($species);
GetOptions( 'species=s'=>\$species,);

my $config = ( YAML::LoadFile('ensembl_lite.conf') )->{$species};

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host   => $config->{database}->{host},
        -user   => $config->{database}->{user},
        -dbname => $config->{database}->{dbname},
        -pass   => $config->{database}->{password},
        -port   => $config->{database}->{port},
    );



my $slice_adaptor = $db->get_SliceAdaptor();
my @slices = @{$slice_adaptor->fetch_all('toplevel')};
foreach my $slice(@slices){
	my $slname=$slice->seq_region_name();
	my $codontable=($slname eq 'CHROMOSOME_MtDNA')?5:1; #elegans hack
	my $table  = Bio::Tools::CodonTable -> new ( -id => $codontable );
	ok($slice,"=> Superlink $slname");
	foreach my $gene (@{$slice->get_all_Genes}){
		ok($gene,"Gene ".$gene->stable_id());
		foreach my $trans (@{$gene->get_all_Transcripts()}) {
			my $protein=$trans->translate()->seq;
			ok($trans,'- Transcript '.$trans->stable_id());
		       	ok(length($protein)>0,'-- minimum one AA '.$trans->stable_id());
		       	ok($table->is_start_codon(substr($trans->translateable_seq,0,3)),'-- translation start of '.$trans->stable_id());
			ok($table->is_ter_codon(substr($trans->translateable_seq,-3,3)),'-- translation stop of '.$trans->stable_id());
			ok(($protein=~tr/\*/\*/)==0,'-- internal stop codons within CDS of '.$trans->stable_id());
		}
	}
}


