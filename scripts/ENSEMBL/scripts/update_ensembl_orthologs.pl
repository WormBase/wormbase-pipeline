#!/usr/bin/perl -w
#===============================================================================
#   updates ortholog connections from the main Compara database
#
#       AUTHOR:   (), <>
#      VERSION:  1.0
#      CREATED:  03/08/06 13:26:19 BST
#     REVISION:  ---
#===============================================================================

use strict;
use IO::File;
use lib '/software/worm/ensembl/ensembl/modules';
use lib '/software/worm/ensembl/ensembl-compara/modules';
use lib '/software/worm/ensembl/bioperl-live';

use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ens-livemirror', -user => 'wormro',-port => 3306,-verbose => 0,-db_version => 50);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;


my %cds2wbgene=%{&get_commondata('/nfs/disk100/wormpub/DATABASES/current_DB/COMMON_DATA/cds2wbgene_id.dat')};

my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Caenorhabditis elegans','core','Slice');
my $member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','Member');
my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','Homology');;

my @slices = @{$slice_adaptor->fetch_all('toplevel')};


foreach my $slice(@slices){
	my $genes = $slice->get_all_Genes;
	while (my $gene = shift @{$genes}){

		my $gid=$cds2wbgene{$gene->stable_id}?$cds2wbgene{$gene->stable_id}:next; # don't import outdated genes

		my $member = $member_adaptor->fetch_by_source_stable_id( 'ENSEMBLGENE',$gene->stable_id());
		next unless $member;
		my $homologies = $homology_adaptor->fetch_all_by_Member_method_link_type( $member, 'ENSEMBL_ORTHOLOGUES' );

		my %homol_ids;
		my %omims;
		my %omim_genes;

		foreach my $homology ( @{$homologies} ) {
			foreach my $ma ( @{ $homology->get_all_Member_Attribute } ) {
				my ( $me, $at ) = @{$ma};
				foreach my $pepm ( @{ $me->get_all_peptide_Members() } ) { 
					if ($pepm->taxon_id != 6239){
						$homol_ids{ $pepm->stable_id } = $pepm; #$homology->description,$homology->subtype];
					}

				}
			}
		}

		next unless %homol_ids;


		print "Gene : \"$gid\"\n";
		while (my ($k,$v)=each(%homol_ids)){
#				print "Ortholog_other EnsEMBL gene $k \"${\$$v[0]->name}\"  From_analysis EnsEMBL-Compara\n";
				print "Ortholog_other ENSEMBL:$k  From_analysis EnsEMBL-Compara\n";
		}
		print "\n";

		while (my ($k,$v)=each(%homol_ids)){
					# put the OMIM orthologs into another hash based on the MIM_GENE xref
					print "Protein : ENSEMBL:$k\n";
					print "Species \"${\$v->taxon->name}\"\n";
					print "DB_info DatabaseEnsEMBL ENSEMBL_proteinID $k\n";
					if ($v->taxon_id == 9606){ # meaning if human
					        # uses an undocumented function of get_all_DBEntries, so lets hope it stays
						map {printf "DB_info Database OMIM gene %s\n",$_->primary_id} @{$v->gene->get_all_DBLinks('MIM_MORBID')};
						map {printf "DB_info Database OMIM disease %s\n",$_->primary_id} @{$v->gene->get_all_DBLinks('MIM_GENE')};
					}
					print "\n";
		}
	}	
}


sub get_commondata {
	my ($name)=@_;
	my $file= new IO::File "< $name";
	$/=undef;
	my $data=<$file>;
	$/="\n";
	$file->close;
	my $VAR1;
	my %genehash;
	eval($data);

	while(my ($k,$v)=each(%{$VAR1})){
		$k=~s/[a-z]$//;
		$genehash{$k}=$v
	}
	return \%genehash;
}
