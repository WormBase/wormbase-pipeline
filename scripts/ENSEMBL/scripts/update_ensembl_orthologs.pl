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
use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ensembldb.ensembl.org', -user => 'anonymous',-port => 3306,-verbose => 1,-db_bversion => 47);

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;


my %cds2wbgene=%{&get_commondata('/nfs/disk100/wormpub/DATABASES/current_DB/COMMON_DATA/cds2wbgene_id.dat')};

my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Caenorhabditis elegans','core','Slice');
my $member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','Member');
my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','Homology');;

my @slices = @{$slice_adaptor->fetch_all('toplevel')};

foreach my $slice(@slices){
	foreach my $gene (@{$slice->get_all_Genes}){

		my $gid=$cds2wbgene{$gene->stable_id}?$cds2wbgene{$gene->stable_id}:next; # don't import outdated genes

		my $member = $member_adaptor->fetch_by_source_stable_id( 'ENSEMBLGENE',$gene->stable_id());
		next unless $member;
		my $homologies = $homology_adaptor->fetch_all_by_Member_method_link_type( $member, 'ENSEMBL_ORTHOLOGUES' );

		my %homol_ids;
		my %omims;

		foreach my $homology ( @{$homologies} ) {
			foreach my $ma ( @{ $homology->get_all_Member_Attribute } ) {
				my ( $me, $at ) = @{$ma};
				foreach my $pepm ( @{ $me->get_all_peptide_Members() } ) { 
					if ($pepm->taxon_id != 6239){
						$homol_ids{ $pepm->gene->stable_id } = [$pepm->taxon,$homology->description,$homology->subtype];
					}

					# put the OMIM orthologs into another hash based on the MIM_GENES xref
					if ($pepm->taxon_id == 9606){
					        # uses an undocumented function of get_all_DBEntries, so lets hope it stays
						$omims{$pepm->gene->stable_id}=1 if $pepm->gene->get_all_DBEntries('MIM_GENES');
					}
				}
			}
		}

		next unless %homol_ids;


		print "Gene : \"$gid\"\n";
		while (my ($k,$v)=each(%homol_ids)){
				print "Ortholog_other EnsEMBL gene $k \"${\$$v[0]->name}\" Inferred_automatically compara\n";
		}
		foreach my $key(keys %omims){
			print "Ortholog_other EnsEMBL gene $key \"Homo sapiens\" Inferred_automatically OMIM_compara\n";
		}

		print "\n";

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
