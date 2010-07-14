#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  get_cgc_orthologs.pl
#
#        USAGE:  ./get_cgc_orthologs.pl
#
#      CREATED:  03/08/06 13:26:19 BST
#      
#  DESCRIPTION:  dumps all elegans orthologs including dN/dS and cgc-names from the database
#===============================================================================

use strict;
use YAML;
use IO::File;

use lib '/software/worm/ensembl/bioperl-live/';
use lib '/software/worm/ensembl/ensembl-compara/modules/';
use lib '/software/worm/ensembl/ensembl/modules/';


use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;


my %species = ( 135651 => 'Caenorhabditis brenneri',54126 => 'Pristionchus pacificus', 6239 => 'Caenorhabditis elegans',6238 => 'Caenorhabditis briggsae', 31234 => 'Caenorhabditis remanei',281687 => 'Caenorhabditis japonica',6279 => 'Brugia malayi');


my $config = ( YAML::LoadFile(glob('~/wormbase/scripts/ENSEMBL/etc/ensembl_lite.conf')) )->{'elegans'};



my %tmp1=%{&get_commondata('/nfs/wormpub/BUILD/autoace/COMMON_DATA/cds2wbgene_id.dat')};
my %wbgene2cgc=reverse %{&get_commondata('/nfs/wormpub/BUILD/autoace/COMMON_DATA/cgc_name2gene.dat',1)};
my %tmp2=%{&get_commondata('/nfs/wormpub/BUILD/brenneri/COMMON_DATA/cds2wbgene_id.dat')};
my %cds2wbgene=(%tmp1,%tmp2);

#my %cds2wbgene=%{&get_commondata('/nfs/wormpub/DATABASES/current_DB/COMMON_DATA/worm_gene2geneID_name.dat')};
#my %wbgene2cgc=reverse %{&get_commondata('/nfs/wormpub/DATABASES/current_DB/COMMON_DATA/cgc_name2gene.dat',1)};


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host   => $config->{database}->{host},
        -user   => $config->{database}->{user},
        -dbname => $config->{database}->{dbname},
        -pass   => $config->{database}->{password},
        -port   => $config->{database}->{port},
    );

my $compara_db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -host   => 'ia64d',        # change that
    -user   => 'wormro',       # and that
    -dbname => 'worm_compara'
);


my $slice_adaptor = $db->get_SliceAdaptor();
my $member_adaptor = $compara_db->get_MemberAdaptor();
my $homology_adaptor = $compara_db->get_HomologyAdaptor();

my @slices = @{$slice_adaptor->fetch_all('toplevel')};
foreach my $slice(@slices){
	foreach my $gene (@{$slice->get_all_Genes}){
		my $member = $member_adaptor->fetch_by_source_stable_id( 'ENSEMBLGENE',$gene->stable_id());
		my $homologies = $homology_adaptor->fetch_all_by_Member_method_link_type( $member, 'ENSEMBL_ORTHOLOGUES' );

		my %ids;
		foreach my $homology ( @{$homologies} ) {

			my $gid=$cds2wbgene{$gene->stable_id}?$cds2wbgene{$gene->stable_id}:$gene->stable_id;
			my $attr= join "\t", ( map { $homology->{"_$_"} } qw(dn ds n s threshold_on_ds description subtype) );

			foreach my $ma ( @{ $homology->get_all_Member_Attribute } ) {
				my ( $me, $at ) = @{$ma};
				map { 
					my $id=($cds2wbgene{$_->stable_id} || $_->stable_id);
					my $dnds=$homology->{'_dn'}/ ($homology->{'_ds'}+0.00000000000000000001);
					printf ("%.4f", $dnds) if $_->taxon_id != 6239;
					print  "\t$gid\t\"",$species{6239},"\"\t$id\t\"", $species{$_->taxon_id}, "\"\t$attr (",($wbgene2cgc{$gid}||''),'|',($wbgene2cgc{$id}||''),")\n" if $_->taxon_id != 6239
				}
				@{ $me->get_all_peptide_Members() };
			}

		}
	}	
}

sub get_commondata {
	my ($name,$unclean)=@_;
	my $file= new IO::File "< $name";
	$/=undef;
	my $data=<$file>;
	$/="\n";
	$file->close;
	my $VAR1;
	my %genehash;
	eval($data);

	return $VAR1 if $unclean;

	while(my ($k,$v)=each(%{$VAR1})){
		$k=~s/[a-z]$//;
		$genehash{$k}=$v
	}
	return \%genehash;
}
