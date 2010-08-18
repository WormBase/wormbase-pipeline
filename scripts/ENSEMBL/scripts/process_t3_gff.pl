#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  process_t3_gff.pl
#
#        USAGE:  cat gff3_file | ./process_t3_gff.pl "Flunky speciesname"
#
#      CREATED:  03/08/06 13:26:19 BST
#      
#  DESCRIPTION:  pulls from mysql the orthologs and from the build the gene/cds 
#                connections to add an Alias to the GFF3 gene name 
#===============================================================================

use strict;
use YAML;
use IO::File;

use lib '/software/worm/ensembl/bioperl-live/';
use lib '/software/worm/ensembl/ensembl-compara/modules/';
use lib '/software/worm/ensembl/ensembl/modules/';

my $le_species=shift;
#$le_species=~/^(\w)\w+\s(\w\w\w)/;
$le_species=~/^(\w)\w+\s(\w\w)/;

my $prefix="$1$2";

print STDERR "\x1b[38;5;31m","slurping gene names from ensembl for $le_species / $prefix\n","\x1b[0m";
my %g2cgc=%{grab_ids($le_species)};

print STDERR  "\x1b[38;5;31m";
printf STDERR "slurped %d genes\n",scalar(keys %g2cgc);
print STDERR "\x1b[0m";

while (<>){
	if(/WormBase\s+gene/){
	  chomp;
	  print $_;
	  /Name=([^;]+);/;
	  print ";Alias=$prefix-$g2cgc{$1}" if $g2cgc{$1};
	  print "\n";
	}
	else {
		print
	}
}

sub grab_ids {
	my ($speci)=@_;

        my %gid2cgc;

        use Bio::EnsEMBL::DBSQL::DBAdaptor;
        use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;


        my %species = ( 135651 => 'Caenorhabditis brenneri',
	        54126  => 'Pristionchus pacificus', 
		6239   => 'Caenorhabditis elegans',
		6238   => 'Caenorhabditis briggsae',
		31234  => 'Caenorhabditis remanei',
		281687 => 'Caenorhabditis japonica',
		6279   => 'Brugia malayi',
		54126  => 'Pristionchus pacificus',
                6305   => 'Meloidogyne hapla',
                6289   => 'Haemonchus contortus',
		96668  => 'Caenorhabditis angaria',
	);


       my $config = ( YAML::LoadFile(glob('~/wormbase/scripts/ENSEMBL/etc/ensembl_lite.conf')) )->{'elegans'};

       my %cds2wbgene=%{&get_commondata('/nfs/wormpub/DATABASES/current_DB/COMMON_DATA/cds2wbgene_id.dat')};
       my %wbgene2cgc=reverse %{&get_commondata('/nfs/wormpub/DATABASES/current_DB/COMMON_DATA/cgc_name2gene.dat',1)};

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
			foreach my $ma ( @{ $homology->get_all_Member_Attribute } ) {
				my ( $me, $at ) = @{$ma};
				map { 
#				    my $id=($cds2wbgene{$_->stable_id} || $_->stable_id);
                                    my $id=($cds2wbgene{$_->stable_id} ||  $_->get_Gene->stable_id);
				    if (($homology->{_description} eq 'ortholog_one2one') && ($species{$_->taxon_id} eq $speci) && $wbgene2cgc{$gid} && $_->taxon_id != 6239){
					    $gid2cgc{$id}=$wbgene2cgc{$gid};
					    print STDERR "$id -> $wbgene2cgc{$gid}\n";
				    }
				}
				@{ $me->get_all_peptide_Members() };
			}
		}
	}	
   }
   return \%gid2cgc
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
