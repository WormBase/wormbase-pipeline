#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  get_all_elegans_orthologues.pl
#
#        USAGE:  ./get_all_elegans_orthologues.pl 
#
#       AUTHOR:   (), <>
#      VERSION:  1.0
#      CREATED:  03/08/06 13:26:19 BST
#     REVISION:  ---
#===============================================================================

use strict;
use YAML;
use FindBin;
use IO::File;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;


my %species = ( 
	6239   => 'Caenorhabditis elegans',
	6238   => 'Caenorhabditis briggsae', 
	31234  => 'Caenorhabditis remanei',
#	135651 => 'Caenorhabditis brenneri',
	6279   => 'Brugia malayi'
);

my $config = ( YAML::LoadFile("$FindBin::Bin/../etc/ensembl_lite.conf") )->{'elegans'};
my %cds2wbgene=%{&get_commondata('cds2wbgene_id')};

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
        -host   => $config->{database}->{host},
        -user   => $config->{database}->{user},
        -dbname => $config->{database}->{dbname},
        -pass   => $config->{database}->{password},
        -port   => $config->{database}->{port},
    );

my $compara_db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -host   => 'ia64b',        # change that
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

		# needs some better way
		my %briggsae_ids;
		my %remanei_ids;
		my %brenneri_ids;
		my %brugia_ids;
		my %elegans_ids;


		foreach my $homology ( @{$homologies} ) {
			foreach my $ma ( @{ $homology->get_all_Member_Attribute } ) {
				my ( $me, $at ) = @{$ma};
				foreach my $pepm ( @{ $me->get_all_peptide_Members() } ) { 
					if ($pepm->taxon_id==6238){
						$briggsae_ids{ $pepm->stable_id } = [$pepm->taxon_id,$homology->description,$homology->subtype] 
					}
					elsif ($pepm->taxon_id==31234) {
						$remanei_ids{ $pepm->stable_id } = [$pepm->taxon_id,$homology->description,$homology->subtype]
					} 
					elsif ($pepm->taxon_id==135651){
						$brenneri_ids{ $pepm->stable_id } = [$pepm->taxon_id,$homology->description,$homology->subtype]
					}
					elsif ($pepm->taxon_id==6279){
                                                $brugia_ids{ $pepm->stable_id } = [$pepm->taxon_id,$homology->description,$homology->subtype]
					}
                        		elsif ($pepm->taxon_id==6239){
						$elegans_ids{ $pepm->stable_id } = [$pepm->taxon_id,$homology->description,$homology->subtype] 
					}

					else {print STDERR "cannot find:".$pepm->taxon_id."\n"}

				}

			}
		}

		my $gid=$cds2wbgene{$gene->stable_id}?$cds2wbgene{$gene->stable_id}:$gene->stable_id;

		next unless %briggsae_ids || %remanei_ids;


		# the next bit is especially ugly. It is possible to make it nicer by moving the species loops into the elegans,
		# but loading speed would suffer.

		# elegans part

		# briggsae
		print "Gene : \"$gid\"\n";
		while (my ($k,$v)=each(%briggsae_ids)){
			        my $sid=$cds2wbgene{$k}?$cds2wbgene{$k}:$k;
				print "Ortholog $sid \"$species{$$v[0]}\" From_analysis WormBase-Compara\n";
		}

		# remanei
		while (my ($k,$v)=each(%remanei_ids)){
			        my $sid=$cds2wbgene{$k}?$cds2wbgene{$k}:$k;
				print "Ortholog $sid \"$species{$$v[0]}\" From_analysis WormBase-Compara\n";
		}
		# brugia
                while (my ($k,$v)=each(%brugia_ids)){
		                my $rid=$cds2wbgene{$k}?$cds2wbgene{$k}:$k;
				print "Ortholog_other Brugia_database gene $rid \"$species{$$v[0]}\" From_analysis WormBase-Compara\n"
		}

		print "\n";

		# briggsae part

		# briggsae	
		while (my ($k,$v)=each(%briggsae_ids)){
			   	my $sid=$cds2wbgene{$k}?$cds2wbgene{$k}:$k;
				print "Gene : \"$sid\"\n";
				print "Ortholog $gid \"${\$species{ $config->{taxon_id}}}\" From_analysis WormBase-Compara\n"; # elegans
				while (my ($r_k,$r_v)=each(%remanei_ids)){                                                     # remanei
				        my $rid=$cds2wbgene{$r_k}?$cds2wbgene{$r_k}:$r_k;
					print "Ortholog $rid \"$species{$$r_v[0]}\" From_analysis WormBase-Compara\n"
				}
                                while (my ($r_k,$r_v)=each(%brugia_ids)){                                                      # brenneri
				        my $rid=$cds2wbgene{$r_k}?$cds2wbgene{$r_k}:$r_k;
					print "Ortholog_other Brugia_database gene $rid \"$species{$$r_v[0]}\" From_analysis WormBase-Compara\n"
				}

				print "\n";
		
		}

		# remanei part
		while (my ($k,$v)=each(%remanei_ids)){
			   	my $sid=$cds2wbgene{$k}?$cds2wbgene{$k}:$k;
				print "Gene : \"$sid\"\n";
				print "Ortholog $gid \"${\$species{ $config->{taxon_id}}}\" From_analysis WormBase-Compara\n"; # elegans
				while (my ($r_k,$r_v)=each(%briggsae_ids)){                                                    # briggsae
				        my $rid=$cds2wbgene{$r_k}?$cds2wbgene{$r_k}:$r_k;
					print "Ortholog $rid \"$species{$$r_v[0]}\" From_analysis WormBase-Compara\n"
				}
                                while (my ($r_k,$r_v)=each(%brugia_ids)){                                                      # brigia
				        my $rid=$cds2wbgene{$r_k}?$cds2wbgene{$r_k}:$r_k;
					print "Ortholog_other Brugia_databse gene $rid \"$species{$$r_v[0]}\" From_analysis WormBase-Compara\n"
				}

				print "\n";
		
		}

	}	
}



# needs a merging step
# /nfs/disk100/wormpub/BUILD/autoace/COMMON_DATA/
# /nfs/disk100/wormpub/BUILD/remanei/COMMON_DATA/
sub get_commondata {
	my ($name)=@_;
	my %genehash;
	my @locations=qw(autoace remanei briggsae);
	my $dir=glob('~wormpub/BUILD/');
	foreach my $loc(@locations) {
		my $file_name="$dir/$loc/COMMON_DATA/$name.dat";
		my $file= new IO::File "< $file_name" || die("@! can't open $file_name");
		$/=undef;
		my $data=<$file>;
		$/="\n";
		$file->close;
		my $VAR1;
		eval($data);

		while(my ($k,$v)=each(%{$VAR1})){
			$genehas{$k}=$v;
			$k=~s/[a-z]$//;
			$genehash{$k}=$v;
		}
	}
	return \%genehash;
}
