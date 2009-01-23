#!/software/worm/perl_510/bin/perl
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
	135651 => 'Caenorhabditis brenneri',
        281687 => 'Caenorhabditis japonica',
	6279   => 'Brugia malayi',
	54126  => 'Pristionchus pacificus',
);

my $config = ( YAML::LoadFile("$FindBin::Bin/../etc/ensembl_lite.conf") )->{'elegans'};
my %cds2wbgene=%{&get_commondata('cds2wbgene_id')};
my %cds2swiss=%{&get_cds2swiss('brugia')};

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

		# needs some better way
		my %briggsae_ids;
		my %remanei_ids;
		my %brenneri_ids;
		my %brugia_ids;
		my %elegans_ids;
		my %pristionchus_ids;
		my %japonica_ids;


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
					elsif ($pepm->taxon_id==54126){
						$pristionchus_ids{$pepm->stable_id} = [$pepm->taxon_id,$homology->description,$homology->subtype]
					} 
					elsif ($pepm->taxon_id==281687){
						$japonica_ids{$pepm->stable_id} = [$pepm->taxon_id,$homology->description,$homology->subtype]
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

                # could probably merge all hashes into one big one
                my %all_ortho = (%briggsae_ids,%remanei_ids,%pristionchus_ids,%japonica_ids,%brenneri_ids);
		
		print "Gene : \"$gid\"\n";
		while (my ($k,$v)=each(%all_ortho)){
			        my $sid=$cds2wbgene{$k}?$cds2wbgene{$k}:$k;
				print "Ortholog $sid \"$species{$$v[0]}\" From_analysis WormBase-Compara\n";
		}

		while (my ($k,$v)=each(%brugia_ids)){                                               # brugia exception
				my $bid=$cds2swiss{$k}?$cds2swiss{$k}:$k;
				print "Ortholog_other $bid From_analysis WormBase-Compara\n";
		}

                undef %all_ortho;
		print "\n";

                
		# briggsae	

		&print_thing(\%briggsae_ids,\%pristionchus_ids,\%remanei_ids,\%japonica_ids,\%brenneri_ids,\%brugia_ids,$gid,\%species);

		# remanei part
	        &print_thing(\%remanei_ids,\%briggsae_ids,\%pristionchus_ids,\%japonica_ids,\%brenneri_ids,\%brugia_ids,$gid,\%species);

		# pristionchus part
		&print_thing(\%pristionchus_ids,\%remanei_ids,\%briggsae_ids,\%japonica_ids,\%brenneri_ids,\%brugia_ids,$gid,\%species);

		# japonica part
                &print_thing(\%japonica_ids,\%pristionchus_ids,\%remanei_ids,\%briggsae_ids,\%brenneri_ids,\%brugia_ids,$gid,\%species);

		# brenneri part
                &print_thing(\%brenneri_ids,\%japonica_ids,\%pristionchus_ids,\%remanei_ids,\%briggsae_ids,\%brugia_ids,$gid,\%species);
	}	
}

################################### REFACTOR ME #######################
sub print_thing {
	my @genes = @_;
	my $gid=$genes[6];
	my %species=%{$genes[7]};
	while (my ($k,$v)=each(%{$genes[0]})){
			   	my $sid=$cds2wbgene{$k}?$cds2wbgene{$k}:$k;
				print "Gene : \"$sid\"\n";
				print "Ortholog $gid \"${\$species{ $config->{taxon_id}}}\" From_analysis WormBase-Compara\n"; # elegans

				while (my ($r_k,$r_v)=each(%{$genes[1]})){                                                   
				        my $rid=$cds2wbgene{$r_k}?$cds2wbgene{$r_k}:$r_k;
					print "Ortholog $rid \"$species{$$r_v[0]}\" From_analysis WormBase-Compara\n"
				}
                                while (my ($r_k,$r_v)=each(%{$genes[2]})){                                                    
				        my $rid=$cds2wbgene{$r_k}?$cds2wbgene{$r_k}:$r_k;
					print "Ortholog $rid \"$species{$$r_v[0]}\" From_analysis WormBase-Compara\n"
				}
                                while (my ($r_k,$r_v)=each(%{$genes[3]})){                                               
				        my $rid=$cds2wbgene{$r_k}?$cds2wbgene{$r_k}:$r_k;
					print "Ortholog $rid \"$species{$$r_v[0]}\" From_analysis WormBase-Compara\n"
				}
				while (my ($r_k,$r_v)=each(%{$genes[4]})){                                          
				        my $rid=$cds2wbgene{$r_k}?$cds2wbgene{$r_k}:$r_k;
					print "Ortholog $rid \"$species{$$r_v[0]}\" From_analysis WormBase-Compara\n"
				}

				while (my ($r_k,$r_v)=each(%{$genes[5]})){                                               # brugia exception
				        my $rid=$cds2swiss{$r_k}?$cds2swiss{$r_k}:$r_k;
					print "Ortholog_other $rid From_analysis WormBase-Compara\n"
				}

				print "\n";
		}
}
#################


# needs a merging step
# /nfs/disk100/wormpub/BUILD/autoace/COMMON_DATA/
# /nfs/disk100/wormpub/BUILD/remanei/COMMON_DATA/
sub get_commondata {
	my ($name)=@_;
	my %genehash;
	my @locations=qw(autoace remanei briggsae pristionchus japonica brenneri);
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
			$genehash{$k}=$v;
			$k=~s/[a-z]$//;
			$genehash{$k}=$v;
		}
	}
	return \%genehash;
}

sub get_cds2swiss {
	my ($name)=@_;
	my %cds2swiss;
	my $dir = glob("~wormpub/DATABASES/$name/swiss2cds.dat");
	my $file = new IO::File "<$dir" ||die ("@! cannot open $dir");
	while (<$file>){
	    chomp;
            my ($swiss,$cds)=split;
	    $cds2swiss{"${cds}A"}=$swiss;
	}
	return \%cds2swiss;
}
