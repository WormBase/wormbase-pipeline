#!/software/worm/perl_510/bin/perl
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

my %species = ( 
    6239   => 'Caenorhabditis elegans',
    6238   => 'Caenorhabditis briggsae', 
    31234  => 'Caenorhabditis remanei',
    135651 => 'Caenorhabditis brenneri',
    281687 => 'Caenorhabditis japonica',
    6279   => 'Brugia malayi',
    54126  => 'Pristionchus pacificus',
    6305   => 'Meloidogyne hapla',
    6289   => 'Haemonchus contortus',
    96668  => 'Caenorhabditis angaria',
);

my %cds2wbgene=%{&get_commondata('cds2wbgene_id')};
my %cds2swiss=%{&get_cds2swiss('brugia')};

$comparadb||='worm_compara';

my $compara_db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -host   => 'farmdb1',        # change that
    -user   => 'wormro',       # and that
    -dbname => $comparadb
);


my $member_adaptor = $compara_db->get_MemberAdaptor();
my $homology_adaptor = $compara_db->get_HomologyAdaptor();
my @members = @{$member_adaptor->fetch_all()};

while( my $member = shift @members){
    
    next if $member->taxon_id == 6279; # skip brugia, as it does not have ?Gene objects
    next if $member->taxon_id == 6305; # skip M.hapla, as it does not have ?Gene objects
    next if $member->taxon_id == 6289; # skip H.contortus, as it does not have ?Gene objects
    next if $member->taxon_id == 96668;# skip C.angaria, because it does not have ?Gene objects

    my @homologies = @{$homology_adaptor->fetch_all_by_Member( $member)};

    # needs some better way
    my %brugia;
    my %hapla;
    my %hcont;
    my %cang;
    my %all;

    foreach my $homology ( @homologies) {
        
        next if $homology->description eq 'between_species_paralog';
        foreach my $ma ( @{ $homology->get_all_Member_Attribute } ) {
            
            my ( $me, $at ) = @{$ma};
            foreach my $pepm ( @{ $me->get_all_peptide_Members() } ) {
                
                if ($pepm->taxon_id==6279){
                    $brugia{ $pepm->stable_id } = [$pepm->taxon_id,$homology->description] 
                }
		elsif ($pepm->taxon_id==6305){
                    $hapla{$pepm->stable_id} = [$pepm->taxon_id,$homology->description] 
                }
                elsif ($pepm->taxon_id==6289){
		    $hcont{$pepm->stable_id} = [$pepm->taxon_id,$homology->description]
		}
	    	elsif ($pepm->taxon_id==96668){
		    $cang{$pepm->stable_id} = [$pepm->taxon_id,$homology->description]
		}
                else {
                    $all{ $pepm->stable_id } = [$pepm->taxon_id,$homology->description]
                }
            }

        }
    }

    my $gid=$cds2wbgene{$member->stable_id}?$cds2wbgene{$member->stable_id}:$member->stable_id;

    next unless (scalar keys %all > 1);

    print "Gene : \"$gid\"\n";
    
    while (my ($k,$v)=each(%all)){
            my $sid=$cds2wbgene{$k}?$cds2wbgene{$k}:$k;
            next if $gid eq $sid;
            
            if ($$v[0] != $member->taxon_id){
                print "Ortholog $sid \"$species{$$v[0]}\" From_analysis WormBase-Compara";
            }
            else{
                print "Paralog $sid \"$species{$$v[0]}\" From_analysis WormBase-Compara";
            }
            print " // Type $$v[1]" if $ENV{debug};
            print "\n";
    }

    my %t3 = (%brugia,%hapla,%hcont,%cang);
    while (my ($k,$v)=each(%t3)){ # brugia exception
            my $bid=$cds2swiss{$k}?$cds2swiss{$k}:$k;
            print "Ortholog_other $bid From_analysis WormBase-Compara\n";
    }
    print "\n";
    

    while (my ($k,$v)=each(%t3)){
	    print "Protein : $k\nSpecies \"$species{$$v[0]}\"\n\n";
    }

}   

#################


# needs a merging step
# /nfs/wormpub/BUILD/autoace/COMMON_DATA/
# /nfs/wormpub/BUILD/remanei/COMMON_DATA/
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
