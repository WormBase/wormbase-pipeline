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

my ($comparadb,$dbhost,$dbuser,$dbport);

GetOptions(
	   'database=s' => \$comparadb,
	   'dbhost=s'   => \$dbhost,
	   'dbuser=s'   => \$dbuser,
	   'dbport=s'   => \$dbport,
) || die("cant parse the command line parameter\n");

$comparadb ||= 'worm_compara';
$dbhost    ||= $ENV{'WORM_DBHOST'};
$dbuser    ||= 'wormro';
$dbport    ||= $ENV{'WORM_DBPORT'};

# let's be fair, that can be probably build from the WormBase.pm
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
    860376 => 'Caenorhabditis angaria',
    6334   => 'Trichinella spiralis',
    886184 => 'Caenorhabditis tropicalis',
    6253   => 'Ascaris suum',
    37862  => 'Heterorhabditis bacteriophora',
    6326   => 'Bursephelenchus xylophilus',
    497829 => 'Caenorhabditis sp.5',
    34506  => 'Strongyloides ratti',
    7209   => 'Loa loa',
    6233   => 'Panagrellus redivivus',
    6287   => 'Dirofilaria immitis',
    6282   => 'Onchocerca volvulus',
    51031  => 'Necator americanus',
);

my %coreSpecies = (
   6239   => 1,
   6238   => 1,
   31234  => 1,
   135651 => 1,
   281687 => 1,
   54126  => 1,
   6279   => 1,
   6282   => 1,
);

my %cds2wbgene=%{&get_commondata('cds2wbgene_id')};

my $compara_db = new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -host   => $dbhost,
    -user   => $dbuser,
    -port   => $dbport,
    -dbname => $comparadb
) or die(@!);

my $member_adaptor = $compara_db->get_GeneMemberAdaptor();
my $homology_adaptor = $compara_db->get_HomologyAdaptor();
my @members = @{$member_adaptor->fetch_all()};

while( my $member = shift @members){
    
    next unless $coreSpecies{$member->taxon_id};

    my @homologies = @{$homology_adaptor->fetch_all_by_Member( $member)};

    my (%t3,%t2);

    foreach my $homology ( @homologies) {
        
        next if $homology->description eq 'between_species_paralog';
        foreach my $ma ( @{ $homology->get_all_Members } ) { # was $homology->get_all_Member_Attribute but this is not in the API any more
            
            foreach my $pepm ( @{ $ma->gene_member()->get_all_SeqMembers() } ) { # was $me->get_all_peptide_Members()

                if ($coreSpecies{$pepm->taxon_id}){ 
                    $t2{$pepm->stable_id} = [$pepm->taxon_id,$homology->description]
                } else {
		    $t3{$pepm->stable_id} = [$pepm->taxon_id,$homology->description,$pepm->sequence]
	        }
            }
        }
    }

    my $gid=$cds2wbgene{$member->stable_id}?$cds2wbgene{$member->stable_id}:$member->stable_id;

    next unless (scalar(keys %t2) + scalar(keys %t3) > 1); # the self-hit should always exist

    print "Gene : \"$gid\"\n";
    
    while (my ($k,$v)=each(%t2)){
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

    while (my ($k,$v)=each(%t3)){
            print "Ortholog_other $k From_analysis WormBase-Compara\n";
    }
    print "\n";
    
    while (my ($k,$v)=each(%t3)){
	    print "Protein : $k\nSpecies \"$species{$$v[0]}\"\nPeptide \"$k\"\n\n";
            print "Peptide : \"$k\"\n$$v[2]\n\n"
    }
}   

#################
# also adds the sequence name of the parent gene
sub get_commondata {
    my ($name)=@_;
    my %genehash;
    my @locations=qw(autoace remanei briggsae pristionchus japonica brenneri brugia ovolvulus);
    my $dir='/nfs/panda/ensemblgenomes/wormbase/BUILD/';
      
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
