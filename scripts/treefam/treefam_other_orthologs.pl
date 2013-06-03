#!/usr/bin/env perl
# create a MISC_DATA acefile

use IO::File;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

my %cds2wbgeneid = %{slurp('~mh6/project/DATABASES/WS238/COMMON_DATA/cds2wbgene_id.dat')};

my %worms = (
     'Caenorhabditis elegans' => 6239,
     'Caenorhabditis briggsae'=> 6238,
     'Caenorhabditis remanei' => 31234,
     'Caenorhabditis brenneri'=> 135651,
     'Caenorhabditis japonica'=> 281687,
     'Pristionchus pacificus' => 5126,
);

my %interesting_species = (
     7227  => 'Drosophila melanogaster',
     9606  => 'Homo sapiens',
     10090 => 'Mus musculus',
     7955  => 'Danio rerio',
);

my $comparadb= new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -host     => 'mysql-treefam-prod',
    -port     => 4401,
    -user     => 'treefam_ro',
    -dbname   => 'treefam_production_9_69',
    -pass     => 'treefam',
);

my $memberAdaptor=$comparadb->get_adaptor('Member');
my $homologyAdaptor=$comparadb->get_adaptor('Homology');

foreach my $taxId(values %worms){
     my $members = $memberAdaptor->fetch_all_by_source_taxon('ENSEMBLGENE',$taxId);
     foreach my $member (@$members){
        my $homologies = $homologyAdaptor->fetch_all_by_Member_method_link_type( $member, 'ENSEMBL_ORTHOLOGUES' );
        my %homol_ids;
        my %paralogs;
        my %orthologs;
        print "processing: ",$member->stable_id , "\n" if $ENV{DEBUG};
        foreach my $homology ( @{$homologies} ) {
          foreach my $ma ( @{ $homology->get_all_Member_Attribute } ) {
            my ( $me, $at ) = @{$ma};
            print "-",$me->stable_id,"\n" if $ENV{DEBUG};
            if ($interesting_species{$me->genome_db->taxon_id}){
               $homol_ids{ $me->stable_id } = $me->taxon->name;
            }
          }
        }
        next unless %homol_ids;
        next unless $cds2wbgeneid{$member->stable_id};
        print "Gene : ${\$cds2wbgeneid{$member->stable_id}}\n";
        while (my ($k,$v)=each(%homol_ids)){
           print "Ortholog_other ENSEMBL:$k  From_analysis TreeFam\n";
        }
        print "\n";

        while (my ($k,$v)=each(%homol_ids)){
           print "Protein : ENSEMBL:$k\n";
           print "Species \"$v\"\n";
           print "DB_info Database EnsEMBL ENSEMBL_proteinID $k\n";
           print "\n";
        }
     }
}

sub slurp{
   my ($file) = @_;
   undef $/;
   my $VAR1;
   my $inf = IO::File->new(glob($file),"r") or die @!;
   my $txt = <$inf>;
   $inf->close;
   eval $txt;
   $/="\n";
   return $VAR1;
}
