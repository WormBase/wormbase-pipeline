#!/usr/bin/env perl
# create a MISC_DATA acefile

use IO::File;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

my %cds2wbgeneid = %{slurp('~mh6/ebi_home/project/DATABASES/WS238/COMMON_DATA/cds2wbgene_id.dat')};

my %worms = (
     'Caenorhabditis elegans' => 6239,
     'Caenorhabditis briggsae'=> 6238,
     'Caenorhabditis remanei' => 31234,
     'Caenorhabditis brenneri'=> 135651,
     'Caenorhabditis japonica'=> 281687,
     'Pristionchus pacificus' => 5126,
);

my $comparadb= new Bio::EnsEMBL::Compara::DBSQL::DBAdaptor(
    -host     => 'mysql-treefam-prod',
    -port     => 4401,
    -user     => 'treefam_ro',
    -dbname   => 'treefam_production_9_69',
    -pass     => 'treefam_ro',
);

my $memberAdaptor=$comparadb->get_adaptor('Member');
my $homologyAdaptor=$comparadb->get_adaptor('Homology');
my %taxID2tierII = reverse %worms;

foreach my $taxId(values %worms){
     my $members = $memberAdaptor->fetch_all_by_source_taxon('ENSEMBLGENE',$taxId);
     foreach my $member (@$members){
        my $homologies = $homologyAdaptor->fetch_all_by_Member_method_link_type( $member, 'ENSEMBL_ORTHOLOGUES' );
        my $paralogies = $homologyAdaptor->fetch_all_by_Member_method_link_type( $member, 'ENSEMBL_PARALOGUES' );

        my %homol_ids;

        my %paralogs;
        my %orthologs;

        print "processing: ",$member->stable_id , "\n" if $ENV{DEBUG};
        foreach my $homology ( @{$homologies} , @$paralogies ) {
          foreach my $ma ( @{ $homology->get_all_Member_Attribute } ) {
            my ( $m, $at ) = @{$ma};
            my $me = $m->gene_member(); 
            print "-",$me->stable_id,"\n" if $ENV{DEBUG};
            if ($taxID2tierII{$me->genome_db->taxon_id}){
               my $realID = $cds2wbgeneid{$me->stable_id};
               next if $realID eq $cds2wbgeneid{$member->stable_id};
               next unless $realID;
               if ($me->genome_db->taxon_id == $taxId){
                 $paralogs{$realID} = $me->taxon->name;
               }else{
                 $orthologs{$realID} = $me->taxon->name;
               }
            }
          }
        }

        next unless (%paralogs || %orthologs);
        next unless $cds2wbgeneid{$member->stable_id};

        print "Gene : ${\$cds2wbgeneid{$member->stable_id}}\n";
        while (my ($k,$v)=each(%orthologs)){
           print "Ortholog $k \"$v\" From_analysis TreeFam\n";
        }
        while (my ($k,$v)=each(%paralogs)){
           print "Paralog $k \"$v\" From_analysis TreeFam\n";
        }

        print "\n";
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
 
   #cds2sequenceNameConversion
   while(my($k,$v)=each %$VAR1){
     $k=~s/[a-z]+$//;
     $$VAR1{$k}||=$v;
   }   
  
   return $VAR1;
}
