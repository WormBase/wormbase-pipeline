#!/usr/bin/perl -w
#===============================================================================
#   updates ortholog connections from the main Compara database
#
#       AUTHOR:   (), <>
#      VERSION:  1.0
#      CREATED:  03/08/06 13:26:19 BST
#     REVISION:  ---
#===============================================================================

use lib '/software/worm/ensembl/bioperl-live';
use lib $ENV{CVS_DIR};

use strict;
use IO::File;
use Getopt::Long;
use Storable;
use Wormbase;



#use lib '/software/worm/ensembl-64/ensembl/modules';
#use lib '/software/worm/ensembl-64/ensembl-compara/modules';

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;


my %speciesOfInterest = (
  'saccharomyces_cerevisiae'  => 'yeast',
  'drosophila_melanogaster' => 'fly',
  'homo_sapiens'            => 'human',
  'mus_musculus'            => 'mouse',
);


my ($wormbase, $debug, $test, $verbose, $store, $species, $use_current_db, %cds2wbgene, %current_genes);

GetOptions (
  "debug=s"    => \$debug,
  "test"       => \$test,
  "verbose"    => \$verbose,
  "store:s"    => \$store,
  "species:s"  => \$species,
  "usecurrentdb" => \$use_current_db,
            );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species
                             );
}

my $species_name = $wormbase->full_name;

Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ens-livemirror', 
                                              -user => 'wormro',
                                              -port => 3306,
                                              -verbose => 0);


my $slice_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species_name,'core','Slice');
my $member_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','Member');
my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor('Multi','compara','Homology');;


if ($use_current_db) {
  %cds2wbgene = %{&get_commondata( "/nfs/DATABASES/current_db/COMMON_DATA/cds2wbgene_id.dat" )};
} else {
  %cds2wbgene = %{&get_commondata( $wormbase->common_data . "/cds2wbgene_id.dat" )};
}
map { $current_genes{$_} = 1 } values %cds2wbgene;


my @slices = @{$slice_adaptor->fetch_all('toplevel')};


foreach my $slice(@slices){
  my $genes = $slice->get_all_Genes;
  while (my $gene = shift @{$genes}){
    my ($db_entry) = @{$gene->get_all_DBLinks('wormbase_gene')};
    my $wb_gene_id = $db_entry->primary_id;

    next if not exists $current_genes{$wb_gene_id};
    
    my $member = $member_adaptor->fetch_by_source_stable_id( 'ENSEMBLGENE',$gene->stable_id());
    next unless $member;
    my $homologies = $homology_adaptor->fetch_all_by_Member_method_link_type( $member, 'ENSEMBL_ORTHOLOGUES' );
    
    my %homol_ids;
    my %omims;
    my %omim_genes;
    
    foreach my $homology ( @{$homologies} ) {
      foreach my $ma ( @{ $homology->get_all_Member_Attribute } ) {
        my ( $me, $at ) = @{$ma};
        if (exists $speciesOfInterest{$me->genome_db->name}){
          my $pepm = $me->get_canonical_Member(); 
          $homol_ids{ $pepm->stable_id } = $pepm; #$homology->description,$homology->subtype];
        }
      }
    }
    
    next unless %homol_ids;
        
    print "Gene : \"$wb_gene_id\"\n";
    while (my ($k,$v)=each(%homol_ids)){
      print "Ortholog_other ENSEMBL:$k  From_analysis EnsEMBL-Compara\n";
    }
    print "\n";

    while (my ($k,$v)=each(%homol_ids)){
      # put the OMIM orthologs into another hash based on the MIM_GENE xref
      print "Protein : ENSEMBL:$k\n";
      print "Species \"${\$v->taxon->name}\"\n";
      print "DB_info Database EnsEMBL ENSEMBL_proteinID $k\n";
      print "\n";
    }

  }	
}


sub get_commondata {
  my ($fname)=@_;
  
  open(my $fh, $fname) or die "could not open $fname for reading\n";
  $/=undef;
  my $data=<$fh>;
  $/="\n";
  close($fh);

  my $VAR1;
  my %genehash;
  eval($data);
  
  while(my ($k,$v)=each(%{$VAR1})){
    $k=~s/[a-z]$//;
    $genehash{$k}=$v
  }
  return \%genehash;
}
