#!/usr/bin/env perl
#===============================================================================
#   updates ortholog connections from the main Compara database
#
#       AUTHOR:   (), <>
#      VERSION:  1.0
#      CREATED:  03/08/06 13:26:19 BST
#     REVISION:  ---
#===============================================================================

use lib $ENV{CVS_DIR};

use strict;
use Getopt::Long;
use Storable;
use Wormbase;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

my ($wormbase, $debug, $test, $verbose, $store, $species, $use_current_db, $compara, %cds2wbgene, %current_genes);

GetOptions (
  "debug=s"      => \$debug,
  "test"         => \$test,
  "verbose"      => \$verbose,
  "store:s"      => \$store,
  "species:s"    => \$species,
  "usecurrentdb" => \$use_current_db,
  "compara=s"    => \$compara,
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
$compara = "multi" if not defined $compara;

Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ensembldb.ensembl.org', 
                                              -user => 'anonymous',
                                              -port => 5306,
                                              -verbose => 0,
    );

my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species_name,'core','Gene');
my $member_adaptor   = Bio::EnsEMBL::Registry->get_adaptor($compara,'compara','GeneMember');
my $homology_adaptor = Bio::EnsEMBL::Registry->get_adaptor($compara,'compara','Homology');

die "Could not fetch compara homology adaptor('$compara')" if not defined $homology_adaptor;
die "Could not fetch compara member adaptor('$compara')" if not defined $member_adaptor;


if ($use_current_db) {
  %cds2wbgene = %{&get_commondata( "/nfs/DATABASES/current_db/COMMON_DATA/cds2wbgene_id.dat" )};
} else {
  %cds2wbgene = %{&get_commondata( $wormbase->common_data . "/cds2wbgene_id.dat" )};
}
map { $current_genes{$_} = 1 } values %cds2wbgene;

my $genes = $gene_adaptor->fetch_all_by_biotype('protein_coding');
foreach my $other_species ('Homo sapiens', 
                           'Mus musculus', 
                           'Drosophila melanogaster', 
                           'Saccharomyces cerevisiae') {
  my %homols_for_species;

  foreach my $gene (@$genes) {
    my $wb_gene_id = $gene->stable_id;

    next if not exists $current_genes{$wb_gene_id};
    
    my $member = $member_adaptor->fetch_by_stable_id( $gene->stable_id());
    next unless $member;
    my $homologies = $homology_adaptor->fetch_all_by_Member_paired_species($member, $other_species, ['ENSEMBL_ORTHOLOGUES']);
    
    my %homols_for_gene;
    
    foreach my $homology ( @{$homologies} ) {
      foreach my $me (grep { $_->stable_id ne $wb_gene_id } @{$homology->get_all_GeneMembers}) {
        $homols_for_gene{$me->stable_id} = 1;
        $homols_for_species{$me->stable_id} = 1;
      }
    }
    
    next unless %homols_for_gene;
        
    print "\nGene : \"$wb_gene_id\"\n";
    while (my ($k,$v)=each(%homols_for_gene)){
      printf "Ortholog \"%s\" \"%s\" From_analysis EnsEMBL-Compara\n", $k, $other_species;
    }
  }

  while (my ($k,$v)=each(%homols_for_species)){
    print "\nGene : \"$k\"\n";
    print "Species \"$other_species\"\n";
    print "DB_info Database EnsEMBL ENSEMBL_geneID $k\n";
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
