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

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'ens-livemirror', 
                                              -user => 'wormro',
                                              -port => 3306,
                                              -verbose => 0);

my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor('homo_sapiens','core','Gene');

foreach my $gene (@{$gene_adaptor->fetch_all_by_biotype('protein_coding')}) {
  my @omim_disease = @{$gene->get_all_DBLinks('MIM_MORBID')};
  my @omim_gene = @{$gene->get_all_DBLinks('MIM_GENE')};

  if (@omim_gene or @omim_disease) {
    my $trans = $gene->canonical_transcript;
    my $trl = $trans->translation; 
        
    printf "Protein : \"ENSEMBL:%s\"\n", $trl->stable_id;
    print "Species \"Homo sapiens\"\n";
    printf "DB_info Database EnsEMBL ENSEMBL_proteinID %s\n", $trl->stable_id;
    printf "DB_info Database EnsEMBL ENSEMBL_transcriptID %s\n", $trans->stable_id;
    printf "DB_info Database EnsEMBL ENSEMBL_geneID %s\n", $gene->stable_id;
    foreach my $og (@omim_gene) {
      printf "DB_info Database OMIM gene %s\n", $og->primary_id;
    }
    foreach my $od (@omim_disease) {
      printf "DB_info Database OMIM disease %s\n", $od->primary_id;
    }
    print "\n";
  }	
}

