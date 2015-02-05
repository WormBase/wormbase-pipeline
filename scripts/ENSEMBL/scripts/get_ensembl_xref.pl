#!/usr/bin/perl -w
# creates gene stubs from the core database on mysql-ps-prod-1
# currently: fish, fly, mouse, human, yeast but you can pass others as commandline

use lib $ENV{CVS_DIR};

use IO::File;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'mysql-ps-prod-1.ebi.ac.uk', 
                                              -user => 'ensro',
                                              -port => 4450,
                                              -verbose => 0);

my @names = ('danio_rerio','drosophila_melanogaster','mus_musculus','homo_sapiens','saccharomyces_cerevisiae');

GetOptions(
     'productionNames=s' => \@names,
)||die(@!);


map {dump_ace_by_production_name($_)} @names;

sub dump_ace_by_production_name{
  my ($prodName)=@_;
  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor($prodName,'core','Gene');
  my $meta_adaptor = Bio::EnsEMBL::Registry->get_adaptor($prodName,'core','MetaContainer');

  my $species = (($meta_adaptor->single_value_by_key('species.scientific_name'))=~/(\w+\s+\w+)/g)[0];

  foreach my $gene (@{$gene_adaptor->fetch_all_by_biotype('protein_coding')}) {

    my $publicName = $gene->display_xref ? $gene->display_xref->display_id : $gene->stable_id;
    my $geneName = $gene->stable_id;

    print "Gene : \"$geneName\"\nPublic_name \"$publicName\"\nSpecies \"$species\"\n";
    print "Database EnsEMBL ENSEMBL_geneID \"$geneName\"\n";
    dbxrefS('UniProt','UniProt_AC',@{$gene->get_all_xrefs('Uniprot/%')});
    dbxrefS('UniProt','UniProtACC',@{$gene->get_all_xrefs('Uniprot/%')});
    dbxrefS('ZFIN','acc',@{ $gene->get_all_xrefs('ZFIN%')});
    dbxrefS('SGD','acc',@{ $gene->get_all_xrefs('SGD%')});
    dbxrefS('MGI','acc',@{ $gene->get_all_xrefs('MGI%')});
    dbxrefS('FLYBASE','FlyBase_gn',@{ $gene->get_all_xrefs('FlybaseCGID_gene')});
    dbxrefS('FLYBASE','FlyBase_ID',@{ $gene->get_all_xrefs('flybase_gene_id')});
    dbxrefS('HGNC','symbol',@{ $gene->get_all_xrefs('HGNC%')});
    dbxrefP('OMIM','gene',@{ $gene->get_all_xrefs('MIM_GENE')});
    dbxrefP('OMIM','disease',@{ $gene->get_all_xrefs('MIM_DISEASE')});
    print "\n";
  }	
}

# display_id xref printing
sub dbxrefS{
  my ($db,$ac,@xrefs)=@_;
  my %seen; 
  map {  my $p = $_->info_type eq 'PROJECTION'?'projected_':'';
         printf "Database $db $p$ac \"%s\"\n",$_->display_id unless $seen{$_->display_id} 
   } @xrefs;
}

# display_id xref printing
sub dbxrefP{
  my ($db,$ac,@xrefs)=@_;
  my %seen; 
  map {  my $p = $_->info_type eq 'PROJECTION'?'projected_':'';
         printf "Database $db $p$ac \"%s\"\n",$_->primary_id unless $seen{$_->primary_id} 
  } @xrefs;
}

