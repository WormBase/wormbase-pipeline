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
    my @xrefs = @{ $gene->get_all_xrefs};

    print "Gene : \"$geneName\"\nPublic_name \"$publicName\"\nSpecies \"$species\"\n";
    map {print "Other_name \"$_\" Inferred_automatically \"uniprot_xrefs_from_ensembl\"\n"} grep {$_->dbname eq 'Uniprot_gn_gene_name'} @xrefs;
    print "Database EnsEMBL ENSEMBL_geneID \"$geneName\"\n";
    dbxrefS('UniProt','UniProt_AC',grep {$_->dbname =~ /^Uniprot\//} @xrefs);
    dbxrefS('UniProt','UniProtACC',grep {$_->dbname =~ /^Uniprot\//} @xrefs);
    dbxrefS('ZFIN','acc',grep {$_->dbname =~ /^ZFIN/} @xrefs);
    dbxrefS('SGD','acc',grep {$_->dbname =~ /^SGD/} @xrefs);
    dbxrefS('MGI','acc',grep {$_->dbname =~ /^MGI/} @xrefs);
    dbxrefS('FLYBASE','FlyBase_gn',grep {$_->dbname eq 'FlybaseCGID_gene'} @xrefs);
    dbxrefS('FLYBASE','FlyBase_ID',grep {$_->dbname eq 'flybase_gene_id'} @xrefs);
    dbxrefS('HGNC','symbol',grep{$_->dbname =~ /^HGNC/} @xrefs);
    dbxrefP('OMIM','gene',grep{$_->dbname eq 'MIM_GENE'} @xrefs);
    dbxrefP('OMIM','disease',grep{$_->dbname eq 'MIM_DISEASE'} @xrefs);
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

