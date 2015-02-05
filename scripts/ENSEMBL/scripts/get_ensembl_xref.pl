#!/usr/bin/perl -w
# creates gene stubs from the core database on mysql-ps-prod-1
# currently: fish, fly, mouse, human, yeast but you can pass others as commandline

use lib $ENV{CVS_DIR};
use 5.010;
use IO::File;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

Bio::EnsEMBL::Registry->load_registry_from_db(-host => 'mysql-ps-prod-1.ebi.ac.uk', 
                                              -user => 'ensro',
                                              -port => 4450,
                                              -verbose => 0);

my @names;
GetOptions(
     'productionNames=s' => \@names,
)||die(@!);
@names = ('danio_rerio','drosophila_melanogaster','mus_musculus','homo_sapiens','saccharomyces_cerevisiae') unless $names[0];


map {dump_ace_by_production_name($_)} @names;

sub dump_ace_by_production_name{
  my ($prodName)=@_;
  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor($prodName,'core','Gene');
  my $meta_adaptor = Bio::EnsEMBL::Registry->get_adaptor($prodName,'core','MetaContainer');

  my $species = (($meta_adaptor->single_value_by_key('species.scientific_name'))=~/(\w+\s+\w+)/g)[0];
  
  my @genes=  @{$gene_adaptor->fetch_all_by_biotype('protein_coding')};
  while( my $gene = shift @genes) {

    my $publicName = $gene->display_xref ? $gene->display_xref->display_id : $gene->stable_id;
    my $geneName = $gene->stable_id;

    print "Gene : \"$geneName\"\nPublic_name \"$publicName\"\nSpecies \"$species\"\n";
    print "Database EnsEMBL ENSEMBL_geneID \"$geneName\"\n";

    my @xrefs =  (@{ $gene->get_all_xrefs},$gene->display_xref); # because some of the displayxrefs are not in the objext xref table :-(
    my %seen; # like wsb2
    while (my $x = shift @xrefs){
      next if $seen{$x->dbID};
      $seen{$x->dbID}=1;
      given($x->dbname){
         when(/^Uniprot\//){dbxrefS('UniProt','UniProt_AC',$x)}
         when('ZFIN_ID'){dbxrefP('ZFIN','primary_acc',$x)}
         when('SGD_GENE'){dbxrefS('SGD','acc',$x)}
         when('MGI'){dbxrefS('MGI','acc',$x)}
         when('FlybaseCGID_gene'){dbxrefS('FLYBASE','FlyBase_gn',$x)}
         when('flybase_gene_id'){dbxrefS('FLYBASE','FlyBase_ID',$x)}
         when('HGNC'){dbxrefS('HGNC','symbol',$x)}
         when('MIM_GENE'){dbxrefP('OMIM','gene',$x)}
         when('MIM_DISEASE'){dbxrefP('OMIM','disease',$x)}
         when('Uniprot_gn'){printf "Other_name \"%s\" Inferred_automatically \"uniprot_xrefs_from_ensembl\"\n",$x->display_id}
         default{}
      }
    }
    print "\n";
  }	
}

#  display_id xref printing
sub dbxrefS{
  my ($db,$ac,@xrefs)=@_;
  map {  my $p = $_->info_type eq 'PROJECTION'?'projected_':'';
         printf "Database $db $p$ac \"%s\"\n",$_->display_id 
   } @xrefs;
}

# primary_id xref printing
sub dbxrefP{
  my ($db,$ac,@xrefs)=@_;
  map {  my $p = $_->info_type eq 'PROJECTION'?'projected_':'';
         printf "Database $db $p$ac \"%s\"\n",$_->primary_id
  } @xrefs;
}

