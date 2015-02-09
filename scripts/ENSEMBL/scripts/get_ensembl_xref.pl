#!/usr/bin/perl -w
# creates gene stubs from the core database on mysql-ps-prod-1
# currently: fish, fly, mouse, human, yeast. but you can pass others as commandline

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

  my $species = (($meta_adaptor->single_value_by_key('species.scientific_name'))=~/(\w+\s+\w+)/g)[0]; # yeast special case
  
  my @genes=  @{$gene_adaptor->fetch_all_by_biotype('protein_coding')};
  while( my $gene = shift @genes) {

    my $publicName = $gene->display_xref ? $gene->display_xref->display_id : $gene->stable_id;
    $publicName =~s/;//;
    my $geneName = $gene->stable_id;

    print "Gene : \"$geneName\"\nPublic_name \"$publicName\"\n";
  
    if ($gene->display_xref && print_evidence($gene->display_xref)){print "CGC_name \"$publicName\" ",print_evidence($gene->display_xref),"\n"}

    print "Species \"$species\"\nDatabase EnsEMBL ENSEMBL_geneID \"$geneName\"\n";

    my @xrefs =  (@{ $gene->get_all_xrefs},$gene->display_xref); # because some of the displayxrefs are not in the objext xref table :-(
    my %seen; # like wsb2
    while (my $xref = shift @xrefs){
      next if $seen{$xref->dbID};
      my $dbline = &print_xref($xref,$publicName);
      print "$dbline\n" if $dbline;
      $seen{$xref->dbID}=1;
    }
    print "\n";
  }	
}
sub print_evidence{
      my($x)=@_;
      if($x->dbname =~/^Uniprot\//){return 'Accession_evidence "UniProt" "'.$x->primary_id.'"'}
      elsif($x->dbname eq 'ZFIN_ID'){return 'Accession_evidence "ZFIN" "'.$x->primary_id.'"'}
      elsif($x->dbname eq 'SGD_GENE'){return 'Accession_evidence "SGD" "'.$x->display_id.'"'}
      elsif($x->dbname eq 'MGI'){return'Accession_evidence "MGI" "'.$x->display_id.'"'}
      elsif($x->dbname eq 'FlyBaseCGID_gene'){return 'Accession_evidence "FLYBASE" "'.$x->display_id.'"'}
      elsif($x->dbname eq 'flybase_gene_id'){return 'Accession_evidence "FLYBASE" "'.$x->display_id.'"'}
      elsif($x->dbname eq 'HGNC'){return 'Accession_evidence "HGNC" "'.$x->display_id.'"'}
      elsif($x->dbname eq 'MIM_GENE'){return 'Accession_evidence "OMIM" "'.$x->primary_id.'"'}
      elsif($x->dbname eq 'MIM_MORBID'){return 'Accession_evidence "OMIM" "'.$x->primary_id.'"'}
      else{
          my $tmp=$x->dbname;
          $tmp=~s/;//;
          return 'Inferred_automatically "'.$x->dbname.'"'
      }
}


sub print_xref{
      my($x,$p)=@_;
      my $xs;
      if($x->dbname =~/^Uniprot\//){$xs = dbxrefP('UniProt','UniProt_AC',$x)}
      elsif($x->dbname eq 'ZFIN_ID'){$xs = dbxrefP('ZFIN','primary_acc',$x)}
      elsif($x->dbname eq 'SGD_GENE'){$xs = dbxrefS('SGD','SGD_acc',$x)}
      elsif($x->dbname eq 'MGI'){$xs = dbxrefS('MGI','MGI_acc',$x)         }
      elsif($x->dbname eq 'FlyBaseCGID_gene'){$xs = dbxrefS('FLYBASE','FlyBase_gn',$x)}
      elsif($x->dbname eq 'flybase_gene_id'){$xs = dbxrefS('FLYBASE','FlyBase_ID',$x)}
      elsif($x->dbname eq 'HGNC'){$xs = dbxrefS('HGNC','symbol',$x)}
      elsif($x->dbname eq 'MIM_GENE'){$xs = dbxrefP('OMIM','gene',$x)}
      elsif($x->dbname eq 'MIM_MORBID'){$xs = dbxrefP('OMIM','disease',$x)}
      elsif($x->dbname eq 'Uniprot_gn'){
           my $oname = $x->display_id;
           $oname=~s/;//; # strip semicolons from the IDs
           if($p && lc($p) ne lc($oname)){
            $xs = sprintf "Other_name \"$oname\" Accession_evidence \"UniProt\" \"$oname\"" if $p;
           }
      }
      return $xs;
}

#  display_id xref printing
sub dbxrefS{
  my ($db,$ac,$xref)=@_;
  my $p = $xref->info_type eq 'PROJECTION'?'projected_':'';
  return sprintf "Database \"$db\" \"$p$ac\" \"%s\"",$xref->display_id 
}

# primary_id xref printing
sub dbxrefP{
  my ($db,$ac,$xref)=@_;
  my $p = $xref->info_type eq 'PROJECTION'?'projected_':'';
  return sprintf "Database \"$db\" \"$p$ac\" \"%s\"",$xref->primary_id;
}

