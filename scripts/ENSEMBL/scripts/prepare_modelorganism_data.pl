#!/usr/bin/env perl -w
# creates gene stubs from the core database on mysql-ps-prod-1
# currently: fish, fly, mouse, human, yeast. but you can pass others as commandline

use lib $ENV{CVS_DIR};

use strict;

use Wormbase;
use Getopt::Long;
use Bio::EnsEMBL::Registry;

my %SPECIES_INFO = (danio_rerio              => { prefix => 'Dre-' },
                    mus_musculus             => { prefix => 'Mmu-' },
                    drosophila_melanogaster  => { prefix => 'Dme-' },
                    saccharomyces_cerevisiae => { prefix => 'Sce-' },
                    homo_sapiens             => { prefix => 'Hsa-' });

my (
  $wb, $debug, $test, $store,
  @server_urls, @mod_species, $outfile,
  $compara_name, $ace_orthologs, $ace_genes, $out_dir,
    );

&GetOptions(
  'urls=s@'       => \@server_urls,
  'compara'       => \$compara_name,
  "debug=s"       => \$debug,
  "test"          => \$test,
  "store=s"       => \$store,
  "aceorthologs"  => \$ace_orthologs,
  "acegenes"      => \$ace_genes,
  "modspecies=s"  => \@mod_species,
  "outfile=s"     => \$outfile,
)||die(@!);

die "You must supply a list of server URLs containing the compara database and comparator species core databases\n" if not @server_urls;

$compara_name = 'wbparasite' if not defined $compara_name;
@mod_species = keys %SPECIES_INFO if not @mod_species;

if ( $store ) {
  $wb = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wb = Wormbase->new( -debug   => $debug,
                       -test    => $test,
      );
}

my %cgc_names = $wb->FetchData('cgc_name2gene');
my %wb_accessors = $wb->species_accessors;
$wb_accessors{$wb->species} = $wb;

foreach my $url (@server_urls) {
  Bio::EnsEMBL::Registry->load_registry_from_url($url);
}

my $genomedb_adap = Bio::EnsEMBL::Registry->get_adaptor($compara_name,'compara','GenomeDB');
my $member_adap = Bio::EnsEMBL::Registry->get_adaptor($compara_name,'compara','GeneMember');
my $mlss_adap =  Bio::EnsEMBL::Registry->get_adaptor($compara_name,'compara','MethodLinkSpeciesSet');
my $homology_adap =  Bio::EnsEMBL::Registry->get_adaptor($compara_name,'compara','Homology');

my $fh;
if ($outfile) {
  open($fh, ">$outfile") or die "Could not open $outfile for writing\n";
} else {
  $fh = \*STDOUT;
}

foreach my $genome (@mod_species) {
    
  &write_ace_genes($genome, $fh) if $ace_genes;
  foreach my $wb_species (sort keys %wb_accessors) {
    &write_ace_orthologs($genome, $wb_species, $fh) if $ace_orthologs;
  }
}



###################################################
sub write_ace_genes {
  my ($genome, $outfh) = @_;

  my $gene_adaptor = Bio::EnsEMBL::Registry->get_adaptor($genome,'core','Gene');
  my $meta_adaptor = Bio::EnsEMBL::Registry->get_adaptor($genome,'core','MetaContainer');

  #my $ens_version = $meta_adaptor->get_schema_version();
  my $species = $meta_adaptor->single_value_by_key('species.scientific_name');
  die "Could not fetch species name for $genome\n" if not $species;
  # hack for yeast (which has sub-species name in its scientific name)
  $species =~ s/(\S+)_\s+(\S+)\s+\S+/$1 $2/; 


  my $gdb = $genomedb_adap->fetch_by_name_assembly($genome);

  my @mem = @{$member_adap->fetch_all_by_GenomeDB($gdb)};

  foreach my $mem (@mem) {
    my $gene = $gene_adaptor->fetch_by_stable_id( $mem->stable_id );
    my $gene_name = $gene->stable_id;

    my $public_name = $gene_name;
    if (defined $gene->display_xref) {
      $public_name =  $gene->display_xref->display_id;
      $public_name =~s/;//;
      if (exists $cgc_names{$public_name}) {
        $public_name = $SPECIES_INFO{$genome}->{prefix} . $public_name;
      }
    }
    
    print $outfh "\nGene : \"$gene_name\"\nPublic_name \"$public_name\"\n";
    
    if ($gene->display_xref) {
      my $evidence = get_evidence($gene->display_xref);
      if ($evidence) {
        print $outfh "CGC_name \"$public_name\" $evidence\n";
      }
    }
    
    print $outfh "Species \"$species\"\nDatabase EnsEMBL ENSEMBL_geneID \"$gene_name\"\n";
    
    # because some of the displayxrefs are not in the objext xref table 
    my @xrefs =  (@{ $gene->get_all_xrefs},$gene->display_xref); 
    my %seen; # like wsb2
    while (my $xref = shift @xrefs){
      next if $seen{$xref->dbID};
      my $dbline = &get_xref($xref,$public_name);
      print $outfh "$dbline\n" if $dbline;
      $seen{$xref->dbID}=1;
    }
    print $outfh "\n";
  }
}


###################################################
sub write_ace_orthologs {
  my ($m_genome, $wb_genome, $outfh) = @_;

  my $meta_adaptor = Bio::EnsEMBL::Registry->get_adaptor($m_genome,'core','MetaContainer');

  my $m_gdb = $genomedb_adap->fetch_by_name_assembly($m_genome);
  die "Could not fetch genome db for $m_genome\n" if not $m_gdb;
  my $m_species = $meta_adaptor->single_value_by_key('species.scientific_name');
  die "Could not fetch species name for $m_genome\n" if not $m_species;
  # hack for yeast (which has sub-species name in its scientific name)
  $m_species =~ s/(\S+)_\s+(\S+)\s+\S+/$1 $2/; 

  # need to deal with WormBase and parasite having different production names for same genome. Grr...
  my $wb_species = $wb_accessors{$wb_genome}->full_name;
  my $prod_1 = lc($wb_species);  $prod_1 =~ s/\s+/_/; 
  my $prod_2 = join("_", $prod_1, lc($wb_accessors{$wb_genome}->ncbi_bioproject));
  
  my $wb_gdb;
  foreach my $prod_name ($prod_2, $prod_1) {
    eval {
      $wb_gdb = $genomedb_adap->fetch_by_registry_name( $prod_name );
    };
    last if not $@ and defined $wb_gdb;
  }
  if (not defined $wb_gdb) {
    die "Could not find genome_db as $prod_1 or $prod_2\n";
  }

  my $mlss = $mlss_adap->fetch_by_method_link_type_GenomeDBs('ENSEMBL_ORTHOLOGUES', [$m_gdb, $wb_gdb]);

  my %live_genes = $wb_accessors{$wb_genome}->FetchData('wbgene_id2cds');

  my %homols;

  my @homols = @{$homology_adap->fetch_all_by_MethodLinkSpeciesSet($mlss)};
  print STDERR "Fetched ", scalar(@homols), " orthologs between $m_genome and $wb_genome\n";
  my $i = 0;
  while (my $hom = shift @homols) {
    my ($m_gm, $wb_gm) =  map { $_->gene_member } @{$hom->get_all_Members};
    
    if ($m_gm->genome_db->dbID != $m_gdb->dbID) {
      ($m_gm, $wb_gm) = ($wb_gm, $m_gm);
    }
    
    next if not exists $live_genes{$wb_gm->stable_id};
    #printf STDERR " Registering homol...%d...\n", $i++;
    $homols{$m_gm->stable_id}->{$wb_species}->{$wb_gm->stable_id} = 1;
    $homols{$wb_gm->stable_id}->{$m_species}->{$m_gm->stable_id} = 1;
  }

  foreach my $g (keys %homols) {
    print $outfh "Gene : \"$g\"\n";
    foreach my $spe (keys %{$homols{$g}}) {
      foreach my $og (keys %{$homols{$g}->{$spe}}) {
        print $outfh "Ortholog \"$og\" \"$spe\" From_analysis ParaSite-Compara\n";
      }
    }
    print $outfh "\n";
  }
}


###################################################
sub get_evidence{
  my($x)=@_;

  if($x->dbname =~/^Uniprot\//){
    return 'Accession_evidence "UniProt" "'.$x->primary_id.'"'
  } elsif($x->dbname eq 'ZFIN_ID'){
    return 'Accession_evidence "ZFIN" "'.$x->primary_id.'"'
  } elsif($x->dbname eq 'SGD_GENE'){
    return 'Accession_evidence "SGD" "'.$x->display_id.'"'
  } elsif($x->dbname eq 'MGI'){
    return'Accession_evidence "MGI" "'.$x->display_id.'"'
  } elsif($x->dbname eq 'FlyBaseCGID_gene'){
    return 'Accession_evidence "FLYBASE" "'.$x->display_id.'"'
  } elsif($x->dbname eq 'flybase_gene_id'){
    return 'Accession_evidence "FLYBASE" "'.$x->display_id.'"'
  } elsif($x->dbname eq 'HGNC'){
    return 'Accession_evidence "HGNC" "'.$x->display_id.'"'
  } elsif($x->dbname eq 'MIM_GENE'){
    return 'Accession_evidence "OMIM" "'.$x->primary_id.'"'
  } elsif($x->dbname eq 'MIM_MORBID'){
    return 'Accession_evidence "OMIM" "'.$x->primary_id.'"'
  } else{
    my $tmp=$x->dbname;
    $tmp=~s/;//;
    return 'Inferred_automatically "'.$x->dbname.'"'
  }
}


sub get_xref{
  my($x,$p)=@_;
  my $xs;
  if($x->dbname =~/^Uniprot\//){
    $xs = dbxrefP('UniProt','UniProt_AC',$x);
  } elsif($x->dbname eq 'ZFIN_ID'){
    $xs = dbxrefP('ZFIN','primary_acc',$x);
  } elsif($x->dbname eq 'SGD_GENE'){
    $xs = dbxrefS('SGD','SGD_acc',$x);
  } elsif($x->dbname eq 'MGI'){
    $xs = dbxrefS('MGI','MGI_acc',$x);
  } elsif($x->dbname eq 'FlyBaseCGID_gene'){
    $xs = dbxrefS('FLYBASE','FlyBase_gn',$x);
  } elsif($x->dbname eq 'flybase_gene_id'){
    $xs = dbxrefS('FLYBASE','FlyBase_ID',$x);
  } elsif($x->dbname eq 'HGNC'){
    $xs = dbxrefS('HGNC','symbol',$x);
  } elsif($x->dbname eq 'MIM_GENE'){
    $xs = dbxrefP('OMIM','gene',$x);
  } elsif($x->dbname eq 'MIM_MORBID'){
    $xs = dbxrefP('OMIM','disease',$x);
  } elsif($x->dbname eq 'EntrezGene') {
    $xs = dbxrefP('EntrezGene', "EntrezGene_id", $x);
  } elsif($x->dbname eq 'Uniprot_gn') {
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

