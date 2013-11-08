#!/usr/bin/env perl

use strict;
use Switch;
use Getopt::Long;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use WormBase2Ensembl;
#use WormFeature::Blat;
#use WormFeature::WublastX;
#use Clone2Acc;
#use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
#use Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


my ($estori, $species, $debug, $verbose,
    $dbhost, $dbname, $dbuser, $dbpass, $dbport,
    $gff_file, $seleno,
    );

&GetOptions('estorientations=s' => \$estori,
            'species=s'         => \$species,
            'host=s'            => \$dbhost,
            'user=s'            => \$dbuser, 
            'pass=s'            => \$dbpass,
            'port=s'            => \$dbport,
            'gff=s'             => \$gff_file,
            'dbname=s'          => \$dbname,
            'debug'             => \$debug,
            'verbose'           => \$verbose,
            'seleno'            => \$seleno,
    );



if (defined $estori and -e $estori) {
  $estori = &read_est_orientations($estori);
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $dbhost,
    -user   => $dbuser,
    -dbname => $dbname,
    -pass   => $dbpass,
    -port   => $dbport,
);

my @analysis_logic_names = @ARGV;

my $analysis_adaptor    = $db->get_AnalysisAdaptor();


my $wb2ens = WormBase2Ensembl->new(
  -species => $species,
  -debug   => $debug,
  -verbose => $verbose,
  -dbh     => $db);


foreach my $analysis_logic_name (@analysis_logic_names) {
  my $analysis = $analysis_adaptor->fetch_by_logic_name($analysis_logic_name);

  if (not $analysis and $analysis_logic_name ne 'wormbase_genes_gff3') {
    die "COULD NOT GET ANALYSIS ADAPTOR FOR $analysis_logic_name";
  }

  switch ($analysis_logic_name) {
    
    # simple features from GFF
    case [ 'operon','fosmid', 'rnai_pcr_product'] {
      my $features;
      switch ($analysis_logic_name) {
        case 'rnai_pcr_product'   { $features = $wb2ens->parse_pcr_products_gff( $gff_file, $analysis, undef, 'PCR_product' ) }
        case 'operon'             { $features = $wb2ens->parse_operons_gff( $gff_file, $analysis, 'operon' ) }
        case 'fosmid'	          { $features = $wb2ens->parse_simplefeatures_gff($gff_file, $analysis ,'Vancouver_fosmid')}
      }
      $wb2ens->write_simple_features( $features  ) if $features;
    }
    # BLASTX from GFF
    case ['slimswissprotx','slimtremblx','ipi_humanx','wormpepx','brigpepx','remaneix','brepepx','jappepx','ppapepx', 'gadflyx','yeastx'] {
      my $tag;
      switch ($analysis_logic_name){
        case 'slimswissprotx' {$tag='SW'}
        case 'slimtremblx'    {$tag='TR'}
        case 'ipi_humanx'     {$tag='ENSEMBL'}
        case 'gadflyx'        {$tag='FLYBASE'}
        case 'yeastx'         {$tag='SGD'}
        case 'wormpepx'       {$tag='WP'}
        case 'brigpepx'       {$tag='BP'}
        case 'remaneix'       {$tag='RP'}
        case 'jappepx'        {$tag='JA'}
        case 'brepepx'        {$tag='CN'}
        case 'ppapepx'        {$tag='PP'}
        case 'bmapepx'        {$tag='BM'}
        
      }
      
      my $hits = $wb2ens->parse_protein_align_features_gff($gff_file, $analysis, $tag);
      $wb2ens->write_protein_align_features($hits);
    }
    
    case ['celegans_mrna', 'cbriggsae_mrna', 'cremanei_mrna', 'cbrenneri_mrna', 'cremanei_mrna', 'cjaponica_mrna','ppacificus_mrna', 'bmalayi_mrna'] {
      my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'BLAT_mRNA_BEST', 'expressed_sequence_match', $estori);
      $wb2ens->write_dna_align_features($hits);
    }
    case ['celegans_est', 'cbriggsae_est', 'cremanei_est', 'cbrenneri_est', 'cremanei_est', 'cjaponica_est','ppacificus_est', 'bmalayi_est'] {
      my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'BLAT_EST_BEST', 'expressed_sequence_match', $estori);
      $wb2ens->write_dna_align_features($hits);
    }      
    case 'celegans_ost' {
      my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'BLAT_OST_BEST', 'expressed_sequence_match', $estori);
      $wb2ens->write_dna_align_features($hits);
    }
    case 'celegans_rst' {
      my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'BLAT_RST_BEST', 'expressed_sequence_match', $estori);
      $wb2ens->write_dna_align_features($hits);
    }
    
    case 'caenorhabditis_mrna' {
      my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'BLAT_Caen_mRNA_BEST', 'expressed_sequence_match', $estori);
      $wb2ens->write_dna_align_features($hits);
    }
    case 'caenorhabditis_est' {
      my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'BLAT_Caen_EST_BEST', 'expressed_sequence_match', $estori);
      $wb2ens->write_dna_align_features($hits);
    }
    case 'nembase_est' {
      my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'NEMBASE_cDNAs-BLAT', 'expressed_sequence_match', $estori);
      $wb2ens->write_dna_align_features($hits);
    }
    case 'washu_est' {
      my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'NEMATODE.NET_cDNAs-BLAT', 'expressed_sequence_match', $estori);
      $wb2ens->write_dna_align_features($hits);
    }
    case 'nematode_est' {
      my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'EMBL_nematode_cDNAs-BLAT', 'expressed_sequence_match', $estori);
      $wb2ens->write_dna_align_features($hits);
    }
    case 'wormbase_non_coding' {
      ## tRNA genes ##
      my $tRNA_scan_genes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis,'tRNA');
      $wb2ens->write_genes( $tRNA_scan_genes);
      
      ## rRNA-genes #
      my $rRNAgenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis, 'rRNA' );
      $wb2ens->write_genes( $rRNAgenes );
      
      ## generic ncRNA genes ##
      my $ncRNAgenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis,'ncRNA' );
      $wb2ens->write_genes( $ncRNAgenes );
      
      ## snRNA-genes
      my $snRNAgenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis,'snRNA' );
      $wb2ens->write_genes( $snRNAgenes );
      
      ## snoRNA-genes
      my $snoRNAgenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis, 'snoRNA' );
      $wb2ens->write_genes( $snoRNAgenes );
      
      ## scRNA-genes
      my $scRNAgenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis,'scRNA' );
      $wb2ens->write_genes( $scRNAgenes );
      
      ## miRNA-genes
      my $miRNAgenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis, 'miRNA_mature', 'miRNA');
      $wb2ens->write_genes( $miRNAgenes );
    }
    case 'wormbase_pseudogene' {
      ## pseudogenes
      my $pseudogenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis, 'Pseudogene', 'pseudogene');
      $wb2ens->write_genes( $pseudogenes );
    }
    case 'wormbase' {
      
      #parse out real worm genes
      my $genes = $wb2ens->parse_protein_coding_gff2( $gff_file, $analysis );
      $wb2ens->write_genes( $genes );
      
      #check translations
      my @genes = @{$db->get_GeneAdaptor->fetch_all_by_biotype('protein_coding')};
      &check_translations(\@genes);
    }
    case 'wormbase_genes_gff3' {
      my $cod_ana    = $db->get_AnalysisAdaptor->fetch_by_logic_name('wormbase');
      my $nc_ana     = $db->get_AnalysisAdaptor->fetch_by_logic_name('wormbase_non_coding');
      my $pseudo_ana = $db->get_AnalysisAdaptor->fetch_by_logic_name('wormbase_pseudogene');

      if (not $cod_ana or not $nc_ana or not $pseudo_ana) {
        die "Analyses 'wormbase', 'wormbase_non_coding' and 'wormbase_pseudogene' must be in the database to use this option\n";
      }
      my $genes = $wb2ens->parse_genes_gff3( $gff_file, $cod_ana, $nc_ana, $pseudo_ana, { WormBase => 1 });
      $verbose and print STDERR "Writing genes to database...\n";
      $wb2ens->write_genes( $genes );

      $verbose and print STDERR "Checking translations...\n";
      my @genes = @{$db->get_GeneAdaptor->fetch_all_by_biotype('protein_coding')};
      &check_translations(\@genes);
    }
    else {
      print "\n ANALYSIS NOT DEFINED ($analysis_logic_name)!\n";
      exit 1;
    }
  }
}


exit 0;


####################
sub read_est_orientations {
  my ($file) = @_;

  my ($data, $VAR1);

  open(my $fh, $file) or die "Could not open $file for reading\n";
  while(<$fh>) {
    chomp;
    $data .= $_;
  }

  eval $data;
  return $VAR1;
}



sub check_translations {
  my ($genes) = @_;

  foreach my $gene (@$genes) {
    my @non_translate;
    foreach my $tran (@{$gene->get_all_Transcripts}) {
      if ($tran->biotype eq 'protein_coding' and not $wb2ens->translation_check($tran)) {
        push @non_translate, $tran;
      }
    }
    
    if (@non_translate) {
      print "Transcripts from " . $gene->stable_id . " do not translate\n";
      if ($seleno) {
        print "-seleno option given, so will replace stops with selenocystein\n";
        foreach my $tran (@non_translate) {
          $wb2ens->translation_fix($tran);
        }
      }
    }
  }
}
