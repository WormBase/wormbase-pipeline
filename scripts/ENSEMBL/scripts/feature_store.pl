#!/usr/bin/env perl

use strict;
use Getopt::Long;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use WormBase2Ensembl;
use Bio::EnsEMBL::DBSQL::DBAdaptor;

my %protein_blast_logics = (
  'slimswissprotx' => 'SW',
  'slimtremblx'    => 'TR',
  'ipi_humanx'     => 'ENSEMBL',
  'gadflyx'        => 'FLYBASE',
  'yeastx'         => 'SGD',
  'wormpepx'       => 'WP',
  'brigpepx'       => 'BP',
  'remapepx'       => 'RP',
  'brepepx'        => 'CN',
  'jappepx'        => 'JA',
  'ppapepx'        => 'PP', 
  'brugpepx'       => 'BM',
  'ovolpepx'       => 'OV',
  'srapepx'        => 'SRP'
    );


my ($estori, $species, $debug, $verbose,
    $dbhost, $dbname, $dbuser, $dbpass, $dbport,
    $gff_file, $seleno, @gff3_sources,
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
            'selenocorrection'  => \$seleno,
            'gff3source=s@'     => \@gff3_sources,

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

  if ($analysis_logic_name eq 'operon' or
      $analysis_logic_name eq 'fosmid' or
      $analysis_logic_name eq 'rnai_pcr_product') {
    my $features;
    if ($analysis_logic_name eq 'operon') {
      $features = $wb2ens->parse_simplefeature_gff( $gff_file, $analysis, 'operon' );
    } elsif ($analysis_logic_name eq 'rnai_pcr_product') {
      $features = $wb2ens->parse_simplefeature_gff( $gff_file, $analysis, undef, 'PCR_product' );
    } elsif ($analysis_logic_name eq 'fosmid') {
      $features = $wb2ens->parse_simplefeature_gff( $gff_file, $analysis ,'Vancouver_fosmid');
    }
    $wb2ens->write_simple_features( $features  ) if $features;
  }
  # BLASTX from GFF
  elsif (exists $protein_blast_logics{$analysis_logic_name}) {
    my $tag = $protein_blast_logics{$analysis_logic_name};
    my $hits = $wb2ens->parse_protein_align_features_gff($gff_file, $analysis, $tag);
    $wb2ens->write_protein_align_features($hits);

  } elsif (grep { $analysis_logic_name eq $_ } ('celegans_mrna', 'cbriggsae_mrna', 'cremanei_mrna', 'cbrenneri_mrna', 'cremanei_mrna', 'cjaponica_mrna','ppacificus_mrna', 'bmalayi_mrna', 'ovolvulus_mrna','sratti_mrna')) {
    my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'BLAT_mRNA_BEST', 'expressed_sequence_match', $estori);
    $wb2ens->write_dna_align_features($hits);

  } elsif (grep { $analysis_logic_name eq $_ } ('celegans_est', 'cbriggsae_est', 'cremanei_est', 'cbrenneri_est', 'cremanei_est', 'cjaponica_est','ppacificus_est', 'bmalayi_est','ovolvulus_est','sratti_est')) {
    my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'BLAT_EST_BEST', 'expressed_sequence_match', $estori);
    $wb2ens->write_dna_align_features($hits);

  } elsif ($analysis_logic_name eq 'celegans_ost') {
    my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'BLAT_OST_BEST', 'expressed_sequence_match', $estori);
    $wb2ens->write_dna_align_features($hits);

  } elsif ($analysis_logic_name eq 'celegans_rst') {
    my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'BLAT_RST_BEST', 'expressed_sequence_match', $estori);
    $wb2ens->write_dna_align_features($hits);

  } elsif ($analysis_logic_name eq 'caenorhabditis_mrna') {
    my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'BLAT_Caen_mRNA_BEST', 'expressed_sequence_match', $estori);
    $wb2ens->write_dna_align_features($hits);

  } elsif ($analysis_logic_name eq 'caenorhabditis_est') {
    my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'BLAT_Caen_EST_BEST', 'expressed_sequence_match', $estori);
    $wb2ens->write_dna_align_features($hits);

  } elsif ($analysis_logic_name eq 'nembase_est') {
    my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'NEMBASE_cDNAs-BLAT', 'expressed_sequence_match', $estori);
    $wb2ens->write_dna_align_features($hits);

  } elsif ($analysis_logic_name eq 'washu_est') {
    my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'NEMATODE.NET_cDNAs-BLAT', 'expressed_sequence_match', $estori);
    $wb2ens->write_dna_align_features($hits);

  } elsif ($analysis_logic_name eq 'nematode_est') {
    my $hits = $wb2ens->parse_dna_align_features_gff($gff_file, $analysis, 'EMBL_nematode_cDNAs-BLAT', 'expressed_sequence_match', $estori);
    $wb2ens->write_dna_align_features($hits);

  } 
  # RepeatMask features
  elsif ($analysis_logic_name eq 'repeatmask') {
    my $rep_feats = $wb2ens->parse_repeatfeatures_gff( $gff_file, $analysis, 'RepeatMasker', 'similarity' );
    $verbose and printf STDERR "Read %d repeat features from GFF\n", scalar(@$rep_feats);
    $wb2ens->write_repeat_features($rep_feats);
  }
  # genes
  elsif ($analysis_logic_name eq 'wormbase_non_coding') {
    my @nc_genes;
    
    ## tRNA genes ##
    my $tRNA_scan_genes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis,'tRNA');
    $verbose and printf STDERR "Read %d tRNA genes from GFF\n", scalar(@$tRNA_scan_genes);
    push @nc_genes, @$tRNA_scan_genes;
    
    ## rRNA-genes #
    my $rRNAgenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis, 'rRNA' );
    $verbose and printf STDERR "Read %d rRNA genes from GFF\n", scalar(@$rRNAgenes);
    push @nc_genes, @$rRNAgenes;      
    
    ## generic ncRNA genes ##
    my $ncRNAgenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis,'ncRNA' );
    $verbose and printf STDERR "Read %d ncRNA genes from GFF\n", scalar(@$ncRNAgenes);
    push @nc_genes, @$ncRNAgenes;
    
    ## snRNA-genes
    my $snRNAgenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis,'snRNA' );
    $verbose and printf STDERR "Read %d snRNA genes from GFF\n", scalar(@$snRNAgenes);
    push @nc_genes, @$snRNAgenes;
    
    ## snoRNA-genes
    my $snoRNAgenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis, 'snoRNA' );
    $verbose and printf STDERR "Read %d snoRNA genes from GFF\n", scalar(@$snoRNAgenes);
    push @nc_genes, @$snoRNAgenes;
      
    ## scRNA-genes
    my $scRNAgenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis,'scRNA' );
    $verbose and printf STDERR "Read %d scRNA genes from GFF\n", scalar(@$scRNAgenes);
    push @nc_genes, @$scRNAgenes;
    
    ## miRNA-genes
    my $miRNAgenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis, 'miRNA_mature', 'miRNA');
    $verbose and printf STDERR "Read %d miRNA genes from GFF\n", scalar(@$miRNAgenes);
      push @nc_genes, @$miRNAgenes;
        
    ## lincRNA-genes
    my $lincRNAgenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis, 'lincRNA');
    $verbose and printf STDERR "Read %d lincRNA genes from GFF\n", scalar(@$lincRNAgenes);
      push @nc_genes, @$lincRNAgenes;
        
    ## asRNA-genes
    my $asRNAgenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis, 'asRNA');
    $verbose and printf STDERR "Read %d asRNA genes from GFF\n", scalar(@$asRNAgenes);
    push @nc_genes, @$asRNAgenes;
        
    ## piRNA-genes
    my $piRNAgenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis, 'piRNA');
    $verbose and printf STDERR "Read %d piRNA genes from GFF\n", scalar(@$piRNAgenes);
    push @nc_genes, @$piRNAgenes;
    
    $verbose and printf STDERR "Read %d ncRNA genes in all. Writing...\n", scalar(@nc_genes); 
    $wb2ens->write_genes( \@nc_genes, 1 );

  } elsif ($analysis_logic_name eq 'wormbase_pseudogene') {
    ## pseudogenes
    my $pseudogenes = $wb2ens->parse_non_coding_genes_gff2( $gff_file, $analysis, 'Pseudogene', 'pseudogene');
    $verbose and printf STDERR "Read %d pseudogenes. Writing...\n", scalar(@$pseudogenes); 
    $wb2ens->write_genes( $pseudogenes, 1);

  } elsif ($analysis_logic_name eq 'wormbase') {      
    #parse out real worm genes
    my $genes = $wb2ens->parse_protein_coding_gff2( $gff_file, $analysis );
    $verbose and printf STDERR "Read %d protein-coding genes. Writing...\n", scalar(@$genes); 
    $wb2ens->write_genes( $genes, 1 );
    
    #check translations
    $verbose and print STDERR "Checking translations...\n";
    my @genes = @{$db->get_GeneAdaptor->fetch_all_by_biotype('protein_coding')};
    &check_translations(\@genes);

  } elsif ($analysis_logic_name eq 'wormbase_genes_gff3') {
    my $cod_ana    = $db->get_AnalysisAdaptor->fetch_by_logic_name('wormbase');
    my $nc_ana     = $db->get_AnalysisAdaptor->fetch_by_logic_name('wormbase_non_coding');
    my $pseudo_ana = $db->get_AnalysisAdaptor->fetch_by_logic_name('wormbase_pseudogene');
    
    if (not $cod_ana or not $nc_ana or not $pseudo_ana) {
      die "Analyses 'wormbase', 'wormbase_non_coding' and 'wormbase_pseudogene' must be in the database to use this option\n";
    }

    my $source_hash;
    if (@gff3_sources) {
      $source_hash = {};
      map { $source_hash->{$_} = 1 } @gff3_sources;
    }

    my $genes = $wb2ens->parse_genes_gff3( $gff_file, $cod_ana, $nc_ana, $pseudo_ana, $source_hash );
    $verbose and printf STDERR "Parsed %d genes from GFF3. Writing genes to database...\n", scalar(@$genes);      
    $wb2ens->write_genes( $genes, 1 );
    
    $verbose and print STDERR "Checking translations...\n";
    my @genes = @{$db->get_GeneAdaptor->fetch_all_by_biotype('protein_coding')};
    &check_translations(\@genes);
  }
  else {
    print "\n ANALYSIS NOT DEFINED ($analysis_logic_name)!\n";
    exit 1;
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
