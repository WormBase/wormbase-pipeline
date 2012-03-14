#!/usr/local/ensembl/bin/perl  -w

use strict;
use Switch;
use Getopt::Long;

use WormBase;
use WormBaseConf;
use WormFeature::Blat;
use WormFeature::WublastX;
#use Clone2Acc;
#use Bio::EnsEMBL::Pipeline::SeqFetcher::Pfetch;
#use Bio::EnsEMBL::Pipeline::SeqFetcher::OBDAIndexSeqFetcher;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


my ($estori);
&GetOptions('estoritentiations=s' => \$estori);

if (defined $estori and -e $estori) {
  $estori = &read_est_orientations($estori);
}

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $WB_DBHOST,
    -user   => $WB_DBUSER,
    -dbname => $WB_DBNAME,
    -pass   => $WB_DBPASS,
    -port   => $WB_DBPORT,
);

my $analysis_logic_name = shift;
my $analysis_adaptor    = $db->get_AnalysisAdaptor();
my $analysis            = $analysis_adaptor->fetch_by_logic_name($analysis_logic_name);

if ( !$analysis && !($analysis_logic_name eq 'selenocystein')){
    die "COULD NOT GET ANALYSIS ADAPTOR FOR $analysis_logic_name";
}

foreach my $chromosome_info ( @{$WB_CHR_INFO} ) {

  if ( $chromosome_info->{'chr_name'} && 
       $chromosome_info->{'gff_file'} ) {
    
    my $gff_file = $chromosome_info->{gff_file};
    my $chr_name = $chromosome_info->{chr_name};
    my $agp_file = $chromosome_info->{agp_file};

    if (not -e $gff_file) {
      $gff_file = join("/", $WB_workDIR, $chromosome_info->{gff_file});
      die "Could not find $gff_file" if not -e $gff_file;
    }
    
    print "Handling $chr_name with file  $gff_file\n" if ($WB_DEBUG);
    
    my $chr = $db->get_SliceAdaptor->fetch_by_region( 'chromosome', $chr_name );
    
    #MtDNA does not have Expr_profile, RNAi or Operon so may cause hassles for these
    next if ( $chr_name eq 'MtDNA' and
              ( $analysis_logic_name eq 'expression_profile' or 
                $analysis_logic_name eq 'rnai' or 
                $analysis_logic_name eq 'operon' or 
                $analysis_logic_name eq 'fosmid'));
    
    #read the features from gff file
    switch ($analysis_logic_name) {
      
      # simple features from GFF
      case [ 'expression_profile', 'rnai', 'operon','fosmid'] {
        my $features;
        switch ($analysis_logic_name) {
          case 'expression_profile' { $features = parse_expr( $gff_file,    $chr, $analysis ) }
          case 'rnai'               { $features = parse_rnai( $gff_file,    $chr, $analysis ) }
          case 'operon'             { $features = parse_operons( $gff_file, $chr, $analysis ) }
          case 'fosmid'	      { $features = parse_simplefeature($gff_file, $chr, $analysis ,'Vancouver_fosmid')}
        }
        print "have " . scalar @$features . " features ($analysis_logic_name).\n" if ($WB_DEBUG);
        
        #save the features to db
        write_simple_features( $features, $db ) if $features;
      }
      # BLASTX from GFF
      case ['Swissprot','TrEMBL','Human-IPR','WormPep','FlyBase','BrigPep','RemaneiPep','SGD'] {
        my $tag;
        switch ($analysis_logic_name){
          case 'Swissprot' {$tag='SW'}
          case 'TrEMBL' {$tag='TR'}
          case 'Human-IPR'{$tag='ENSEMBL'}
          case 'WormPep'{$tag='WP'}
          case 'FlyBase'{$tag='FLYBASE'}
          case 'BrigPep'{$tag='BP'}
          case 'RemaneiPep'{$tag='RP'}
          case 'SGD'{$tag='SGD'}
        }
        
        my $wublastx = WublastX->new( $chr, $gff_file, $analysis ,$tag);
        my $hits=  scalar( @{ $wublastx->{hits} } );
        print "has $hits $analysis_logic_name hits\n" if $WB_DEBUG;
        $wublastx->save($db) if $hits >0;
      }
      
      case 'celegans_mrna' {
        my $blat = Blat->new( $chr, $gff_file, $analysis ,{-feature => 'BLAT_mRNA_BEST', -source => 'cDNA_match'}, $estori);
        my $hits=  scalar( @{ $blat->{hits} } );
        print "has $hits  BLAT hits\n" if $WB_DEBUG;
        $blat->save($db) if $hits >0;
      }
      case 'celegans_est' {
        my $blat = Blat->new( $chr, $gff_file, $analysis ,{-feature => 'BLAT_EST_BEST', -source => 'EST_match'}, $estori);
        my $hits=  scalar( @{ $blat->{hits} } );
        print "has $hits  BLAT hits\n" if $WB_DEBUG;
        $blat->save($db) if $hits >0;
      }      
      case 'ost' {
        my $blat = Blat->new( $chr, 
                              $gff_file, 
                              $analysis ,
                              {-feature => 'BLAT_OST_BEST', -source => 'expressed_sequence_match'}, 
                              $estori);
        my $hits=  scalar( @{ $blat->{hits} } );
        print "has $hits  BLAT hits\n" if $WB_DEBUG;
        $blat->save($db) if $hits >0;
      }
      case 'rst' {
        my $blat = Blat->new( $chr, 
                              $gff_file, 
                              $analysis ,
                              {-feature => 'BLAT_RST_BEST', -source => 'expressed_sequence_match' }, 
                              $estori);
        my $hits=  scalar( @{ $blat->{hits} } );
        print "has $hits  BLAT hits\n" if $WB_DEBUG;
        $blat->save($db) if $hits >0;
      }

      case 'caenorhabditis_mrna' {
        my $blat = Blat->new( $chr, 
                              $gff_file, 
                              $analysis ,
                              {-feature => 'BLAT_Caen_mRNA_BEST', -source => 'expressed_sequence_match'},
                              $estori);
        my $hits=  scalar( @{ $blat->{hits} } );
        print "has $hits  BLAT hits\n" if $WB_DEBUG;
        $blat->save($db) if $hits >0;
      }
      case 'caenorhabditis_est' {
        my $blat = Blat->new( $chr, 
                              $gff_file, 
                              $analysis ,
                              {-feature => 'BLAT_Caen_EST_BEST', -source => 'expressed_sequence_match' },
                              $estori);
        my $hits=  scalar( @{ $blat->{hits} } );
        print "has $hits  BLAT hits\n" if $WB_DEBUG;
        $blat->save($db) if $hits >0;
      }
      case 'nembase_contig' {
        my $blat = Blat->new( $chr, 
                              $gff_file, 
                              $analysis ,
                              {-feature => 'BLAT_NEMBASE', -source => 'translated_nucleotide_match' },
                              $estori);
        my $hits=  scalar( @{ $blat->{hits} } );
        print "has $hits  BLAT hits\n" if $WB_DEBUG;
        $blat->save($db) if $hits >0;
      }
      case 'washu_contig' {
        my $blat = Blat->new( $chr, 
                              $gff_file, 
                              $analysis ,
                              {-feature => 'BLAT_WASHU', -source => 'translated_nucleotide_match' },
                              $estori);
        my $hits=  scalar( @{ $blat->{hits} } );
        print "has $hits  BLAT hits\n" if $WB_DEBUG;
        $blat->save($db) if $hits >0;
      }
      case 'other_est' {
        my $blat = Blat->new( $chr, 
                              $gff_file, 
                              $analysis ,
                              {-feature => 'BLAT_NEMATODE', -source => 'translated_nucleotide_match' },
                              $estori);
        my $hits=  scalar( @{ $blat->{hits} } );
        print "has $hits  BLAT hits\n" if $WB_DEBUG;
        $blat->save($db) if $hits >0;
      }
      case 'wormbase_non_coding' {
        ## tRNA genes ##
        my $tRNA_scan_genes = parse_pseudo_files( $gff_file, $chr, $analysis,'tRNAscan-SE-1.23', 'tRNA' );
        print "have " . scalar @$tRNA_scan_genes . " genomic tRNA genes.\n" if ($WB_DEBUG);
        &write_genes( $tRNA_scan_genes, $db );
        
        ## rRNA-genes #
        my $rRNAgenes = parse_pseudo_files( $gff_file, $chr, $analysis, 'rRNA' );
        print "have " . scalar @$rRNAgenes . " rRNA genes.\n" if ($WB_DEBUG);        
        &write_genes( $rRNAgenes, $db );

        ## generic ncRNA genes ##
        my $ncRNAgenes = parse_pseudo_files( $gff_file, $chr, $analysis,'ncRNA' );
        print "have " . scalar @$ncRNAgenes . " ncRNA genes.\n" if ($WB_DEBUG);
        &write_genes( $ncRNAgenes, $db );
        
        ## snRNA-genes
        my $snRNAgenes = parse_pseudo_files( $gff_file, $chr, $analysis,'snRNA_mature_transcript', 'snRNA' );
        print "have " . scalar @$snRNAgenes . " snRNA genes.\n" if ($WB_DEBUG);        
        &write_genes( $snRNAgenes, $db );

        ## snoRNA-genes
        my $snoRNAgenes = parse_pseudo_files( $gff_file, $chr, $analysis,'snoRNA_mature_transcript', 'snoRNA' );
        print "have " . scalar @$snoRNAgenes . " snoRNA genes.\n" if ($WB_DEBUG);
        &write_genes( $snoRNAgenes, $db );

        ## snlRNA-genes
        my $snlRNAgenes = parse_pseudo_files( $gff_file, $chr, $analysis,'snlRNA' );
        print "have " . scalar @$snlRNAgenes . " snlRNA genes.\n" if ($WB_DEBUG);
        &write_genes( $snlRNAgenes, $db );

        ## miRNA-genes
        my $miRNAgenes = parse_pseudo_files( $gff_file, $chr, $analysis,'curated_miRNA', 'miRNA' );
        print "have " . scalar @$miRNAgenes . " miRNA genes.\n" if ($WB_DEBUG);
        &write_genes( $miRNAgenes, $db );
      }
      case 'wormbase_pseudogene' {
        ## pseudogenes
        my $pseudogenes = parse_pseudo_files( $gff_file, $chr, $analysis, 'Pseudogene', 'pseudogene');
        print "have " . scalar @$pseudogenes . " pseudogenes / tRNA genes.\n" if ($WB_DEBUG);
        write_genes( $pseudogenes, $db, 1 );
      }
      case 'wormbase' {
        
        #parse out real worm genes
        my $genes = parse_gff( $gff_file, $chr, $analysis );
        print "have " . scalar @$genes . " genes.\n" if ($WB_DEBUG);
        
        #store genes
        &write_genes( $genes, $db );
        
        #check translations
        my @genes = @{ $chr->get_all_Genes };
        open( TRANSLATE, "+>>" . $WB_NON_TRANSLATE ) or die "couldn't open " . $WB_NON_TRANSLATE . " $!";
        TRANSLATION: foreach my $gene (@genes) {
          my $translation = &translation_check($gene);
          if ($translation) {
            next TRANSLATION;
          }
          else {
            print TRANSLATE $gene->stable_id . " from $chr_name does not translate\n";
            next TRANSLATION;
          }
        }
        close(TRANSLATE);
      }
      case 'selenocystein' {
        #check translations

        my @genes = grep {$_->biotype eq 'protein_coding'} @{ $chr->get_all_Genes };

        print "GENES = ". scalar(@genes), "\n";
        open( TRANSLATE, "+>>" . $WB_NON_TRANSLATE ) or die "couldn't open " . $WB_NON_TRANSLATE . " $!";
        foreach my $gene (@genes) {
          my $translation = &translation_check($gene,$db);
          next if ($translation);
          print TRANSLATE $gene->stable_id . " from $chr_name does not translate\n";
        }
        close(TRANSLATE);
        
      }
      
      else {
        print "\n ANALYSIS NOT DEFINED ($analysis_logic_name)!\n";
        exit 1;
      }
    }
  }
  else { 
    print "skipping chromosome.\n" if ($WB_DEBUG); 
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
