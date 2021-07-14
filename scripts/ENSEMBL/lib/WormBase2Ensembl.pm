#
# You may distribute this module under the same terms as perl itself
#
# POD documentation - main docs before the code
# mod fsk

=pod 

=head1 NAME

WormBase

=head1 SYNOPSIS

Class for converting WormBase GFF (v2, mostly) to Ensembl

=cut

package WormBase2Ensembl;

use strict;

use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::SimpleFeature;
use Bio::EnsEMBL::DnaDnaAlignFeature;
use Bio::EnsEMBL::DnaPepAlignFeature;
use Bio::EnsEMBL::FeaturePair;
use Bio::EnsEMBL::Analysis::Tools::FeatureFactory;

#use Bio::EnsEMBL::Analysis::Tools::GeneBuildUtils::ExonUtils;


########################
sub new {
  my $class = shift;

  my %args = @_;

  my $self = { };
  bless $self, $class;

  if (exists $args{"-species"}) {
    $self->species($args{"-species"});
  }
  if (exists $args{"-debug"}) {
    $self->debug($args{"-debug"});
  }
  if (exists $args{"-recognise_sources"}) {
    $self->recognise_sources($args{"-recognise_sources"});
  }
  if (exists $args{"-verbose"}) {
    $self->verbose($args{"-verbose"});
  }
  if (exists $args{"-dbh"}) {
    $self->database_handle($args{"-dbh"});
  }
  if (exists $args{"-slices"}) {
    $self->slices($args{"-slices"});
  }
  if (exists $args{"-ignoregffphases"}) {
    $self->ignore_gff_phases($args{"-ignoregffphases"});
  }

  return $self;
}


sub parse_protein_coding_gff2 {
  my ( $self, $file, $analysis ) = @_;
  
  my $fh = $self->_open_file($file); 
  my $genes = $self->parse_protein_coding_gff2_fh($fh, $analysis);
  
  return $genes;
}


sub parse_protein_coding_gff2_fh {
  my ( $self, $fh, $analysis ) = @_;
    
  my ( $coding_isoforms, $non_coding_isoforms, $five_prime, $three_prime, $parent_seqs ) = 
      $self->_process_protein_coding_file( $fh );
  
  my ( $processed_transcripts, $five_start, $three_end, $trans_start_exon, $trans_end_exon ) =
      $self->_generate_protein_coding_transcripts( $coding_isoforms, $analysis, $five_prime, $three_prime, $parent_seqs );
    
  my $genes1 = $self->_create_protein_coding_transcripts( $processed_transcripts, $five_start, $three_end, $trans_start_exon, $trans_end_exon );
  
  my $nc_trans = $self->_process_pseudo_transcripts($non_coding_isoforms, $analysis);    
  my $genes2   = $self->_create_pseudo_transcripts($nc_trans, 'ncRNA');

  my @genes;
  foreach my $gene_id ( keys(%$genes1) ) {
    my @all_transcripts = @{$genes1->{$gene_id}};
    if (exists $genes2->{$gene_id}) {
      push @all_transcripts, @{$genes2->{$gene_id}};
    }
    my $gene    = $self->_create_gene( \@all_transcripts, $gene_id, 'protein_coding' );
    push( @genes, $gene );
  }
  
  return \@genes;
}

sub parse_genes_gff3 {
  my ($self, $file, $coding_ana, $nc_ana, $pseudo_ana, $source_hash) = @_;

  my $fh = $self->_open_file($file);
  my $genes = $self->parse_genes_gff3_fh($fh,  $coding_ana, $nc_ana, $pseudo_ana, $source_hash);

  return $genes;
}

sub parse_genes_gff3_fh {
  my ($self, $fh, $coding_ana, $nc_ana, $pseudo_ana, $source_hash) = @_;

  $self->verbose and print STDERR "Reading GFF stream for genes...\n";

  my (%genes, %transcripts, %gene_names);

  while(<$fh>) {
    chomp;
    next unless $_;
    next if /^\#/;
    my @l = split(/\t+/, $_);

    next if defined $source_hash and not exists $source_hash->{$l[1]};
    
    die "pseudogene: should instead be gene + pseudogenic_transcript"
      if $l[2] eq 'pseudogene';
    die "nontranslating_CDS should be gene + nontranslating_transcript"
      if $l[2] eq 'nontranslating_CDS';
    next if (
      $l[2] ne 'gene' and 
      $l[2] ne 'rRNA' and 
      $l[2] ne 'tRNA' and 
      $l[2] ne 'miRNA' and 
      $l[2] ne 'pre_miRNA' and 
      $l[2] ne 'scRNA' and
      $l[2] ne 'snoRNA' and 
      $l[2] ne 'snRNA' and 
      $l[2] ne 'ncRNA' and 
      $l[2] ne 'lincRNA' and 
      $l[2] ne 'piRNA' and 
      $l[2] ne 'antisense_RNA' and 
      $l[2] ne 'circular_RNA' and 
      $l[2] ne 'mRNA' and 
      $l[2] ne 'nc_primary_transcript' and
      $l[2] ne 'miRNA_primary_transcript' and
      $l[2] ne 'protein_coding_primary_transcript' and  
      $l[2] ne 'pseudogenic_transcript' and 
      $l[2] ne 'pseudogenic_rRNA' and
      $l[2] ne 'pseudogenic_tRNA' and
      $l[2] ne 'nontranslating_transcript' and 
      $l[2] ne 'CDS' and 
      $l[2] ne 'exon');

    my ($id, $name,  %parents);

    foreach my $el (split(/;/, $l[8])) {
      my ($k, $v) = $el =~ /(\S+)=(.+)/;
      if ($k eq 'ID') {
        ($id) = $v =~ /^(\S+)/; 
      } elsif ($k eq 'Name') {
        ($name) = $v =~ /^(\S+)/;
      } elsif ($k eq 'Parent') {
        foreach my $p (split(/,/, $v)) {
          my ($par) = $p =~ /^(\S+)/; 
          $parents{$par} = 1;
        }
      } 
    }

    if ($l[2] eq 'gene') {
      die "Every gene feature must have the Name attribute defined\n" if not $name;
      $gene_names{$id} = $name;
      next;
    }
    
    if ($l[2] eq 'exon') {
      # parent is mRNA 
      foreach my $parent (keys %parents) {
        push @{$transcripts{$parent}->{exons}}, {
          seq       => $l[0],
          start     => $l[3],
          end       => $l[4],
          strand    => $l[6],
          phase     => -1,
          end_phase => -1,
        };
      }
    } elsif ($l[2] eq 'CDS') {
      foreach my $parent (keys %parents) {
        push @{$transcripts{$parent}->{cds}}, {
          seq        => $l[0],
          start      => $l[3],
          end        => $l[4],
          strand     => $l[6],
          gff_phase  => $l[7],
        };
      }
    } elsif ($l[2] eq 'mRNA' or 
             $l[2] eq 'tRNA' or
             $l[2] eq 'rRNA' or
             $l[2] eq 'rRNA' or
             $l[2] eq 'tRNA' or
             $l[2] eq 'miRNA' or
             $l[2] eq 'pre_miRNA' or
             $l[2] eq 'miRNA_primary_transcript' or
             $l[2] eq 'scRNA' or
             $l[2] eq 'snoRNA' or
             $l[2] eq 'ncRNA' or
             $l[2] eq 'snRNA' or
             $l[2] eq 'lincRNA' or
             $l[2] eq 'piRNA' or
             $l[2] eq 'antisense_RNA' or
             $l[2] eq 'circular_RNA' or
             $l[2] eq 'pseudogenic_transcript' or
	     $l[2] eq 'nontranslating_transcript' or
             $l[2] eq 'pseudogenic_rRNA' or
             $l[2] eq 'pseudogenic_tRNA' or
             $l[2] eq 'nc_primary_transcript' or 
             $l[2] eq 'protein_coding_primary_transcript') {
      $transcripts{$id}->{source} = $l[1];
      $transcripts{$id}->{type}   = $l[2];
      $transcripts{$id}->{seq}    = $l[0];
      $transcripts{$id}->{start}  = $l[3];
      $transcripts{$id}->{end}    = $l[4];
      $transcripts{$id}->{strand} = $l[6];
      $transcripts{$id}->{name}   = $name;

      foreach my $parent (keys %parents) {
        push @{$genes{$parent}}, $id;
      }
    }
  }

  $self->verbose and print STDERR "Finished reading. Making Gene objects...\n";
   
  my @all_ens_genes;

  GENE: foreach my $gid (sort keys %genes) {
    my $sid = exists $gene_names{$gid} ? $gene_names{$gid} : $gid;
    my @tids = sort @{$genes{$gid}};

    my $gene = Bio::EnsEMBL::Gene->new(
      -stable_id => $sid,
      -source => 'WormBase',
        );
    my (%gene_biotypes);
    my $tcount = 1;
    
    # sanity check: all transcripts have the same strand
    my @strands_unique = do { my %seen; grep { !$seen{$transcripts{$_}->{strand}}++ } @tids }; 
    if(scalar(@strands_unique) > 1 ) {
       die "Bailing on gene $gid - transcripts not all on the same strand";
    }
    TRAN: foreach my $tid (sort @tids) {

      my $tran = $transcripts{$tid};
      my $gff_type = $tran->{type};
      my $strand = $tran->{strand};
      my ($tsid, $gff_source);
      if ($self->recognise_sources){
        $gff_source = $tran->{source};
      }
      else{
        $gff_source = "WormBase"; # to maintain previous behaviour
      }
      if ($tran->{name}) {
        $tsid = $tran->{name};
      } elsif (scalar(@tids) > 1) {
        $tsid = sprintf("%s.%d", $sid, $tcount++);
      } else {
        $tsid = $sid;
      }

      if (not exists $tran->{exons}) {
        # For some single-exon features, e.g. ncRNAs, exons are not compulsory in the GFF3. 
        # Add that exon dynamically here then
        push @{$tran->{exons}}, {
          seq       => $tran->{seq},
          start     => $tran->{start},
          end       => $tran->{end},
          strand    => $tran->{strand},
          phase     => -1,
          end_phase => -1,
        };
      }

      my @exons = sort { $a->{start} <=> $b->{start} } @{$tran->{exons}};

      # sanity check: exons have same strand as transcript
      foreach my $e (@exons) {
        if ($e->{strand} ne $strand) {
          die "Bailing on transcript $tid - at least one exon has different strand\n";
        }
      }

      # sanity check: exons do not overlap
      for(my $i=0; $i < scalar(@exons) - 1; $i++) {
        if ($exons[$i]->{end} >= $exons[$i+1]->{start}) {
          die "Bailing on transcript $tid - exons overlap\n";
        }
      }


      my $slice = $self->slices->{$exons[0]->{seq}};
      die "Could not find slice for $exons[0]->{seq}\n" if not defined $slice;

      @exons = reverse @exons if $strand eq '-';

      if (exists $tran->{cds}) {
        my @cds = sort { $a->{start} <=> $b->{start} } @{$tran->{cds}};

        # sanity check: CDS segments do not overlap
        for(my $i=0; $i < scalar(@cds) - 1; $i++) {
          if ($cds[$i]->{end} >= $cds[$i+1]->{start}) {
            die "Bailing on transcript $tid - CDS segments overlap\n";
          }
        }

        #
        # match up the CDS segments with exons
        #
        @cds = reverse @cds if $strand eq '-';

        foreach my $exon (@exons) {
          my ($matching_cds, @others) = grep { $_->{start} <= $exon->{end} and $_->{start} >= $exon->{start} } @cds;
          
          if (@others) {
            die "Bailing on transcript $tid: multiple matching CDSs for a single exon segment\n";
          }
          if (defined $matching_cds) {
            $exon->{cds_seg} = $matching_cds;
          } 
        }

        # Rarely, UTRs are annotated on CDSs with non-zero start/end phases, which does not
        # make sense. In this case, we clip off the 5' UTR
        my $phase_5prime;
        map { $phase_5prime = $_->{cds_seg}->{gff_phase} if exists $_->{cds_seg} and not defined $phase_5prime } @exons;

        if ($phase_5prime != 0) {
          my @keep_exons;
          while(my $exon = shift @exons) {
            # skip 5' exons that are wholly UTR
            if (not exists $exon->{cds_seg}) {
              push @keep_exons, $exon if @keep_exons;
              next;
            }

            if ($strand eq '-' and $exon->{end} != $exon->{cds_seg}->{end}) {
              # clip the exon back to the CDS end
              $exon->{end} = $exon->{cds_seg}->{end};
            } elsif ($strand eq '+' and $exon->{start} != $exon->{cds_seg}->{start}) {
              # clip the exon forward to the CDS start
              $exon->{start} = $exon->{cds_seg}->{start};
            }
            push @keep_exons, $exon;
          }

          @exons = @keep_exons;
        }
      }

      #
      # sort out phases.
      # 
      for(my $i=0; $i < @exons; $i++) {
        my $exon = $exons[$i];
        
        if (exists $exon->{cds_seg}) {

          if (($strand eq '+' and $exon->{start} == $exon->{cds_seg}->{start}) or
              ($strand eq '-' and $exon->{end} == $exon->{cds_seg}->{end})) {
            # phase should always match end_phase of previous exon, if there is one
            if ($i > 0 and exists($exons[$i-1]->{end_phase})) {
              $exon->{phase} = $exons[$i-1]->{end_phase};
            } elsif ($self->ignore_gff_phases) {
              # we have been told that the GFF phases are unreliable; therefore assume a start phase of 0
              $exon->{phase} = 0;
            } else {
              $exon->{phase} = (3 - $exon->{cds_seg}->{gff_phase}) % 3;
            }
          }
          
          if (($strand eq '+' and $exon->{end} == $exon->{cds_seg}->{end}) or
              ($strand eq '-' and $exon->{start} == $exon->{cds_seg}->{start})) {
            my $cds_seg_len = $exon->{cds_seg}->{end} - $exon->{cds_seg}->{start} + 1;
            if ($exon->{phase} > 0) {
              $cds_seg_len += $exon->{phase};
            }
            $exon->{end_phase} = $cds_seg_len % 3;
          }

          if ($strand eq '+') {
            $exon->{cds_exon_start} = $exon->{cds_seg}->{start} - $exon->{start} + 1;
            $exon->{cds_exon_end}   = $exon->{cds_seg}->{end} - $exon->{start} + 1; 
          } else {
            $exon->{cds_exon_start} = $exon->{end} - $exon->{cds_seg}->{end} + 1;
            $exon->{cds_exon_end}   = $exon->{end} - $exon->{cds_seg}->{start} + 1;
          }
        }
      }
      
      my $transcript = Bio::EnsEMBL::Transcript->new(
        -stable_id => $tsid,
        -source => $gff_source,
          );
      $transcript->version(undef);


      my ($tr_st_ex, $tr_en_ex, $tr_st_off, $tr_en_off);
      my $e_count = 1;

      foreach my $ex (@exons) {
        
        my $ens_ex = Bio::EnsEMBL::Exon->new(
          -slice     => $slice,
          -start     => $ex->{start},
          -end       => $ex->{end},
          -strand    => ($ex->{strand} eq '+') ? 1 : -1,
          -phase     => $ex->{phase},
          -end_phase => $ex->{end_phase},
          -stable_id => sprintf("%s.e%d", $tsid, $e_count++),
            );
        $ens_ex->version(undef);
        $transcript->add_Exon($ens_ex);

        if ($ex->{cds_seg}) {
          if (not defined $tr_st_ex) {
            $tr_st_ex = $ens_ex;
            $tr_st_off = $ex->{cds_exon_start};
          }
          $tr_en_ex = $ens_ex;
          $tr_en_off = $ex->{cds_exon_end};
        }
      }
  
      if (defined $tr_st_ex and defined $tr_en_ex) {
        my $tr = Bio::EnsEMBL::Translation->new();
        $tr->start_Exon($tr_st_ex);
        $tr->end_Exon($tr_en_ex);
        $tr->start($tr_st_off);
        $tr->end($tr_en_off);
        $tr->stable_id($tsid);
        $tr->version(undef);

        $transcript->translation($tr);
        $transcript->biotype('protein_coding');
        $transcript->analysis($coding_ana);
        $gene_biotypes{protein_coding}++;
      } else {
        if ($gff_type eq 'nc_primary_transcript') {
          # non-coding transcript in a coding gene
          $transcript->analysis($coding_ana);
          $transcript->biotype('ncRNA');
          $gene_biotypes{ncRNA}++;
        } elsif ($gff_type =~ /pseudo/) {
          $transcript->analysis($pseudo_ana);
          my $biotype = "pseudogene";
          if ($gff_type =~ /pseudogenic_(\S+)/ and $1 ne 'transcript') {
            $biotype = "${1}_pseudogene";
          }
          $transcript->biotype($biotype);
          $gene_biotypes{pseudogene}++;
       	} elsif ($gff_type =~ /nontranslating_transcript/){
	  $transcript->analysis($coding_ana); # analysis description to match the protein coding genes.
	  $transcript->biotype('nontranslating_CDS');
	  $gene_biotypes{nontranslating_CDS}++;
	} elsif ( $gff_type eq 'mRNA'){
          # mRNA feature with no corresponding CDS. Barf
          die "Transcript $tid is an mRNA, but could not get a valid CDS for it. Aborting\n";
        } elsif ( $gff_type eq 'pre_miRNA' or
                  $gff_type eq 'miRNA' or
                  $gff_type eq 'miRNA_primary_transcript') {
          $transcript->analysis($nc_ana);
          $transcript->biotype($gff_type);
          $gene_biotypes{"miRNA"}++;
        } elsif ( $gff_type eq 'scRNA') {
          # not acknowledged as a biotype by Ensembl; change to default ncRNA
          $transcript->analysis($nc_ana);
          my $bt = "ncRNA";
          $transcript->biotype($bt);
          $gene_biotypes{$bt}++;
        } elsif ( $gff_type eq 'circular_RNA') {
          # not acknowledged as a biotype by Ensembl; change to default ncRNA
          $transcript->analysis($nc_ana);
          my $bt = "ncRNA";
          $transcript->biotype($bt);
          $gene_biotypes{$bt}++;
        } else {
          $transcript->analysis($nc_ana);
          my $bt = ($gff_type =~ /RNA/) ? $gff_type : 'ncRNA';
          $transcript->biotype($bt);
          $gene_biotypes{$bt}++;
        }
      }

      $gene->add_Transcript($transcript);
    }

    my ($first_bt, @other_bts) = sort keys %gene_biotypes;
    if (not @other_bts) {
      $gene->biotype($first_bt);
    } else {
      if (exists $gene_biotypes{protein_coding}) {
        $gene->biotype('protein_coding');
      } else {
        # default to generic nc-RNA
        $gene->biotype('ncRNA');
      }
    }

    # propagate analysis for first transcript up to the gene
    $gene->analysis( $gene->get_all_Transcripts->[0]->analysis );
    $gene->source( $gene->get_all_Transcripts->[0]->source );

    push @all_ens_genes, $gene;
  }

  $self->verbose and print STDERR "Returning ", scalar(@all_ens_genes),  " gene.\n";

  return \@all_ens_genes;
}


sub translation_check {
  my ( $self, $t ) = @_;
  
  my $pep = $t->translate->seq;
  if ( $pep =~ /\*/g ) {
    print STDERR "transcript " . $t->stable_id . " doesn't translate\n";
    if ($self->verbose) {
      print STDERR "translation start " . $t->translation->start . " end " . $t->translation->end . "\n";
      print STDERR "start exon coords " . $t->translation->start_Exon->start . " " . $t->translation->start_Exon->end . "\n";
      print STDERR "end exon coords " . $t->translation->end_Exon->start . " " . $t->translation->end_Exon->end . "\n";
      print STDERR "peptide " . $pep . "\n";
      
      $self->_display_exons( @{ $t->get_all_Exons } );
      $self->_display_non_translate($t);
    }
    return 0;
  } else {
    return 1;
  }
}


sub translation_fix {
  my ( $self, $t ) = @_;
  
  my $pep = $t->translate->seq;
  my $fixed = 0;
  if ( $pep =~ /\*/g ) {
    #add Selenocysteine to translation. There seems to be only on Selenoc. in our worm...
    my $pos = pos($pep);
    print STDERR "transcript "
        . $t->stable_id
        . " doesn't translate. Adding Selenocystein at position $pos.\n"
        . "Please beware of problems during further analysis steps because of this.\n";
    $self->_selenocysteine( $t, $pos );
    $fixed = 1;
  }
  return $fixed;
}

# Sometimes submitters give us partial genes - the models are incomplete because of scaffold gap
# We can rescue such genes by correcting the phase of the first exon.
# Also: if an exon of just one phase is wrong, and it gets fixed, that's good
sub phase_fix {
  my ($self, $g_ref) = @_;

  $self->verbose and print STDERR "Fixing phases...\n";

  foreach my $g (@$g_ref) {
    T:
    foreach my $t (@{$g->get_all_Transcripts}) {
      next T if $t->biotype ne 'protein_coding';
      
      my @exons = @{$t->get_all_Exons}; 
      
      next T if $exons[0]->phase < 0;

      next T if seq_ok($t);

      for my $i (0..$#exons){
        my $original_phase = $exons[$i]->phase;

        my %stops_by_phase;
        for my $phase (grep {$_ ne $original_phase} (0, 1, 2)) {
          $exons[$i]->phase($phase);
        
          if(seq_ok($t)){
              printf STDERR "Fixed transcript %s: changed phase of the exon %d phase from %d to %d\n", $t->stable_id, $i, $original_phase, $phase
                if $self->verbose;
              next T;
          }
        }
        $exons[$i]->phase($original_phase);

        my $original_start = $exons[$i]->start;
        for my $d (-1, 1, -2, 2){
          $exons[$i]->start($original_start+$d);
          if (seq_ok($t)){
             printf STDERR "Fixed transcript %s: changed start of the exon %d from %d to %d\n", $t->stable_id, $i, $original_start, $original_start+$d
               if $self->verbose;
             next T;
          }
        }
        $exons[$i]->start($original_start); 

        my $original_end = $exons[$i]->end;
        for my $d (-1, 1, -2, 2){
          $exons[$i]->end($original_end+$d);
          if (seq_ok($t)){
             printf STDERR "Fixed transcript %s: changed end of the exon %d from %d to %d\n", $t->stable_id, $i, $original_end, $original_end+$d
               if $self->verbose;
             next T;
          }
        }
        $exons[$i]->end($original_end);

      }
    }
  }
}
sub seq_ok {
  my ($t) = @_;
  my $seq = $t->translate->seq;
  my ($stops) = $seq =~ tr/\*/\*/;
  return not $stops;
}

sub parse_simplefeature_gff {
  my ($self, $file, $analysis, $source, $type) = @_;

  my $fh = $self->_open_file($file);
  my $feats = $self->parse_simplefeature_gff_fh($fh, $analysis, $source, $type);

  return $feats;
}


sub parse_simplefeature_gff_fh {
  my ($self, $fh, $analysis, $given_source, $given_type ) = @_;
    
  my @features;
  while (<$fh>) {
    next if /^\#/;
    my @l = split(/\t/, $_);
    my ($chr, $source, $type, $start, $end, $strand, $attr) = 
        @l[0,1,2,3,4,6,8];
    
    next if defined $given_source and $source ne $given_source;
    next if defined $given_type and $type ne $given_type;

    my $id;
    if ($attr =~ /^\S+\s+\"(\S+)\"/) {
      # GFF2
      $id = $1;
    } else {
      # assume GFF3
      foreach my $k_v (split(/;/, $attr)) {
        my ($k, $v) = split(/=/, $k_v);
        if ($k eq 'Name') {
          $id = $v;
          last;
        }
      }      
    }

    die "Could not identify feature name/id from line ($_)\n" if not $id; 
    die "Could not find slice for $chr\n" if not exists $self->slices->{$chr};
    push @features, $self->_create_simple_feature( $start, 
                                                   $end, 
                                                   ($strand eq '-') ? -1 : 1,
                                                   $id, 
                                                   $self->slices->{$chr}, 
                                                   $analysis );
    
  }
  return \@features;
}


sub parse_non_coding_genes_gff2 {
  my ($self, $file, $analysis, $type, $biotype ) = @_;

  my $fh = $self->_open_file($file);
  my $genes = $self->parse_non_coding_genes_gff2_fh($fh, $analysis, $type, $biotype);

  return $genes;
}


sub parse_non_coding_genes_gff2_fh {
  my ($self, $fh, $analysis, $types, $biotype ) = @_;

  my @genes;
  my ($transcripts, $tran2seq) = $self->_process_pseudo_files( $fh, $types );
  $self->verbose and print STDERR "There are " . keys(%$transcripts) . " raw transcripts for source $types\n";
  
  my ($processed_transcripts) = $self->_process_pseudo_transcripts( $transcripts, $analysis );
  $self->verbose and print STDERR "There are " . keys(%$processed_transcripts) . " processed transcripts for source $types\n";
  
  $biotype = $types if not defined $biotype; 
  
  my $genes = $self->_create_pseudo_transcripts($processed_transcripts, $biotype);
  $self->verbose and print STDERR "There are " . keys(%$genes) . " processed genes for source $types\n";
  
  foreach my $gene_id ( keys(%$genes) ) {
    my $transcripts = $genes->{$gene_id};
    my $gene = $self->_create_gene( $transcripts, $gene_id, $biotype );
    push( @genes, $gene );
  }
  return \@genes;
}


sub parse_protein_align_features_gff {
  my ($self, $gff_file, $analysis, $tag) = @_;

  my $fh = $self->_open_file($gff_file);
  my $hits = $self->parse_protein_align_features_gff($fh, $analysis, $tag);

  return $hits;
}


sub parse_protein_align_features_gff_fh {
  my ($self, $fh, $analysis, $tag) = @_;

  my @hits;
  while (<$fh>) {
    next if not /BLASTX\s+protein_match/;
    my @col = split(/\t/, $_);

    my $chr = $col[0];

    die "Could not find $chr in slice hash\n" if not exists $self->slices->{$chr};

    my $strand = $col[6] eq '+' ? 1 : -1;            # convert to ensembl style

    my ($type, $hitid, $hitst, $hiten);
    if ($col[8] =~ /\"Protein:(.+):(\S+)\"\s+(\d+)\d+(\d+)/) {
      # GFF2
      ($type, $hitid, $hitst, $hiten) = ($1, $2, $3, $4);
    } elsif ($col[8] =~ /Target=(.+):(\S+)\s+(\d+)\s+(\d+)/) {
      # assume GFF3
      ($type, $hitid, $hitst, $hiten) = ($1, $2, $3, $4);
    } else {
      $self->verbose and print STDERR "Could not extract sensible hit info from '$col[8]', so skippinh\n";
      next;
    }

    next if $type ne $tag;

    my $length = $col[4] - $col[3];
    ( $hitst, $hiten ) = sort { $a <=> $b } ( $hitst, $hiten );
    if ($hiten - $hitst != $length) {
      # fake it; unimportant
      $hiten = $hitst + $length;
    }

    push @hits, Bio::EnsEMBL::DnaPepAlignFeature->new(
      -start         => $col[3],
      -end           => $col[4],
      -strand        => $strand,
      -slice         => $self->slices->{$chr},
      -analysis      => $analysis,
      -score         => $col[5],
      -p_value       => 1,
      -percent_id    => 100,
      -hseqname      => $hitid,
      -hstart        => $hitst,
      -hend          => $hiten,
      -cigar_string  => "${length}M");        
  }

  return \@hits;
}


sub parse_dna_align_features_gff {
  my ($self, $gff_file, $analysis, $source, $feature, $orimap, $group) = @_;

  my $fh = $self->_open_file($gff_file);
  my $hits = $self->parse_dna_align_features_gff_fh($fh, $analysis, $source, $feature, $orimap, $group);

  return $hits;
}


sub parse_dna_align_features_gff_fh {
  my ($self, $fh, $analysis, $source, $feature, $orimap, $group) = @_;

  my (@hits, %hits_by_align_id, $align_count);
  while (<$fh>) {
    /^\#/ and next;

    my @col = split(/\t/, $_);

    next if $col[1] ne $source or $col[2] ne $feature;

    my $chr = $col[0];
    
    die "Could not find $chr in slice hash\n" if not exists $self->slices->{$chr};
    
    my $strand = $col[6] eq '+' ? 1 : -1;            # convert to ensembl style
    my ($hitid, $hitst, $hiten, $align_id);

    if ($col[8] =~ /\"Sequence:(\S+)\"\s+(\d+)\d+(\d+)/) {
      ($hitid, $hitst, $hiten) = ($1,$2,$3);
      $align_id = ++$align_count;
    } elsif ($col[8] =~ /Target=(\S+)\s+(\d+)\s+(\d+)/) {
      # GFF3
      ($hitid, $hitst, $hiten) = ($1, $2, $3);
      if ($col[8] =~ /ID=([^;]+)/) {
        $align_id = $1;
      } else {
        $align_id = ++$align_count;
      }
    } else {
      $self->verbose and print STDERR "Could not extract hit info from $col[8], so skipping\n";
      next;
    }
    
    my $length = $col[4] - $col[3];
    ( $hitst, $hiten ) = sort { $a <=> $b } ( $hitst, $hiten );
    if ($hiten - $hitst != $length) {
      # fake it; unimportant
      $hiten = $hitst + $length;
    }

    if (defined $orimap and 
        exists $orimap->{$hitid} and
        $orimap->{$hitid} eq '3') {
      $strand *= -1;
    }
    my $fp = new Bio::EnsEMBL::FeaturePair;
    $fp->score($col[5]);
    $fp->percent_id(100);
    $fp->slice($self->slices->{$chr});
    $fp->seqname($col[0]);
    $fp->start($col[3]);
    $fp->end($col[4]);
    $fp->strand($strand);
    $fp->hseqname($hitid);
    $fp->hstart($hitst);
    $fp->hend($hiten);
    $fp->hstrand(1);
    $fp->analysis($analysis);
    
    push @{$hits_by_align_id{$align_id}}, $fp;
  }

  foreach my $align_id (sort keys %hits_by_align_id) {
    my @feats = sort { $a->start <=> $b->start } @{$hits_by_align_id{$align_id}};

    my %strands;
    map { $strands{$_->strand} += $_->end - $_->start + 1 } @feats;
    if (scalar(keys %strands) > 1) {
      # sometimes 1-bp features can be dumped on the wrong strand
      my ($cons_strand) = sort { $strands{$b} <=> $strands{$a} } keys %strands;
      map { $_->strand($cons_strand) } @feats;
    }

    if ($group) {
      push @hits, Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => \@feats);
    } else {
      foreach my $f (@feats) {
        push @hits, Bio::EnsEMBL::DnaDnaAlignFeature->new(-features => [$f]);
      }
    }
    delete $hits_by_align_id{$align_id};
  }

  return \@hits;
}


sub parse_repeatfeatures_gff {
  my ($self, $file, $analysis, $source, $type) = @_;

  my $fh = $self->_open_file($file);
  my $feats = $self->parse_repeatfeatures_gff_fh($fh, $analysis, $source, $type);

  return $feats;
}


sub parse_repeatfeatures_gff_fh {
  my ($self, $fh, $analysis, $source, $type ) = @_;

  my $ff = Bio::EnsEMBL::Analysis::Tools::FeatureFactory->new();

  my @rfs;

  while(<$fh>){
    /^\#/ and next;
    
    my @col = split(/\t/, $_);
    
    next if $col[1] ne $source or $col[2] ne $type;
    
    my $chr = $col[0];
    
    die "Could not find $chr in slice hash\n" if not exists $self->slices->{$chr};

    my ($rep_name, $rep_start, $rep_end);

    if ($col[8] =~ /Target\s+\"(.+)\"\s+(\d+)\s+(\d+)/) {
      ($rep_name, $rep_start, $rep_end) = ($1, $2, $3);
    } elsif ($col[8] =~ /Target=(\S+)\s+(\d+)\s+(\d+)/) {
      ($rep_name, $rep_start, $rep_end) = ($1, $2, $3);
    } else {
      die "Could not get repeat details from GFF line ($col[8])\n";
    }

    $rep_name =~ s/Motif://; 
    my $rep_class = "Unknown";

    my $rc = $ff->create_repeat_consensus($rep_name, "Unknown");

    my $rf = $ff->create_repeat_feature( $col[3],
                                         $col[4],
                                         $col[6] eq '-' ? -1 : 1,
                                         $col[5],
                                         $rep_start, 
                                         $rep_end,
                                         $rc,
                                         $self->slices->{$chr}->seq_region_name,
                                         $self->slices->{$chr},
                                         $analysis,
        );

    push @rfs, $rf;
    
  }

  return \@rfs;
}

###########################################################

sub write_genes {
  my ($self, $genes, $exon_stable_ids ) = @_;

  $self->verbose and print STDERR "Writing ", scalar(@$genes), " genes...\n";

  my $gene_adaptor = $self->database_handle->get_GeneAdaptor;
  foreach my $g (@$genes) {
    eval {
      $gene_adaptor->store($g);
    };
    if ($@) {
      die("couldn't store " . $g->stable_id . " problems=:$@:\n");
    }
  }
  
  $self->verbose and print STDERR "Finished writing genes.\n";

  if ($exon_stable_ids) {
    my $sth = $self->database_handle->dbc->prepare("UPDATE exon set stable_id = ? WHERE exon_id = ?");

    foreach my $g (@$genes) {
      my @exons = sort { $a->start <=> $b->start } @{$g->get_all_Exons};
      if ($g->strand < 0) {
        @exons = reverse @exons;
      }
      for (my $edx=0; $edx < @exons; $edx++) {
        my $ex = $exons[$edx];
        my $db_id = $ex->dbID;
        my $e_stable_id = sprintf("%s.e%d", $g->stable_id, $edx+1);
        $sth->execute($e_stable_id, $db_id);
      }
    }
  }
}


sub write_simple_features {
  my ($self, $features ) = @_;

  eval { print "\n check 1: " . $$features[0]->display_label };
  eval { print "\n check 2: " . $$features[0]->start . " - " . $$features[0]->end . "\n" };
  
  my $feature_adaptor = $self->database_handle->get_SimpleFeatureAdaptor;
  
  eval { $feature_adaptor->store(@$features); };
  if ($@) {
    die "Could not store simple feature; problems: $@\n";
  }
}

sub write_dna_align_features {
  my ( $self, $feats ) = @_;

  $self->verbose and print STDERR "Storing " . scalar(@$feats) . " DnaDnaAlignFeatures\n";

  if (@$feats) {
    eval { print "\n check 1: " . $feats->[0]->display_label };
    eval { print "\n check 2: " . $feats->[0]->start . " - " . $feats->[0]->end . "\n" };
  
    eval {
      $self->database_handle->get_DnaAlignFeatureAdaptor->store(@$feats);
    };
    if ($@) {
      die "couldn't store features , problems: " . $@;
    }
  }
  return 1;
}


sub write_protein_align_features {
  my ($self, $feats) = @_;
  
  $self->verbose and print STDERR "Storing " . scalar(@$feats) . " DnaPepAlignFeatures\n";

  if (@$feats)  {
    eval { print "\n check 1: " . $feats->[0]->display_label };
    eval { print "\n check 2: " . $feats->[0]->start . " - " . $feats->[0]->end . "\n" };
  
    eval {
      $self->database_handle->get_ProteinAlignFeatureAdaptor->store(@$feats);
    };
    if ($@) {
      die "couldn't store features , problems: " . $@;
    }
  }
  return 1;
}

sub write_repeat_features {
  my ($self, $feats) = @_;
  
  $self->verbose and print STDERR "Storing " . scalar(@$feats) . " RepeatFeatures\n";

  if (@$feats)  {
    eval { print "\n check 1: " . $feats->[0]->display_label };
    eval { print "\n check 2: " . $feats->[0]->start . " - " . $feats->[0]->end . "\n" };
  
    eval {
      $self->database_handle->get_RepeatFeatureAdaptor->store(@$feats);
    };
    if ($@) {
      die "couldn't store features , problems: " . $@;
    }
  }
  return 1;
}



############################################
#
# Get/Set
#
############################################

#####################################
sub species {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_species} = $val;
  }
  if (not exists $self->{_species} or not $self->{_species}) {
    return 'elegans';
  } else {
    return $self->{_species};
  }
}

#####################################
sub debug {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_debug} = $val;
  }
  if (not exists $self->{_debug}) {
    return 0;
  } else {
    return $self->{_debug};
  }
}

#####################################
sub recognise_sources {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_recognise_sources} = $val;
  }
  if (not exists $self->{_recognise_sources}) {
    return 0;
  } else {
    return $self->{_recognise_sources};
  }
}



#####################################
sub verbose {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_verbose} = $val;
  }
  if (not exists $self->{_verbose}) {
    return 0;
  } else {
    return $self->{_verbose};
  }
}

#####################################
sub database_handle {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_dbh} = $val;
  }
  if (not exists $self->{_dbh}) {
    return 0;
  } else {
    return $self->{_dbh};
  }
}

#####################################
sub ignore_gff_phases {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_ignore_gff_phases} = $val;
  }
  if (not exists $self->{_ignore_gff_phases}) {
    return 0;
  } else {
    return $self->{_ignore_gff_phases};
  }
}



#####################################
sub slices {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{_slices} = $val;
  }
  if (not exists $self->{_slices}) {
    my $slice_hash = {};
    foreach my $slice (@{$self->database_handle->get_SliceAdaptor->fetch_all('toplevel')}) {
      $slice_hash->{$slice->seq_region_name} = $slice;
      if ($self->species eq 'elegans') {
        my $other_name;
        if ($slice->seq_region_name !~ /^CHROMOSOME/) {
          $other_name = "CHROMOSOME_" . $slice->seq_region_name; 
        } else {
          $other_name = $slice->seq_region_name;
          $other_name =~ s/^CHROMOSOME_//;
        }
        $slice_hash->{$other_name} = $slice;
      } elsif ($self->species eq 'briggsae') {
        my $other_name;
        if ($slice->seq_region_name !~ /^chr/) {
          $other_name = "chr" . $slice->seq_region_name; 
        } else {
          $other_name = $slice->seq_region_name;
          $other_name =~ s/^chr//;
        }
        $slice_hash->{$other_name} = $slice;
      }
    }

    $self->{_slices} = $slice_hash;
  }
  return $self->{_slices};

}


############################################
#
# Internal methods
#
############################################

sub _process_protein_coding_file {
  my ($self, $fh) = @_;
  
  $self->verbose and print STDERR "Reading GFF fh...\n";

  my (%coding_isoforms, %non_coding_isoforms, %five_prime, %three_prime, %parent_seqs);

  while (<$fh>) {

    #CHROMOSOME_I    curated three_prime_UTR 11696828        11697110        .       -       .       Transcript "T22H2.5a"
    #CHROMOSOME_I    curated three_prime_UTR 11697157        11697230        .       -       .       Transcript "T22H2.5a"
    #CHROMOSOME_I    curated five_prime_UTR  11698944        11698954        .       -       .       Transcript "T22H2.5a"

    /^\#/ and next;

    chomp;
    my $element = $_;
    my ( $chr, $status, $type, $start, $end, $score, $strand, $frame, $attr ) = split(/\t/, $element);
    my ($isoform) = $attr =~ /^\S+\s+\"(\S+)\"/;

    my $line = $status . " " . $type;

    next if $status ne 'Coding_transcript' and $status ne 'curated' and $status ne 'Non_coding_transcript';
    
    if ( $line eq 'Coding_transcript protein_coding_primary_transcript') {      
      my $transcript = $isoform;
      
      if ($self->species and $self->species !~ /elegans/) {
        $isoform =~ s/^(\w+)\.\d+$/$1/;
      } else {
        $isoform =~ s/(\.\w+)\.\d+$/$1/;
      }
      
      if (not $five_prime{$isoform}) {
        $five_prime{$isoform} = {};
      }
      if (not $five_prime{$isoform}{$transcript}) {
        $five_prime{$isoform}{$transcript} = [];
      }
      
      if (not $three_prime{$isoform}) {
        $three_prime{$isoform} = {};
      }
      if (not $three_prime{$isoform}{$transcript}) {
        $three_prime{$isoform}{$transcript} = [];
      }
    
    } elsif ( ( $line eq 'Coding_transcript five_prime_UTR' ) or ( $line eq 'Coding_transcript three_prime_UTR' ) ) {
      my $transcript = $isoform;
      
      #remove transcript-specific part: Y105E8B.1a.2
      if ($self->species and $self->species !~ /elegans/) {
        $isoform =~ s/^(\w+)\.\d+$/$1/;
      } else {
        $isoform =~ s/(\.\w+)\.\d+$/$1/;
      }
      
      my ($position) = $type;
      if ( $position =~ /^five/ ) {
        
        if ( !$five_prime{$isoform} ) {
          $five_prime{$isoform} = {};
          if ( !$five_prime{$isoform}{$transcript} ) {
            $five_prime{$isoform}{$transcript} = [];
          }
          push( @{ $five_prime{$isoform}{$transcript} }, $element );
        }
        else {
          if ( !$five_prime{$isoform}{$transcript} ) {
            $five_prime{$isoform}{$transcript} = [];
          }
          push( @{ $five_prime{$isoform}{$transcript} }, $element );
        }
      }
      elsif ( $position =~ /^three/ ) {
        
        if ( !$three_prime{$isoform} ) {
          $three_prime{$isoform} = {};
          if ( !$three_prime{$isoform}{$transcript} ) {
            $three_prime{$isoform}{$transcript} = [];
          }
          push( @{ $three_prime{$isoform}{$transcript} }, $element );
        }
        else {
          if ( !$three_prime{$isoform}{$transcript} ) {
            $three_prime{$isoform}{$transcript} = [];
          }
          push( @{ $three_prime{$isoform}{$transcript} }, $element );
        }
      }
    } elsif ( $line eq 'curated coding_exon') {
      if ( !$coding_isoforms{$isoform} ) {
        $coding_isoforms{$isoform} = [];
        $parent_seqs{$isoform} = $chr;
      }
      push( @{ $coding_isoforms{$isoform} }, $element );
    } elsif ($line eq 'Non_coding_transcript exon') {
      if ( !$non_coding_isoforms{$isoform} ) {
        $non_coding_isoforms{$isoform} = [];
        $parent_seqs{$isoform} = $chr;
      }
      push( @{ $non_coding_isoforms{$isoform} }, $element );
    }
  }

  #
  # Add "empty" UTRs for all of the genes that did not have one; downstream code needs this
  #
  foreach my $isoform (keys %coding_isoforms) {
    if (not $five_prime{$isoform} ) {
      $five_prime{$isoform} = {};
      $five_prime{$isoform}{$isoform} = [];
    }
    if (not $three_prime{$isoform}) {
      $three_prime{$isoform} = {};
      $three_prime{$isoform}{$isoform} = [];
    }
  }

  $self->verbose and printf STDERR "Finished reading the GFF file for protein-coding genes; got %d CDS isoforms and %d non-coding isoforms\n", scalar(keys %coding_isoforms), scalar(keys %non_coding_isoforms);

  return (\%coding_isoforms, \%non_coding_isoforms, \%five_prime, \%three_prime, \%parent_seqs);
}


sub _generate_protein_coding_transcripts {
    my ($self,  $genesRef, $analysis, $five_prime, $three_prime, $parent_seqs ) = @_;
    my %genes;
    my %transcripts;
    my %temp_transcripts;
    my %five_trans_start;
    my %three_trans_end;
    my %trans_start_exon;
    my %trans_end_exon;
    my $translation_end;
    my $genecount = 0;
    my @global_exons;
    my %overlapcheck;

    #go through all genes
    GENE: foreach my $gene_name ( keys(%$genesRef) ) {
      if (not exists $parent_seqs->{$gene_name}) {
        die "Could not find parent sequence for $gene_name\n";
      }
      my $parent_seq = $parent_seqs->{$gene_name};

      if (not exists $self->slices->{$parent_seq}) {
        die "Could not find slice object for $parent_seq\n";
      }
      my $slice = $self->slices->{$parent_seq};

      #create gene-hash entry
      $genes{$gene_name} = [];
      my $transcriptcount = 0;
      %temp_transcripts = ();
      
      ## collect all "curated_coding_exons" for this gene ##
      my @lines        = @{ $$genesRef{$gene_name} };    #is this right?
      my @global_exons = ();
      my %three_prime_exons;
      my %five_prime_exons;
      
      foreach my $line (@lines) {
        my ( $chr, $status, $type, $start, $end, $score, $strand, $frame, $attr ) = split(/\t/, $line);

        my $exon  = new Bio::EnsEMBL::Exon;

        $exon->start($start);
        $exon->end($end);
        $exon->analysis($analysis);
        $exon->slice($slice);
        if ($frame eq '.') {
          $exon->phase(-1);
        } else {
          my $phase = ( 3 - $frame ) % 3;
          $exon->phase($phase);
          my $end_phase = ( $phase + ( $exon->end - $exon->start ) + 1 ) % 3;       
          $exon->end_phase($end_phase);
        }
        if ( $strand eq '+' ) {
          $exon->strand(1);
        }
        else {
          $exon->strand(-1);
        }
        
        #$exon->score(100);
        push( @global_exons, $exon );
      }
      
      #sort exons for this gene
      if ( $global_exons[0]->strand == -1 ) {
        @global_exons = sort { $b->start <=> $a->start } @global_exons;
      }
      else {
        @global_exons = sort { $a->start <=> $b->start } @global_exons;
      }
      
      #collect 5' UTRs
      foreach my $transcript_name ( keys( %{ $$five_prime{$gene_name} } ) ) {
        
        my @five_prime_exons = ();
        %overlapcheck = ();
        
        #more than one transcript at 5 prime level
        $temp_transcripts{$transcript_name} = 1;
        
        #get UTR lines
        foreach my $line ( @{ $$five_prime{$gene_name}{$transcript_name} } ) {
          my ( $chr, $status, $type, $start, $end, $score, $strand, $frame, $attr ) = split( /\t/, $line );
          my ($gene) = $attr =~ /^\S+\s+\"(\S+)\"/;

          #avoid saving two identical exons
          if ( defined $overlapcheck{$start} ) {
            next;
          }
          $overlapcheck{$start} = 1;
          my $exon  = new Bio::EnsEMBL::Exon;
          my $phase = -1;
          $exon->start($start);
          $exon->end($end);
          $exon->analysis($analysis);
          $exon->slice($slice);
          $exon->phase($phase);
          my $end_phase = -1;
          $exon->end_phase($end_phase);
          
          if ( $strand eq '+' ) {
            $exon->strand(1);
          }
          else {
            $exon->strand(-1);
          }
          push( @five_prime_exons, $exon );
        }
        
        #sort exons for this transcript
        if ( @five_prime_exons and $five_prime_exons[0]->strand == -1 ) {
          @five_prime_exons = sort { $b->start <=> $a->start } @five_prime_exons;
        }
        else {
          @five_prime_exons = sort { $a->start <=> $b->start } @five_prime_exons;
        }
        
        #save them to transcript
        $five_prime_exons{$transcript_name} = \@five_prime_exons;
      }
      
      ## collect 3' UTRs ##
      foreach my $transcript_name ( keys( %{ $$three_prime{$gene_name} } ) ) {
        
        my @three_prime_exons = ();
        %overlapcheck = ();
        
        #more than one transcript at the 3 prime level, save the name
        $temp_transcripts{$transcript_name} = 1;
        
        #get UTR lines
        foreach my $line ( @{ $$three_prime{$gene_name}{$transcript_name} } ) {
          my ( $chr, $status, $type, $start, $end, $score, $strand, $frame, $attr ) = split(/\t/, $line);
          my ($gene) = $attr =~ /^\S+\s+\"(\S+)\"/;

          #avoid saving two identical exons
          if ( defined $overlapcheck{$start} ) {
            next;
          }
          $overlapcheck{$start} = 1;
          my $exon  = new Bio::EnsEMBL::Exon;
          my $phase = -1;
          $exon->start($start);
          $exon->end($end);
          $exon->analysis($analysis);
          $exon->slice($slice);
          $exon->phase($phase);
          my $end_phase = -1;
          $exon->end_phase($end_phase);
          
          if ( $strand eq '+' ) {
            $exon->strand(1);
          }
          else {
            $exon->strand(-1);
          }
          push( @three_prime_exons, $exon );
        }
        
        #sort exons for this transcript
        if ( @three_prime_exons and $three_prime_exons[0]->strand == -1 ) {
          @three_prime_exons = sort { $b->start <=> $a->start } @three_prime_exons;
        }
        else {
          @three_prime_exons = sort { $a->start <=> $b->start } @three_prime_exons;
        }

        #save them to transcript
        $three_prime_exons{$transcript_name} = \@three_prime_exons;
      }
      
      ## combine exons, 5' and 3' for every transcript ##
      foreach my $transcript_name ( keys %temp_transcripts ) {
        $transcriptcount++;
        my @exons = ();
        foreach my $temp_exon (@global_exons) {
          push @exons, $self->_clone_Exon($temp_exon);
        }
        my $translation_start = 1;
        my $first             = 1;
        
        #set default translation range
        $trans_start_exon{$transcript_name} = 0;
        $trans_end_exon{$transcript_name}   = $#exons;
        #print "\ntrans-exons: " . $trans_start_exon{$transcript_name} . " - " . $trans_end_exon{$transcript_name} . " (" . scalar @exons . ")";
        
        #check 5' exons
        if ( @{$five_prime_exons{$transcript_name}} ) {
          my @five_prime_exons = @{ $five_prime_exons{$transcript_name} };
          
          #is last 5' UTR exon part of first coding exon?
          
          FIVEUTR: while ( my $five_prime_exon = shift(@five_prime_exons) ) {

            my $start  = $five_prime_exon->start;
            my $end    = $five_prime_exon->end;
            my $strand = $five_prime_exon->strand;
            
            if ( $exons[ $trans_start_exon{$transcript_name} ]->strand == 1 and $strand == 1 ) {
              
              #forward strand
              if ( $start > $end ) {
                next FIVEUTR;
              }
              if ( $end == ( $exons[ $trans_start_exon{$transcript_name} ]->start ) - 1 ) {
                
                #combine exons, adjust translation start
                $translation_start = $exons[ $trans_start_exon{$transcript_name} ]->start - $start + 1;
                
                $five_trans_start{$transcript_name} = $translation_start;
                $exons[ $trans_start_exon{$transcript_name} ]->start($start);
                $exons[ $trans_start_exon{$transcript_name} ]->phase(-1);
              }
              elsif ( $end < $exons[ $trans_start_exon{$transcript_name} ]->start - 1 ) {
                
                #additional non-coding exon
                #add to exon array, keep translation start on last coding exon
                $trans_start_exon{$transcript_name}++;
                $trans_end_exon{$transcript_name}++;
                unshift( @exons, $five_prime_exon );
                $self->debug and print "\nadditional non-coding exon (+) " . $start . " - " . $end;
              }
              else {
                $self->debug and print STDERR "\n>>$transcript_name strange 5' UTR exon (+): $start - $end with 1.exons at "
                    . $exons[ $trans_start_exon{$transcript_name} ]->start . " - "
                    . $exons[ $trans_start_exon{$transcript_name} ]->end;
                next FIVEUTR;
              }
            }
            elsif ( $exons[ $trans_start_exon{$transcript_name} ]->strand == -1 and $strand == -1 ) {
              
              #reverse strand
              if ( $start > $end ) {
                next FIVEUTR;
              }
              if ( $start == ( $exons[ $trans_start_exon{$transcript_name} ]->end ) + 1 ) {
                
                #combine exons, adjust translation start
                $translation_start = ( $end - $exons[ $trans_start_exon{$transcript_name} ]->end + 1 );
                
                $five_trans_start{$transcript_name} = $translation_start;
                $exons[ $trans_start_exon{$transcript_name} ]->end($end);
                $exons[ $trans_start_exon{$transcript_name} ]->phase(-1);
              }
              elsif ( $start > ( $exons[ $trans_start_exon{$transcript_name} ]->end ) + 1 ) {
                
                #additional non-coding exon
                #add to exon array, keep translation start on last coding exon
                $trans_start_exon{$transcript_name}++;
                $trans_end_exon{$transcript_name}++;
                unshift( @exons, $five_prime_exon );
                
              }
              else {
                $self->debug and print STDERR "\n>>$transcript_name strange 5' UTR exon (-): $start - $end with 1.exons at "
                    . $exons[ $trans_start_exon{$transcript_name} ]->start . " - "
                    . $exons[ $trans_start_exon{$transcript_name} ]->end;
                next FIVEUTR;
              }
            }
            else {
              $self->debug and print STDERR "\n>>strand switch in UTR / coding!";
            }
          }
        }
        
        #check 3' exons
        if ( @{$three_prime_exons{$transcript_name}} ) {
          my @three_prime_exons = @{ $three_prime_exons{$transcript_name} };
          
          #is first 3' UTR exon part of last coding exon?
          
          THREEUTR: while ( my $three_prime_exon = shift(@three_prime_exons) ) {
          
            my $start  = $three_prime_exon->start;
            my $end    = $three_prime_exon->end;
            my $strand = $three_prime_exon->strand;
            
            if ( $exons[ $trans_end_exon{$transcript_name} ]->strand == 1 and $strand == 1 ) {
              
              #forward strand
              if ( $start > $end ) {
                $self->debug and print STDERR "\n>>$transcript_name strange 3' UTR (+) exon: " . $start . " - " . $end;
                next THREEUTR;
              }
              if ( $start == ( ( $exons[ $trans_end_exon{$transcript_name} ]->end ) + 1 ) ) {
                
                #combine exons, adjust translation start
                $translation_end =
                    ( ( $exons[ $trans_end_exon{$transcript_name} ]->end - $exons[ $trans_end_exon{$transcript_name} ]->start ) + 1 );
                
                $three_trans_end{$transcript_name} = $translation_end;
                $exons[ $trans_end_exon{$transcript_name} ]->end($end);
                $exons[ $trans_end_exon{$transcript_name} ]->end_phase(-1);
                
              }
              elsif ( $start > ( ( $exons[ $trans_end_exon{$transcript_name} ]->end ) + 1 ) ) {
                
                #additional non-coding exon
                #add to exon array
                push( @exons, $three_prime_exon );
                
              }
              else {
                $self->debug and print STDERR "\n$transcript_name strange 3' UTR exon (+): $start - $end with 1.exons at "
                    . $exons[ $trans_end_exon{$transcript_name} ]->start;
                next THREEUTR;
              }
            }
            elsif ( $exons[ $trans_end_exon{$transcript_name} ]->strand == -1 and $strand == -1 ) {
              
              #reverse strand
              if ( $start > $end ) {
                $self->debug and print STDERR "\n>>$transcript_name strange 3' UTR (-) exon: " . $start . " - " . $end;
                next THREEUTR;
              }
              if ( $end == ( ( $exons[ $trans_end_exon{$transcript_name} ]->start ) - 1 ) ) {
                
                #combine exons, keep translation start
                $translation_end =
                    ( ( $exons[ $trans_end_exon{$transcript_name} ]->end - $exons[ $trans_end_exon{$transcript_name} ]->start ) + 1 );
                
                $three_trans_end{$transcript_name} = $translation_end;
                $exons[ $trans_end_exon{$transcript_name} ]->start($start);
                $exons[ $trans_end_exon{$transcript_name} ]->end_phase(-1);
              }
              elsif ( $end < ( ( $exons[ $trans_end_exon{$transcript_name} ]->start ) - 1 ) ) {
                
                #additional non-coding exon
                #add to exon array
                push( @exons, $three_prime_exon );
                
              }
              else {
                $self->debug and print STDERR "\n$transcript_name strange 3' UTR exon (-): $start - $end with 1.exons at "
                    . $exons[ $trans_end_exon{$transcript_name} ]->start;
                next THREEUTR;
              }
            }
          }
          
        }
        
        #add exon data to transcript
        $transcripts{$transcript_name} = \@exons;
        
      }
    }

    return ( \%transcripts, \%five_trans_start, \%three_trans_end, \%trans_start_exon, \%trans_end_exon );
}



sub _create_protein_coding_transcripts {
  my ($self,  $transcriptsRef, $five_start, $three_end, $trans_start_exon, $trans_end_exon ) = @_;
  
  my %transcripts = %$transcriptsRef;
  my @non_translate;
  my %genes;
  my $gene_name;
  my $transcript_id;
  foreach my $transcript ( keys(%transcripts) ) {
    my $time  = time;
    my @exons = @{ $transcripts{$transcript} };
    
    #get the gene-name

    if ($self->species and $self->species !~ /elegans/) {
      $gene_name= ( $transcript =~ /^([A-Za-z]+\d+)/ ) ? $1 : $transcript;
    } else {
      $gene_name= ( $transcript =~ /^([A-Za-z0-9_]+\.\d+)/ ) ? $1 : $transcript; 
    }
    
    $transcript_id = $transcript;

    my $transcript  = new Bio::EnsEMBL::Transcript;
    my $translation = new Bio::EnsEMBL::Translation;
    my @sorted_exons;
    if ( $exons[0]->strand == 1 ) {
      @sorted_exons = sort { $a->start <=> $b->start } @exons;
    }
    else {
      @sorted_exons = sort { $b->start <=> $a->start } @exons;
    }
    my $exon_count = 1;
    my $phase      = 0;
    foreach my $exon (@sorted_exons) {
      $exon->stable_id( $transcript_id . "." . $exon_count );
      $exon_count++;
      eval {
        $transcript->add_Exon($exon);
      };
      if ($@) { 
        print STDERR "\n>>$transcript_id EXON ERROR: " . $@ . "\n";
      }
    }
    my $start_exon_ind;
    if ( exists( $trans_start_exon->{$transcript_id} ) ) {
      $self->debug and 
          print STDERR "\n adjusting coding exons to " . $trans_start_exon->{$transcript_id} . " - " . $trans_end_exon->{$transcript_id};
      $translation->start_Exon( $sorted_exons[ $trans_start_exon->{$transcript_id} ] );
      $translation->end_Exon( $sorted_exons[ $trans_end_exon->{$transcript_id} ] );
      $start_exon_ind = $trans_start_exon->{$transcript_id};
    }
    else {
      $self->debug and print STDERR  "\n no UTRs - setting translation to span whole transcript\n";
      $translation->start_Exon( $sorted_exons[0] );
      $translation->end_Exon( $sorted_exons[$#sorted_exons] );
      $start_exon_ind = 0;
    }
    
    if ( exists( $five_start->{$transcript_id} ) ) {
      $self->debug and print STDERR "1 setting translation start on transcript " . $transcript_id . " to " . $five_start->{$transcript_id} . "\n";
      $translation->start( $five_start->{$transcript_id} );
    }
    else {
      $self->debug and print STDERR "1 setting translation start on transcript $transcript_id to 1\n";
      $translation->start(1);
    }
    
    if ( ( !defined( $translation->start ) ) or ( $translation->start <= 0 ) ) {
      $self->debug and do {
        print STDERR ">> no translation start info for " . $transcript_id;
        print STDERR ".." . $five_start->{$transcript_id} . "\n";
      };
      die();
    }
    
    if ( exists( $three_end->{$transcript_id} ) ) {
      $self->debug and print STDERR "2 setting translation end on transcript " . $transcript_id . " to " . $three_end->{$transcript_id} . " (1)\n";
      $translation->end( $three_end->{$transcript_id} );
    }
    else {
      if ( defined( $trans_end_exon->{$transcript_id} ) ) {
        $translation->end( $sorted_exons[ $trans_end_exon->{$transcript_id} ]->end - $sorted_exons[ $trans_end_exon->{$transcript_id} ]->start + 1 );
        $self->debug and print STDERR "2 setting translation end on transcript "
            . $transcript_id . " to "
            . ( $exons[ $trans_end_exon->{$transcript_id} ]->end - $exons[ $trans_end_exon->{$transcript_id} ]->start + 1 );
      }
      else {
        $translation->end( $sorted_exons[$#sorted_exons]->end - $sorted_exons[$#sorted_exons]->start + 1 );
        $self->debug and print STDERR "2 setting translation end on transcript "
            . $transcript_id . " to "
            . ( $sorted_exons[$#sorted_exons]->end - $sorted_exons[$#sorted_exons]->start + 1 );
      }
    }
    
    $translation->stable_id($transcript_id);
    $translation->version(undef);
    $transcript->translation($translation);
    $transcript->stable_id($transcript_id);
    $transcript->biotype('protein_coding');
    if ( !$genes{$gene_name} ) {
      $genes{$gene_name} = [];
      push( @{ $genes{$gene_name} }, $transcript );
    }
    else {
      push( @{ $genes{$gene_name} }, $transcript );
    }
  }
  return \%genes;
}



sub _create_gene {
  my ($self,  $transcripts, $name, $biotype ) = @_;
  my $time     = time;
  my $gene     = new Bio::EnsEMBL::Gene;
  my $analysis = $transcripts->[0]->get_all_Exons->[0]->analysis;
  
  $gene->analysis($analysis);
  $gene->biotype( $biotype ) if defined $biotype;
  $gene->stable_id($name);
  $gene->source("WormBase");
  foreach my $transcript (@$transcripts) {
    $transcript->analysis($analysis) if not $transcript->analysis;
    $gene->add_Transcript($transcript);
  }
  
  return $gene;
}


sub _create_simple_feature {
  my ($self, $start, $end, $strand, $id, $seq, $analysis ) = @_;
  
  #warn "first: $start, $end, $strand, $id...";
  
  my $simple_feature = Bio::EnsEMBL::SimpleFeature->new();
  $simple_feature->start($start);
  $simple_feature->strand($strand);
  $simple_feature->end($end);
  $simple_feature->display_label($id);
  $simple_feature->slice($seq);
  $simple_feature->analysis($analysis);
  
  return $simple_feature;
}


sub _process_pseudo_files {
  my ($self, $fh, $source_to_look_for ) = @_;
  my (%transcripts);

  while (<$fh>) {
    /^\#/ and next;

    chomp;

    my $line = $_;
    my ( $chr, $source, $type, $start, $end, $score, $strand, $frame, $attr ) 
        = split(/\t/, $_);
    
    next if $source ne $source_to_look_for or $type ne 'exon';

    my ($gene) = $attr =~ /^\S+\s+\"(\S+)\"/;
    
    if ( !$transcripts{$gene} ) {
      $transcripts{$gene} = [];
      push( @{ $transcripts{$gene} }, $line );
    }
    else {
      push( @{ $transcripts{$gene} }, $line );
    }
  }
  return \%transcripts;
}

sub _process_pseudo_transcripts {
  my ($self, $transcripts, $analysis ) = @_;
  
  my %genes;
  my %transcripts = %$transcripts;
  my @names       = keys(%transcripts);
  
  foreach my $name (@names) {
    my @lines = @{ $transcripts{$name} };
    $transcripts{$name} = [];
    my @exons;
    foreach my $line (@lines) {
      
      my ( $chr, $status, $type, $start, $end, $score, $strand, $frame, $sequence, $gene ) = split /\s+/, $line;
      
      if ( $start == $end ) {
        next;
      }
      
      die "Could not not find slice for '$chr' in hash\n" if not exists $self->slices->{$chr};
      
      my $exon = new Bio::EnsEMBL::Exon;
      
      # non-coding/pseudogene exons always have a phase -1/-1
      
      my $phase = -1; 
      my $end_phase = -1;
      
      $exon->start($start);
      $exon->end($end);
      $exon->analysis($analysis);
      $exon->slice($self->slices->{$chr});
      $exon->phase($phase);
      $exon->end_phase($end_phase);
      if ( $strand eq '+' ) {
        $exon->strand(1);
      }
      else {
        $exon->strand(-1);
      }
      
      #$exon->score(100);
      push( @exons, $exon );
    }
    if ( $exons[0]->strand == -1 ) {
      @exons = sort { $b->start <=> $a->start } @exons;
    }
    else {
      @exons = sort { $a->start <=> $b->start } @exons;
    }
    my $phase = 0;
    foreach my $e (@exons) {
      push( @{ $transcripts{$name} }, $e );

    }
  }
  
  return ( \%transcripts );
}

sub _create_pseudo_transcripts {
  my ($self, $transcripts, $biotype) = @_;
  
  my %transcripts = %$transcripts;
  my %genes;
  my $gene_name;
  my $transcript_id;

  foreach my $transcript ( keys(%transcripts) ) {
    my $time  = time;
    my @exons = @{ $transcripts{$transcript} };


    if ($self->species and $self->species !~ /elegans/) {
      $gene_name= ( $transcript =~ /^([A-Za-z]+\d+)/ ) ? $1 : $transcript;
    } else {
      $gene_name= ( $transcript =~ /^([A-Za-z0-9_]+\.\d+)/ ) ? $1 : $transcript; 
    }

    $transcript_id = $transcript;

    my $transcript = new Bio::EnsEMBL::Transcript;
    my @sorted_exons;
    if ( $exons[0]->strand == 1 ) {
      @sorted_exons = sort { $a->start <=> $b->start } @exons;
    }
    else {
      @sorted_exons = sort { $b->start <=> $a->start } @exons;
    }
    my $exon_count = 1;
    my $phase      = 0;
    foreach my $exon (@sorted_exons) {
      $exon->stable_id( $transcript_id . "." . $exon_count );
      $exon_count++;
      $transcript->add_Exon($exon);
    }
    $transcript->stable_id($transcript_id);
    $transcript->biotype($biotype) if defined $biotype;
    $transcript->source("WormBase");

    if ( !$genes{$gene_name} ) {
      $genes{$gene_name} = [];
      push( @{ $genes{$gene_name} }, $transcript );
    }
    else {
      push( @{ $genes{$gene_name} }, $transcript );
    }
  }
  return \%genes;
  
}


sub _selenocysteine {
  my ($self, $transcript, $pos ) = @_;
  
  my $seq_edit = Bio::EnsEMBL::SeqEdit->new(
    -CODE    => '_selenocysteine',
    -NAME    => 'Selenocysteine',
    -DESC    => 'Selenocysteine',
    -START   => $pos,
    -END     => $pos,
    -ALT_SEQ => 'U'                  #the one-letter symbol for selenocysteine
      );
  my $attribute         = $seq_edit->get_Attribute();
  my $translation       = $transcript->translation();
  my $attribute_adaptor = $self->database_handle->get_AttributeAdaptor();
  $attribute_adaptor->store_on_Translation( $translation, [$attribute] );
}


sub _display_exons {
  my ($self, @exons) = @_;
  
  @exons = sort { $a->start <=> $b->start || $a->end <=> $b->end } @exons if ( $exons[0]->strand == 1 );
  
  @exons = sort { $b->start <=> $a->start || $b->end <=> $a->end } @exons if ( $exons[0]->strand == -1 );
  
  foreach my $e (@exons) {
    print STDERR $e->stable_id . "\t " . $e->start . "\t " . $e->end . "\t " . $e->strand . "\t " . $e->phase . "\t " . $e->end_phase . "\n";
  } 
}


sub _display_non_translate {
  my ($self, @transcripts) = @_;
  
  foreach my $t (@transcripts) {
    
    my @exons = @{ $t->get_all_translateable_Exons };
    
    foreach my $e (@exons) {
      print "exon " . $e->stable_id . " " . $e->start . " " . $e->end . " " . $e->strand . "\n";
      if ($e->end - $e->start + 1 < 5) {
        print " Exon too short; not attempting to translate\n";
        next;
      }
      my $seq  = $e->seq;
      my $pep0 = $seq->translate( '*', 'X', 0 );
      my $pep1 = $seq->translate( '*', 'X', 1 );
      my $pep2 = $seq->translate( '*', 'X', 2 );
      print "exon sequence :\n" . $e->seq->seq . "\n\n";
      print $e->seqname . " " . $e->start . " : " . $e->end . " translation in 0 frame\n " . $pep0->seq . "\n\n";
      print $e->seqname . " " . $e->start . " : " . $e->end . " translation in 1 phase\n " . $pep2->seq . "\n\n";
      print $e->seqname . " " . $e->start . " : " . $e->end . " translation in 2 phase\n " . $pep1->seq . "\n\n";
      print "\n\n";
      
    }
  }
}

sub _clone_Exon {
  my ($self, $exon) = @_;

  my $newexon = Bio::EnsEMBL::Exon->new();
  $newexon->start      ($exon->start);
  $newexon->end        ($exon->end);
  $newexon->phase      ($exon->phase);
  $newexon->end_phase  ($exon->end_phase);
  $newexon->strand     ($exon->strand);
  $newexon->analysis   ($exon->analysis);
  $newexon->dbID       ($exon->dbID);
  $newexon->slice      ($exon->slice);
  $newexon->stable_id  ($exon->stable_id);
  $newexon->version    ($exon->version);

  return $newexon;
}


#############################################
#
#############################################
sub _open_file {
  my ($self, $file) = @_;

  my $fh;
  if ($file =~ /\.gz$/) {
    open($fh, "gunzip -c $file |") or die "Could not open gunzip stream to $file\n";
  } else {
    open($fh, $file) or die "Could not open $file for reading\n";
  }

  return $fh;
}




1;
