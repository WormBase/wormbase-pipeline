#!/usr/bin/env perl

use strict;

use Getopt::Long;
use Bio::EnsEMBL::Registry;

my ($species, 
    $reg_conf,
    %p2g_map, $p2g_file, 
    $pname_file, %wb_oset_names,
    $choose_best,
    $acefile, $ace_fh,
    $verbose);

&GetOptions(
  'reg_conf=s'      => \$reg_conf,
  'species=s'       => \$species,
  'probe2genemap=s' => \$p2g_file,
  'wbprobenames=s'  => \$pname_file,
  'choosebest'      => \$choose_best,
  'verbose'         => \$verbose,
  'acefile=s'       => \$acefile,
    );


die "You must supply a valid Registry file\n" if not defined $reg_conf or not -e $reg_conf;
die "You must supply a species\n" if not defined $species;

Bio::EnsEMBL::Registry->load_all($reg_conf);

if (defined $p2g_file) {
  &parse_probe_to_gene_map($p2g_file, \%p2g_map);
}
if (defined $pname_file) {
  &parse_probename_file($pname_file, \%wb_oset_names);
}


my $core_db = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
my $fgen_db = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');


my $a_adap = $fgen_db->get_ArrayAdaptor;
my $p_adap = $fgen_db->get_ProbeAdaptor;
my $sl_adap = $core_db->get_SliceAdaptor;


my %arrays_by_vendor;

map { push @{$arrays_by_vendor{$_->vendor}}, $_ } @{$a_adap->fetch_all};

if ($acefile) {
  open($ace_fh, ">$acefile") or die "Could not open $acefile for writing\n";
} else {
  $ace_fh = \*STDOUT;
}


foreach my $vendor (sort { $b cmp $a } keys %arrays_by_vendor) {
  my (%results);

  print STDERR "Processing $vendor...\n" if $verbose;

  if ($vendor ne 'AFFY') {    
    my (%probes_by_dbID, %done_wb_probes);
    
    foreach my $array (@{$arrays_by_vendor{$vendor}}) {
      foreach my $probe (@{$array->get_all_Probes}) {
        push @{$probes_by_dbID{$probe->dbID}}, $probe;
      }
    }

    foreach my $probe_group_ref (values %probes_by_dbID) {
      my $probe = $probe_group_ref->[0];

      my %wb_names;

      foreach my $probe (@$probe_group_ref) {
        my ($array) = @{$probe->get_all_Arrays};
        foreach my $pname (@{$probe->get_all_probenames($array->name)}) {
          if (exists $wb_oset_names{$pname}) {
            $wb_names{$pname} = 1;
          } elsif( exists($wb_oset_names{sprintf("%s_%s", $array->name, $pname)})) {
            $wb_names{sprintf("%s_%s", $array->name, $pname)} = 1;
          }
        }
      }

      if (scalar(keys %wb_names) == 0) {
        warn sprintf("Could not determine any WormBase names for probe id %d\n", $probe->dbID);
        next;
      } elsif (scalar(keys %wb_names) > 1) {
        warn "Multiple WormBase names for probe id ", $probe->dbID, "\n";
      }

      my @pf = @{$probe->get_all_ProbeFeatures};

      next if not @pf;
      
      if ($choose_best) {
        my @clust = map { [$_] } @pf;
        my ($example_wb_pname) = keys %wb_names;
        my $best_clust =  &choose_best_cluster($example_wb_pname, \@clust);
        @pf = ($best_clust->[0]);
      }

      foreach my $pf (@pf) {
        my @blocks = &get_feature_blocks($pf);
        
        my $chr_name = $pf->slice->seq_region_name;
        my $chr_strand = $pf->strand;
        
        foreach my $wb_pname (sort keys %wb_names) {
          if (not exists $done_wb_probes{$wb_pname}) {
            push @{$results{$chr_name}->{$wb_pname}->{$chr_strand}}, \@blocks;
          }
          $done_wb_probes{$wb_pname} = 1;
        }
      }      
    }
  } else {
    foreach my $array (@{$arrays_by_vendor{$vendor}}) {
      # we want to choose a globally best location for each probe based on 
      # co-localty across all probes in the set. 
      
      my @probesets = @{$array->get_all_ProbeSets};
      
      foreach my $probeset (@probesets) {
        my (%feats_by_chr, @clusts);
        
        my $pname = $probeset->name;
        
        foreach my $probe (@{$probeset->get_all_Probes}) {
          foreach my $pf (@{$probe->get_all_ProbeFeatures}) {
            push @{$feats_by_chr{$pf->slice->seq_region_name}}, $pf;
          }
        }
        
        foreach my $chr (keys %feats_by_chr) {
          foreach my $f (sort { $a->start <=> $b->start } @{$feats_by_chr{$chr}}) {
            if (not @clusts or 
                $f->slice->seq_region_name ne $clusts[-1]->[-1]->slice->seq_region_name or
                $f->start - $clusts[-1]->[-1]->start > 10000) {
              push @clusts, [];
            }
            push @{$clusts[-1]}, $f;
          }
        }
        
        next if not @clusts;
        
        if ($choose_best) {
          my $best_clust = &choose_best_cluster($pname, \@clusts);
          @clusts = ($best_clust);
        }
        
        #
        # collapse the probe feature mappings for a probe_set down into 
        # a set of non-overlapping segments. This will make us us lose coordinates
        # of the individual probes (both genome and target), so we will fudge the
        # the target coorindates to span the same region as the genome segments
        #
        foreach my $cl (@clusts) {
          my (@blocks, @nr_segs);
          foreach my $f (@$cl) {
            push @blocks, &get_feature_blocks($f);
          }
          
          foreach my $block (sort { $a->{g_st} <=> $b->{g_st} } @blocks) {
            if (not @nr_segs or $block->{g_st} > $nr_segs[-1]->{g_en} + 1) {
              push @nr_segs, $block;
            } elsif ($block->{g_en} > $nr_segs[-1]->{g_en}) {
              $nr_segs[-1]->{g_en} = $block->{g_en};
            }
          }
          
          if ($cl->[0]->strand < 0) {
            @nr_segs = reverse @nr_segs;
          }
          
          my $count = 1;
          foreach my $seg (@nr_segs) {
            my $len = $seg->{g_en} - $seg->{g_st} + 1;
            $seg->{p_st} = $count;
            $seg->{p_en} = $count + $len - 1;
            $count += $len;
          }
          
          my $chr_name = $cl->[0]->slice->seq_region_name;
          my $chr_strand = $cl->[0]->strand;
          
          push @{$results{$chr_name}->{$pname}->{$chr_strand}}, \@nr_segs;
        }
        
      }
    }
  }
  
  foreach my $chr (sort keys %results) {
    my $chr_len = $sl_adap->fetch_by_region('toplevel', $chr)->length;
    
    &write_ace($vendor, $chr, $chr_len, $results{$chr});
    
  }
}


#########################
sub write_ace {
  my ($platform, $chr, $chr_len, $res_hash) = @_;

  if ($species =~ 'elegans' and grep { $_ eq $chr } ('I', 'II', 'III', 'IV', 'V', 'X', 'MtDNA')) {
    $chr = "CHROMOSOME_${chr}";
  } elsif ($species =~ 'briggsae' and grep { $_ eq $chr } ('I', 'II', 'III', 'IV', 'V', 'X')) {
    $chr = "chr${chr}";
  }

  my (@batches, @homol_data);
  
  foreach my $id (sort keys %$res_hash) {
    if (not @batches or scalar(@{$batches[-1]}) >= 5000) {
      push @batches, [];
    }
    push @{$batches[-1]}, $id;
  }
  
  for(my $idx = 0; $idx < @batches; $idx++) {
    my $tile = "";
    
    if (@batches > 1) {
      $tile .= "_";
      $tile .= $idx + 1;
    }
    
    my $homol_data_name = "Oligo_set:$chr:${platform}${tile}";
    
    print $ace_fh "\nHomol_data : \"$homol_data_name\"\n";
    foreach my $pid (@{$batches[$idx]}) {
      foreach my $strand (keys %{$res_hash->{$pid}}) {
        foreach my $seggroup (@{$res_hash->{$pid}->{$strand}}) {
          foreach my $seg (@$seggroup) {
            
            my ($gst, $gend) = ($seg->{g_st}, $seg->{g_en});
            my ($pst, $pend) = ($seg->{p_st}, $seg->{p_en});
            if ($gend > $chr_len) {
              warn sprintf("Segment extends beyond end of %s : %d %d (%s %d %d), so clipping\n", 
                           $chr, $seg->{g_st}, $seg->{g_en}, $pid, $seg->{p_st}, $seg->{p_en}); 
              
              my $overhang = $gend - $chr_len;
              $gend -= $overhang;
              if ($strand < 0) {
                $pst += $overhang;
              } else {
                $pend -= $overhang;
              }
            }
            if ($gst < 1) {
              warn sprintf("Segment extends beyond start of %s : %d %d (%s %d %d), so clipping\n", 
                           $chr, $seg->{g_st}, $seg->{g_en}, $pid, $seg->{p_st}, $seg->{p_en}); 
              my $overhang = abs($gst) + 1;
              $gst += $overhang;
              if ($strand < 0) {
                $pend -= $overhang;
              } else {
                $pst += $overhang;
              }
            }
            printf($ace_fh "Oligo_set_homol %s Oligo_set_mapping 100.0 %d %d %d %d\n", 
                   $pid, 
                   $strand < 0 ? $gend : $gst, 
                   $strand < 0 ? $gst : $gend,
                   $pst, 
                   $pend);
            
          }
        }
      }
    }
    
    push @homol_data, $homol_data_name;
  }
  
  print $ace_fh "\nSequence : \"$chr\"\n";
  foreach my $hd (@homol_data) {
    print $ace_fh "Homol_data $hd 1 $chr_len\n";
  }
}


#########################
sub choose_best_cluster {
  my ($pn_name, $clust_list) = @_;

  if (scalar(@$clust_list) == 1) {
    #print "$pn_name: Mapped to only has one memberegion to start with\n";
    return $clust_list->[0];
  }

  my @clusts = sort { 
    scalar(@$b) <=> scalar(@$a) or
        $a->[0]->slice->seq_region_name cmp $b->[0]->slice->seq_region_name or
        $a->[0]->start <=> $b->[0]->start
  } @$clust_list;

  #
  # first, cull the clusters that have fewer probes matched than the one(s) with the most
  #
  @clusts = grep { scalar(@$_) == scalar(@{$clusts[0]}) } @clusts;

  if (scalar(@clusts) == 1) {
    #print "$pn_name: Mapped to one region that had more features than the others\n";
    return $clusts[0];
  }

  #
  # Next, if the probe(set) was previously associated with a gene, try to 
  # associate it with the same gene again
  #
  if (exists $p2g_map{$pn_name}) {
    my $matching_clust;

    CLUST: foreach my $clust (@clusts) {
    
      my $sl = $sl_adap->fetch_by_region('toplevel',
                                         $clust->[0]->slice->seq_region_name,
                                         $clust->[0]->start,
                                         $clust->[0]->end);
      
      foreach my $g (@{$sl->get_all_Genes}) {
        if (exists $p2g_map{$pn_name}->{$g->stable_id}) {
          $matching_clust = $clust;
          last CLUST;
        }
      }
    }

    if (defined $matching_clust) {
      #print "$pn_name: Mapped to a region with a gene that was previously associated by Caltech\n";
      return $matching_clust;
    }
  }


  #
  # Finally, sort by:
  #
  # - highest coverage of the probe
  # - shortest span
  # - number of mismatches
  my @new_clust;
  
  foreach my $clust (@clusts) {
    my ($min_st, $max_en, $probe_cov, $mm_count);
    
    foreach my $pf (@{$clust}) {
      $min_st    = $pf->start if not defined $min_st or $pf->start < $min_st;
      $max_en    = $pf->end   if not defined $max_en or $pf->end > $max_en;
      $mm_count += $pf->mismatchcount;
      
      my @blocks = &get_feature_blocks($pf);
      map { $probe_cov += $_->{p_cov} } @blocks;
      
    }
    push @new_clust, [$clust, $probe_cov, $max_en - $min_st + 1, $mm_count];
  }
  
  @new_clust = sort {
    $b->[1] <=> $a->[1] or
        $a->[2] <=> $b->[2] or
        $a->[3] <=> $b->[3]
  } @new_clust;

    #
  # Otherwise, we just have to return the first one
  #

  return $new_clust[0]->[0];
    
}

##########################
sub get_feature_blocks {
  my ($f) = @_;

  my $fstart = $f->start;
  my $fend = $f->end;
  my $fstrand = $f->strand;
  my $cigar = $f->cigar_string;

  my (@new_blocks, @merged_blocks);
  
  my $g_st = $fstart;
  my $p_st = 1;
  
  my @blocks = ($cigar =~ /(\d+[=DIXS])/g);
    
  foreach my $block (@blocks) {
    my ($num, $sym) = $block =~ /^(\d+)(\S)$/;
    
    if ($sym eq '=' or $sym eq 'X' or $sym eq 'S') {
      my $g_en = $g_st + $num - 1;
      my $p_en = $p_st + $num - 1;
      push @new_blocks, {
        g_st => $g_st,
        g_en => $g_en,
        p_st => $p_st,
        p_en => $p_en,
        p_cov => ($p_en - $p_st + 1),
        mm => ($sym eq '=' or $sym eq 'S') ? 0 : $num,
      };

      $g_st = $g_en + 1;
      $p_st = $p_en + 1;
    } elsif ($sym eq 'D') {
      # deletion from probe; insertion (intron) in genome
      $g_st += $num;
    } elsif ($sym eq 'I') {
      # inserion into probe; deletion (skip-over) in genome
      $p_st += $num;
    }
  }

  if ($new_blocks[-1]->{g_en} != $fend) {
    die("Expected $fend, but calculated ". $new_blocks[-1]->{g_en} . " ($fstart $fend $fstrand $cigar\n");
  }


  foreach my $block (@new_blocks) {
    if (not @merged_blocks or $block->{g_st} > $merged_blocks[-1]->{g_en} + 1) {
      push @merged_blocks, $block;
    } else {
      $merged_blocks[-1]->{g_en} = $block->{g_en};
      $merged_blocks[-1]->{p_en} = $block->{p_en};
      $merged_blocks[-1]->{p_cov} += $block->{p_cov};
      $merged_blocks[-1]->{mm} += $block->{mm};
    }
  }

  if ($fstrand < 0) {
    my $plen = $merged_blocks[-1]->{p_en};
    foreach my $bl (@merged_blocks) {
      my $n_st = $plen - $bl->{p_en} + 1;
      my $n_en = $plen - $bl->{p_st} + 1;
      $bl->{p_st} = $n_st;
      $bl->{p_en} = $n_en;
    }
  }

  return @merged_blocks;
}


###########################
sub parse_probe_to_gene_map {
  my ($file, $hash) = @_;

  open(my $fh, $file) or die("Could not open $file for reading\n");
  while(<$fh>) {
    /^(\S+)\s+(\S+)/ and $hash->{$1}->{$2} = 1;
  }
}

#########################
sub parse_probename_file {
  my ($file, $hash) = @_;

  open(my $fh, $file) or die("Could not open $file for reading\n");
  while(<$fh>) {
    /^(\S+)/ and $hash->{$1} = 1;
  }

}
