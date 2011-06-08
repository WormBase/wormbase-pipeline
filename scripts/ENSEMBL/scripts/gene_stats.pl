#! /usr/local/ensembl/bin/perl -w

# hacked from my original gene_stats_otter.pl

# gene_stats.pl : calculates basic gene stats.

use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Getopt::Long;

use strict;

my (@dbname,
    $dbhost,
    $dbport, 
    $dbuser,
    $pep_stats,
    $simple,
    $only_biotype, 
    $list_size);

my $strand = 1;

&GetOptions( 
             'dbhost=s' =>  \$dbhost,
	     'dbport=s' => \$dbport,
             'dbuser=s' => \$dbuser,
             'listsize=s' => \$list_size,
             'biotype=s'  => \$only_biotype,
             'simple'    => \$simple,
	     'pepstats' => \$pep_stats);

die "You must give a dbhost\n" unless $dbhost;

$list_size = 1 if not $list_size;

foreach my $dbname (@ARGV) {
  my $res = {
    title => $dbname,
  };


  my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new( -host => $dbhost,
                                                -dbname => $dbname,
                                                -port => $dbport,
                                                -user => "ensro");
  # get gene adaptor
  my $g_ad = $db->get_GeneAdaptor;
  
  my @slices = @{$db->get_SliceAdaptor->fetch_all('toplevel')};
  
  my %all_genes;
  foreach my $slice (@slices) {
    my ($chr_name, $chr_start, $chr_end) = ($slice->seq_region_name, $slice->start, $slice->end);
    
    print STDERR "Doing segment $chr_name start $chr_start end $chr_end\n";
    
    my $slice_for_genes = $db->get_SliceAdaptor->fetch_by_region('chromosome', 
                                                                 $chr_name);
    
    print STDERR "  Fetching genes for slice ".$slice_for_genes->name."...\n";
    
    foreach my $g (@{$g_ad->fetch_all_by_Slice($slice_for_genes, undef, 1)}) {
      push @{$all_genes{$g->biotype}}, $g;
    }
  }
  
  
  foreach my $biotype (keys %all_genes) {
    next if defined $only_biotype and $only_biotype ne $biotype;

    my @this_gene_list = @{$all_genes{$biotype}};

    foreach my $gene (@this_gene_list) {
      
      my $gene_name = $gene->stable_id;
      my $gene_strand = $gene->strand;
      
      my $sten = &get_gene_start_ends( $gene );
      push @{$res->{gene_start_ends}}, $sten;
      my $gene_len =  $sten->[1] - $sten->[0] + 1;
      
      if (not exists($res->{longest_gene}) or
          @{$res->{longest_gene}} < $list_size or
          $res->{longest_gene}->[-1]->[1] < $gene_len) {
        
        # insert into correct place
        push @{$res->{longest_gene}}, [$gene_name, $gene_len];
        $res->{longest_gene} = [sort {$b->[1] <=> $a->[1]} @{$res->{longest_gene}}];
        if (@{$res->{longest_gene}} > $list_size) {
          pop @{$res->{longest_gene}};
        }
        
      }
      if (not exists($res->{shortest_gene}) or
          @{$res->{shortest_gene}} < $list_size or
          $res->{shortest_gene}->[-1]->[1] > $gene_len) {
      
        # insert into correct place
        push @{$res->{shortest_gene}}, [$gene_name, $gene_len];
        $res->{shortest_gene} = [sort {$a->[1] <=> $b->[1]} @{$res->{shortest_gene}}];
        if (@{$res->{shortest_gene}} > $list_size) {
          pop @{$res->{shortest_gene}};
        }
        
      }
      
      my @trans = @{$gene->get_all_Transcripts};
      
      if (not exists($res->{most_transcripts}) or
          @{$res->{most_transcripts}} < $list_size or
          $res->{most_transcripts}->[-1]->[1] < scalar(@trans)) {
        
        # insert into correct place
        push @{$res->{most_transcripts}}, [$gene_name, scalar(@trans)];
        $res->{most_transcripts} = [sort {$b->[1] <=> $a->[1]} @{$res->{most_transcripts}}];
        if (@{$res->{most_transcripts}} > $list_size) {
          pop @{$res->{most_transcripts}};
        }        
      }
      
      push @{$res->{transcript_counts}}, scalar(@trans);
      
      foreach my $tran (@trans) {
        
        my $t_name = $tran->stable_id;
        
        my @exons = sort {$a->start <=> $b->start} @{$tran->get_all_Exons}; 
        my $num_exons = scalar(@exons);
        
        #
        # Record longest and shortest exons
        #
        
        for (my $i = 0; $i < @exons; $i++ ) {		
          my $ex_num = ($gene_strand < 0) ? $num_exons - $i: $i+1; 
          my $ex_len = $exons[$i]->end - $exons[$i]->start + 1;
          
          if (not exists($res->{longest_exon}) or
              @{$res->{longest_exon}} < $list_size or
              $res->{longest_exon}->[-1]->[2] < $ex_len) {
            
            push @{$res->{longest_exon}}, [$t_name, $ex_num, $ex_len];
            $res->{longest_exon} = [sort {$b->[2] <=> $a->[2]} @{$res->{longest_exon}}];
            if (@{$res->{longest_exon}} > $list_size) {
              pop @{$res->{longest_exon}};
            }
          }
          
          if ($pep_stats) {
            if ($exons[$i]->phase != -1 or $exons[$i]->end_phase != -1) {
              if (not exists($res->{longest_coding_exon}) or
                  @{$res->{longest_coding_exon}} < $list_size or
                  $res->{longest_coding_exon}->[-1]->[2] < $ex_len) {
                
                push @{$res->{longest_coding_exon}}, [$t_name, $ex_num, $ex_len];
                $res->{longest_coding_exon} = [sort {$b->[2] <=> $a->[2]} @{$res->{longest_coding_exon}}];
                if (@{$res->{longest_coding_exon}} > $list_size) {
                  pop @{$res->{longest_coding_exon}};
                }
              }
            }
          }
          
          # dont use initial or terminal exons to find shortest exon
          if ($i != 0 and $i != $#exons) {
            
            if (not exists($res->{shortest_exon}) or
                @{$res->{shortest_exon}} < $list_size or
                $res->{shortest_exon}->[-1]->[2] > $ex_len) {
              
              push @{$res->{shortest_exon}}, [$t_name, $ex_num, $ex_len];
              $res->{shortest_exon} = [sort {$a->[2] <=> $b->[2]} @{$res->{shortest_exon}}];
              if (@{$res->{shortest_exon}} > $list_size) {
                pop @{$res->{shortest_exon}};
              }
            }
            
            if ($pep_stats) {
              if ($exons[$i]->phase != -1 or $exons[$i]->end_phase != -1) {
                if (not exists($res->{shortest_coding_exon}) or
                    @{$res->{shortest_coding_exon}} < $list_size or
                    $res->{shortest_coding_exon}->[-1]->[2] < $ex_len) {
                  
                  push @{$res->{shortest_coding_exon}}, [$t_name, $ex_num, $ex_len];
                  $res->{shortest_coding_exon} = [sort {$a->[2] <=> $b->[2]} @{$res->{shortest_coding_exon}}];
                  if (@{$res->{shortest_coding_exon}} > $list_size) {
                    pop @{$res->{shortest_coding_exon}};
                  }
                }
              }
            }
          }
        }
        
        #
        # Record longest and shortest introns
        #
        for (my $i = 0; $i < @exons - 1; $i++) {
          my $ex_num = sprintf("%d/%d", 
                               ($gene_strand < 0) ? $num_exons - ($i+1): $i+1,
                               ($gene_strand < 0) ? $num_exons - $i : $i+2 );; 
          
          my $intron_start = $exons[$i]->end + 1;
          my $intron_end = $exons[$i+1]->start - 1;
          my $intron_len = $intron_end - $intron_start + 1;
          
          if (not exists($res->{longest_intron}) or
              @{$res->{longest_intron}} < $list_size or
              $res->{longest_intron}->[-1]->[2] < $intron_len) {
            
            push @{$res->{longest_intron}}, [$t_name, $ex_num, $intron_len];
            $res->{longest_intron} = [sort {$b->[2] <=> $a->[2]} @{$res->{longest_intron}}];
            if (@{$res->{longest_intron}} > $list_size) {
              pop @{$res->{longest_intron}};
            }
          }
          if (not exists($res->{shortest_intron}) or
              @{$res->{shortest_intron}} < $list_size or
              $res->{shortest_intron}->[-1]->[2] > $intron_len) {
            
            push @{$res->{shortest_intron}}, [$t_name, $ex_num, $intron_len];
            $res->{shortest_intron} = [sort {$a->[2] <=> $b->[2]} @{$res->{shortest_intron}}];
            if (@{$res->{shortest_intron}} > $list_size) {
              pop @{$res->{shortest_intron}};
            }
          }
        }
        
        if ($pep_stats) {
          my ($coding_trans, $complete_cds_5, $complete_cds_3) = (0, 0, 0);
          my ($complete_utr_3, $complete_utr_5, $complete_utr_both) = (0, 0, 0);
          
          my @peptides;
          
          foreach my $tran (@trans) {
            
            my $translation = $tran->translation;
            
            if (defined($translation)) {
              push @peptides, $tran->translate->seq;
              $coding_trans++;
            }		    
          }
          
          push @{$res->{coding_transcript_counts}}, $coding_trans;
        }
      }
      
      # remaining stats are based upon the longest transcript;
      my $longest_tran = &get_gene_longest_transcript($gene);
      my $t_name = $longest_tran->stable_id;
      my @exons = sort {$a->start <=> $b->start} @{$longest_tran->get_all_Exons}; 
      
      push @{$res->{exon_counts}}, scalar(@exons);
      
      if (not exists($res->{most_exons}) or
          @{$res->{most_exons}} < $list_size or
          $res->{most_exons}->[-1]->[1] < scalar(@exons)) {
        
        push @{$res->{most_exons}}, [$t_name, scalar(@exons)];
        $res->{most_exons} = [sort {$b->[1] <=> $a->[1]} @{$res->{most_exons}}];
        if (@{$res->{most_exons}} > $list_size) {
          pop @{$res->{most_exons}};
        }
      }
      
      my $tlen = 0;
      for (my $i = 0; $i < @exons; $i++) {
        my ($st, $en) = ($exons[$i]->start, $exons[$i]->end);
        if ($en > $st) {
          push @{$res->{exon_start_ends}}, [$st, $en];
          $tlen += ($en - $st + 1);
        }
      }
      push @{$res->{spliced_lengths}},  [1, $tlen];
      
      if (@exons == 1) {
        $res->{single_exon_genes}++;
      }
      else {
        for (my $i = 0; $i < @exons - 1; $i++) {
          my ($st, $en) = ($exons[$i]->end + 1,$exons[$i+1]->start - 1);		    
          
          if ($en > $st) {
            push @{$res->{intron_start_ends}}, [$st, $en];
          }
          
          if ($exons[$i]->end_phase != -1) {
            $res->{intron_phases}->{$exons[$i]->end_phase}++;
          }		
        }
      }
    }   
    if ($simple) {
      &process_results_simple($res);
    } else {
      &process_results($res);
    }
  }
}



sub get_gene_start_ends {
  my $g = shift;
  
  my ($st, $en);
  foreach my $e (@{$g->get_all_Exons}) {
    if (not defined $st or $e->start < $st) {
      $st = $e->start;
    } 
    if (not defined $en or $e->end > $en) {
      $en = $e->end;
    }
  }
  
  return [$st, $en];
}
	

sub get_gene_longest_transcript {
  my $g = shift;
  
  my $longest;
  foreach my $t (@{$g->get_all_Transcripts()}) {
    if (not defined $longest or $longest->length < $t->length) {
      $longest = $t;
    }
  }
  
  return $longest;
}


sub project_start_ends {
  my $start_ends = shift;
  
  my @newlist;
  foreach my $sten (sort {$a->[0] <=> $b->[0]} @$start_ends) {
    my $this_sten = [$sten->[0], $sten->[1]];
    
    if (not @newlist or $newlist[-1]->[1] < $this_sten->[0]) {
      push @newlist, $this_sten;
    }
    else {
      if ($this_sten->[1] > $newlist[-1]->[1]) {
        $newlist[-1]->[1] = $this_sten->[1];
      }
    }
  }
  
  my $total = 0;
  foreach my $sten (@newlist) {
    $total += $sten->[1] - $sten->[0] + 1;
  }
  
  return $total;
}


sub calc_mean {
  my @list = @_;
  
  my $tot = 0;
  foreach my $num (@list) {
    $tot += $num;
  }
  
  return $tot / scalar(@list);
}


sub calc_median {
  my @list = @_;
  
  my @sorted = sort {$a <=> $b} @list;
  my $size = scalar(@sorted);
  if ($size % 2) {
    return $sorted[($size - 1)/ 2];
  }
  else {
    my $index = $size / 2;
    return ($sorted[$index - 1] + $sorted[$index]) / 2;
  }
}



sub process_results {
  my $res = shift;
  
  if (not exists ($res->{gene_start_ends})) {
    printf "Total: 0\n";
    return;
  }

  printf "#################\n";
  printf "####### %s ######\n", $res->{title};
  printf "#################\n";

  printf "Total genes : %d\n\n", scalar(@{$res->{gene_start_ends}});
  printf "Mean   length : %d\n", &calc_mean (map { $_->[1] - $_->[0] + 1 } @{$res->{gene_start_ends}});;
  printf "Median length : %d\n", &calc_median (map { $_->[1] - $_->[0] + 1 } @{$res->{gene_start_ends}});;
  printf("Total transcribed bp : %d\n\n", &project_start_ends($res->{gene_start_ends})); 
  
  my $this_total = 0; foreach my $t (@{$res->{transcript_counts}}) { $this_total += $t; }
  printf("Transcripts :  %8s %8s %8s\n", "Total",  "Mean", "Median");
  printf("               %8d %8.1f %8.1f\n",  	   
         $this_total,
         &calc_mean( @{$res->{transcript_counts}} ),
         &calc_median( @{$res->{transcript_counts}} ));
  printf("  Top $list_size genes with most transcripts:\n");
  foreach my $en (@{$res->{most_transcripts}}) {
    printf("    %-30s %d\n", $en->[0], $en->[1]); 
  }
  print "\n";

  if ($pep_stats) {
    $this_total = 0; foreach my $t (@{$res->{coding_transcript_counts}}) { $this_total += $t; }
    printf("Total coding transcripts : %d\n", $this_total);    
    printf("Genes with one more than one transcript : %d\n", scalar(grep { $_ > 1 } @{$res->{transcript_counts}}));

    my $total_complete_cds = 0;
    my $total_different_pep = 0;
    for(my $i = 0; $i < @{$res->{transcript_counts}}; $i++) {
      if ($res->{transcript_counts}->[$i] > 1) {
        if ($res->{complete_cds_transcript_counts}->[$i] > 1) {
          $total_complete_cds++;
        }
        if ($res->{unique_complete_peptide_counts}->[$i] > 1) {
          $total_different_pep++;
        }
      }		
    }
    printf ("   ... with >1 complete CDS : %d\n", $total_complete_cds);
    printf ("      ... with >1 different translations : %d\n\n", $total_different_pep);
  }
	
  printf("LONGEST GENES (top $list_size):\n");
  foreach my $en (@{$res->{longest_gene}}) {
    printf("  %-30s (%d bp)\n", $en->[0], $en->[1]);

  }
  printf("SHORTEST GENES (top $list_size):\n");
  foreach my $en (@{$res->{shortest_gene}}) {
    printf("  %-30s (%d bp)\n", $en->[0], $en->[1]);

  }

  printf("LONGEST EXONS (top $list_size):\n");
  foreach my $en (@{$res->{longest_exon}}) {
    printf("  %-30s (%d bp)\n", $en->[0] . " (ex " . $en->[1] . ") ", $en->[2]);
  }
  printf("SHORTEST EXONS (top $list_size):\n");
  foreach my $en (@{$res->{shortest_exon}}) {
    printf("  %-30s (%d bp)\n", $en->[0] . " (ex " . $en->[1] . ") ", $en->[2]);
  }

  if (exists ($res->{longest_coding_exon})) {
    printf("LONGEST CODING EXONS (top $list_size):\n");
    foreach my $en (@{$res->{longest_coding_exon}}) {
      printf("  %-30s (%d bp)\n", $en->[0] . " (ex " . $en->[1] . ") ", $en->[2]);
    }
  }
  if (exists ($res->{shortest_coding_exon})) {
    printf("SHORTEST CODING EXONS (top $list_size):\n");
    foreach my $en (@{$res->{shortest_coding_exon}}) {
      printf("  %-30s (%d bp)\n", $en->[0] . " (ex " . $en->[1] . ") ", $en->[2]);
    }
  }

  printf("LONGEST INTRONS (top $list_size):\n");
  foreach my $en (@{$res->{longest_intron}}) {
    printf("  %-30s (%d bp)\n", $en->[0] . " (ex " . $en->[1] . ") ", $en->[2]);
  }
  printf("SHORTEST INTRONS (top $list_size):\n");
  foreach my $en (@{$res->{shortest_intron}}) {
    printf("  %-30s (%d bp)\n", $en->[0] . " (ex " . $en->[1] . ") ", $en->[2]);
  }

  printf("\nStats based on the longest transcript from each gene:\n\n");
    
  printf("          %8s %8s %12s %12s\n", "Total", "Total-bp", "Mean-len", "Median-len" );
  if (exists($res->{spliced_lengths})) {
    printf("SPLICED   %8d %8d %12d %12d\n",
           scalar(@{$res->{spliced_lengths}}),
           0,
           &calc_mean (map { $_->[1] - $_->[0] + 1 } @{$res->{spliced_lengths}}),
           &calc_median (map { $_->[1] - $_->[0] + 1 } @{$res->{spliced_lengths}}));
  }

  if (exists($res->{exon_start_ends})) {
    printf("EXONS     %8d %8d %12d %12d\n",
           scalar(@{$res->{exon_start_ends}}),
           &project_start_ends( $res->{exon_start_ends}),
           &calc_mean (map { $_->[1] - $_->[0] + 1 } @{$res->{exon_start_ends}}),
           &calc_median (map { $_->[1] - $_->[0] + 1 } @{$res->{exon_start_ends}}));
  }
  if (exists($res->{intron_start_ends})) {
    printf("INTRONS   %8d %8d %12d %12d\n\n",
           scalar(@{$res->{intron_start_ends}}),
           &project_start_ends( $res->{intron_start_ends}),
           &calc_mean (map { $_->[1] - $_->[0] + 1 } @{$res->{intron_start_ends}}),
           &calc_median (map { $_->[1] - $_->[0] + 1 } @{$res->{intron_start_ends}}));
  }
  printf("Exons per gene :  %8s %8s\n", "Mean", "Median");
  printf("                  %8.1f %8.1f\n",  	   
         &calc_mean( @{$res->{exon_counts}} ),
         &calc_median( @{$res->{exon_counts}}));
  printf("  Top $list_size genes with most exons in longest transcript:\n");
  foreach my $en (@{$res->{most_exons}}) {
    printf("    %-30s (%d exons)\n", $en->[0], $en->[1]);
  }
  
  printf("\nSingle exon genes : %d\n", $res->{single_exon_genes}); 
  foreach my $i ('0', '1', '2') {
    printf("Phase $i introns : %d\n", $res->{intron_phases}->{$i});
  }
  print "\n";
}


sub process_results_simple {
  my $res = shift;
  
  my @list = ($res->{title});

  push @list, scalar(@{$res->{gene_start_ends}});
  my $this_total = 0; foreach my $t (@{$res->{transcript_counts}}) { $this_total += $t; }
  push @list, $this_total;

  push @list, sprintf("%d (%d)", 
                      &calc_mean (map { $_->[1] - $_->[0] + 1 } @{$res->{gene_start_ends}}), 
                      &calc_median (map { $_->[1] - $_->[0] + 1 } @{$res->{gene_start_ends}}));

  push @list, sprintf("%.1f (%.1f)", 
                      &calc_mean( @{$res->{exon_counts}} ),
                      &calc_median( @{$res->{exon_counts}}));

  push @list, sprintf("%d (%d)", 
                      &calc_mean (map { $_->[1] - $_->[0] + 1 } @{$res->{exon_start_ends}}),
                      &calc_median (map { $_->[1] - $_->[0] + 1 } @{$res->{exon_start_ends}}));

  push @list, sprintf("%d (%d)",
                      &calc_mean (map { $_->[1] - $_->[0] + 1 } @{$res->{intron_start_ends}}),
                      &calc_median (map { $_->[1] - $_->[0] + 1 } @{$res->{intron_start_ends}}));


  push @list, sprintf("%d (%d)", 
              &calc_mean (map { $_->[1] - $_->[0] + 1 } @{$res->{spliced_lengths}}),
              &calc_median (map { $_->[1] - $_->[0] + 1 } @{$res->{spliced_lengths}}));

  print join("\t", @list), "\n";
}

