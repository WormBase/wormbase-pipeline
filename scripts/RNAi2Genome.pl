#!/software/bin/perl -w

# Version: $Version: $
# Last updated by: $Author: klh $
# Last updated on: $Date: 2014-05-22 16:29:31 $

use strict;
use warnings;
use Getopt::Long;

use lib $ENV{CVS_DIR};

use Wormbase;
use Log_files;
use Modules::WormSlurm;
use Ace;
use Coords_converter;
use Sequence_extract;
use Bio::Tools::Blat;
use Bio::SeqIO;

my $BLAT = "blat";
my $BLAT_OPTIONS = "-minIdentity=80 -maxIntron=10000 -noHead";

my $PRIMARY_LEN = 100;
my $PRIMARY_QUAL = 95;
my $SECONDARY_LEN = 200;
my $SECONDARY_QUAL = 80;
my $MAX_HITS = 10;
my $CHECK_BLOCK_SIZE = 1;

################################
# command-line options         #
################################

my ($debug,
    $test,
    $store, 
    $db,
    $workdir,
    $query,
    $target, 
    $acefile,
    $pslfile,
    $noload,
    $species, 
    $no_load,
    $keep_best,
    $verbose,
    $only_unmapped,
    );

GetOptions(
  "debug=s"        => \$debug,
  "test"           => \$test,
  "species:s"      => \$species,
  'store=s'        => \$store,
  "workdir=s"      => \$workdir,
  "database:s"     => \$db,
  "target:s"       => \$target,
  "query:s"        => \$query,
  "acefile:s"      => \$acefile,
  "pslfile:s"      => \$pslfile,
  "noload"         => \$no_load,
  "keepbest=s"     => \$keep_best,
  "verbose"        => \$verbose,
  "onlyunmapped"   => \$only_unmapped,
  "blatexe"        => \$BLAT,
    );


############################
# recreate configuration   #
############################
my $wormbase;
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("cant restore wormbase from $store\n"); 
}
else { 
  $wormbase = Wormbase->new( -debug => $debug, 
                             -test => $test, 
                             -organism => $species); 
}

$db = $wormbase->autoace if not $db;

my $log = Log_files->make_build_log($wormbase);
###########################################

$target = $wormbase->genome_seq if not $target;
$db = $wormbase->autoace if not $db;
$acefile = $wormbase->acefiles . "/RNAi_homols.ace" if not $acefile;
$workdir = $wormbase->blat if not defined $workdir;
$pslfile = "$workdir/RNAi_homols.psl" if not $pslfile;

my $coords = Sequence_extract->invoke($db, 0, $wormbase);

if (not $query) {
  ###### Fetch the queries (if not supplied)
  my (@rnai, @query, @rnai_chunks, @out_ace_files, @out_slurm_files, @out_psl_files, %best_hit_only, %seen);

  @rnai = &get_rnai_sequences();

  my $chunk_size = 2000;      
  my $total_rnai = scalar(@rnai);
  my $chunk_count = int(scalar(@rnai) / $chunk_size);
  
  if (not defined $keep_best) {
    $keep_best =  "$workdir/rnai.best_hit_only.txt";
  }
  open my $pcf, ">$keep_best" 
      or $log->log_and_die("Could not open $keep_best for writing\n");

  # randomise order for better distribution
  while(my $rnai = shift @rnai) {
    my $chunk_id = int(rand($chunk_count));
    
    push @{$rnai_chunks[$chunk_id]}, $rnai;
    if ($rnai->[1] =~ /^mv/ or $rnai->[1] =~ /^yk/) {
      print $pcf $rnai->[0], "\n";
    }
  }
  close($pcf);
  
  $log->write_to(sprintf("Fetched %d RNAi objects in %d chunks\n", $total_rnai, scalar(@rnai_chunks)));
  
  for(my $i=0; $i < @rnai_chunks; $i++) {
    my $file = "$workdir/rnai_query.$i.fa";
    open my $fh, ">$file" or $log->log_and_die("Could not open $file for writing\n");
    foreach my $rnai (@{$rnai_chunks[$i]}) {
      printf $fh ">%s\n%s\n", $rnai->[0], $rnai->[2];
    }
    close($fh);
    push @query, $file;
  }

  my %slurm_jobs;
  for(my $i=0; $i < @query; $i++) {
      my $query = $query[$i];
      my $jname = "worm_rna2genome.$i.$$";
      my $out_file = "$workdir/rnai_out.$i.ace";
      my $psl = "$workdir/rnai_out.$i.psl";
      my $slurm_out = "$workdir/rnai_out.$i.slurm_out";
      my $slurm_err = "$workdir/rnai_out.$i.slurm_err";
      
      my $cmd = ($store) 
	  ? $wormbase->build_cmd_line("RNAi2Genome.pl -query $query -keepbest $keep_best -acefile $out_file -pslfile $psl", $store)
	  : $wormbase->build_cmd("RNAi2Genome.pl -query $query -keepbest $keep_best -acefile $out_file -pslfile $psl");

      my $job_id = WormSlurm::submit_job_with_name($cmd, 'production', '2600m', '06:00:00', $slurm_out, $slurm_err, $jname);
      $slurm_jobs{$job_id} = $cmd;

      push @out_ace_files, $out_file;
      push @out_psl_files, $psl;
      push @out_slurm_files, $slurm_out;
  }     
    
  WormSlurm::wait_for_jobs(%slurm_jobs);
  
  $log->write_to("All RNAi2Genome batch jobs have completed!\n");
  for my $job_id (keys %slurm_jobs) {
      $log->error("Slurm job $job_id (" . $slurm_jobs{$job_id} . ") exited non zero\n") if WormSlurm::get_exit_code($job_id) != 0;
  }
  
  my (%parent_seqs);
  open(my $out_fh, ">$acefile") or $log->log_and_die("Could not open $acefile for writing\n");
  foreach my $ace (@out_ace_files) {
    open my $this_fh, $ace or $log->log_and_die("Could not open $ace for reading\n"); 
    while(<$this_fh>) {
      print $out_fh $_;
      if (/^Homol_data\s+:\s+\"(\S+):RNAi\"/) {
        $parent_seqs{$1} = 1;
      }
    }
  }

  foreach my $parent (sort keys %parent_seqs) {
    printf $out_fh "\nSequence : \"%s\"\n", $parent;
    printf $out_fh "Homol_data %s:RNAi %d %d\n", $parent, 1, $coords->Superlink_length($parent); 
  }
  close($out_fh);
     
  if (not $log->report_errors and not $test) {
    unlink @out_ace_files, @out_psl_files, @out_slurm_files, @query;
  }

  if (not $no_load) {
    $wormbase->load_to_database($db, $acefile, "RNAi_mapping", $log);
  }
} else {

  my (%keep_best_clones, %qlengths, $primary, $secondary, %results, %results_by_parent, %has_hit);
  
  if ($keep_best) {
    open my $pf, $keep_best or $log->log_and_die("Could not open $keep_best\n");
    while(<$pf>) {
      /^(\S+)/ and $keep_best_clones{$1} = 1;
    }
  }
  
  my $seqio = Bio::SeqIO->new(-format => 'fasta',
                              -file => $query);
  while(my $seq = $seqio->next_seq) {
    $qlengths{$seq->id} = $seq->length;
  }
  
  my $command = "$BLAT $target $query $pslfile $BLAT_OPTIONS";
  
  print STDERR "$command\n";
  system("$command") and $log->error("Command '$command' failed\n");
  
  open my $out_ace, ">$acefile";

  $log->write_to("Processing PSL file $pslfile...\n");
    
  open(my $fh, "<$pslfile") or $log->log_and_die("cant open blat output $pslfile : $!\n");
  my $parser = Bio::Tools::Blat->new('-fh' => $fh);
  
  {
    my %results;

    while(my $result = $parser->next_result) {
      
      my ($percent) = $result->get_tag_values('percent_id');
      
      my ($total_aligned, $max_block_size) = (0, 0);
      foreach my $sf ($result->get_SeqFeatures) {
        my $len = ($sf->feature1->end - $sf->feature1->start + 1);
        $total_aligned += $len;
        if ($len > $max_block_size) {
          $max_block_size = $len;
        }
      }
      
      # note that filter hits will check that minimum length and quality
      # criteria are met. However, we weed out some stuff here that we know
      # immediately will fail the filter, to save memory.
      if ($keep_best_clones{$result->seq_id}) {
        if ($percent >= $PRIMARY_QUAL and
            ($total_aligned >= $PRIMARY_LEN or $qlengths{$result->seq_id} == $total_aligned) and 
            (not exists $results{$result->seq_id} or
             $total_aligned > $results{$result->seq_id}->[0]->[1])) {
          $results{$result->seq_id} = [[$percent, $total_aligned, $max_block_size, $result]];
        }
      } else {
        if ( ( ($total_aligned == $qlengths{$result->seq_id} or $total_aligned >= $PRIMARY_LEN) and $percent >= $PRIMARY_QUAL) or
             ($total_aligned >= $SECONDARY_LEN and $percent > $SECONDARY_QUAL) ) {
          push @{$results{$result->seq_id}}, [$percent, $total_aligned, $max_block_size, $result];
        }
      }
    }
    
    ($primary, $secondary) = &filter_hits(\%results, \%qlengths, $CHECK_BLOCK_SIZE);
  }

  foreach my $type_pair ([$primary, 'RNAi_primary'], [$secondary, 'RNAi_secondary']) {
    my ($list, $method) = @$type_pair;
    
    while( my $result = shift @$list) {
      
      my @b = $result->get_SeqFeatures;
      
      my ($tid, $qid, $tstart, $tend, $strand) = ($b[0]->seq_id, 
                                                  $b[0]->hseq_id, 
                                                  $result->start, 
                                                  $result->end, 
                                                  $result->strand);
      $has_hit{$method}->{$qid} = 1;
      
      my ($percent) = $result->get_tag_values('percent_id');
      
      my (@blocks, @projected_blocks, %proj_parents);

      foreach my $block (@b) {
        my $seg_tstart = $block->feature1->start;
        my $seg_tend   = $block->feature1->end;
        my $seg_qstart = $block->feature2->start;
        my $seg_qend   = $block->feature2->end;
        
        my ($seg_m_tid, $seg_m_tstart, $seg_m_tend) = $coords->LocateSpan($tid,$seg_tstart,$seg_tend);

        if ($strand < 0) {
          ($seg_tstart, $seg_tend) = ($seg_tend, $seg_tstart);
          ($seg_m_tstart, $seg_m_tend) = ($seg_m_tend, $seg_m_tstart);
        }
        push @blocks, [$seg_tstart, $seg_tend, $seg_qstart, $seg_qend];
        push @projected_blocks, [$seg_m_tstart, $seg_m_tend, $seg_qstart, $seg_qend];
        $proj_parents{$seg_m_tid} = 1;
      }

      $qid =~ s/\.DNA_text_\d+$//;
      $qid =~ s/\.PCR_product_\d+$//;
      
      # only use projected mapping if all segments map down to same parent
      if (scalar(keys %proj_parents) == 1) {
        my ($parent) = keys %proj_parents;
        push @{$results_by_parent{$parent}->{$qid}->{$method}}, [$percent, \@projected_blocks];
      } else {
        push @{$results_by_parent{$tid}->{$qid}->{$method}}, [$percent, \@blocks];
      }
    }
  }
  
  foreach my $parent (keys %results_by_parent) {
    printf $out_ace "\nHomol_data : \"%s:RNAi\"\n", $parent;
    print $out_ace "Sequence $parent\n";
    my $t_res = $results_by_parent{$parent};
    
    foreach my $query (keys %$t_res) {
      my $t_q_res = $t_res->{$query};

      foreach my $meth (keys %$t_q_res) {
        my @hits = @{$t_q_res->{$meth}};
        
        foreach my $hit (@hits) {
          my ($score, $blocks) = @$hit;
          foreach my $block (@$blocks) {
            printf $out_ace "RNAi_homol %s %s %s %d %d %d %s\n", $query, $meth, $score, @$block;
          }
        }  
      }
    }
  }
  
  close($out_ace);

  foreach my $qid (keys %qlengths) {
    if (not exists $has_hit{RNAi_primary}->{$qid} and
        not exists $has_hit{RNAi_secondary}->{$qid}) {
      $log->write_to("Did not find any acceptable hits (primary or secondary) for $qid\n");
    } elsif (not exists $has_hit{RNAi_primary}->{$qid} and
             exists $has_hit{RNA_secondary}->{$qid}) {
      $log->write_to("Only found secondary hits for $qid\n");
    }
  }

}


$log->mail;
exit;



################################################################
sub get_rnai_sequences {

  my (@rnai, %seen, $tb_file, $tb_cmd, $tace_fh);
  
  my $tace = $wormbase->tace;

  #
  # First fetch the objects with DNA_text attached (old Caltech pipeline)
  #
  $tb_file = &query_with_dna_text();
  if (not -e $tb_file) {
    $log->log_and_die("Could not find Table-maker definition file $tb_file\n");
  }
  $tb_cmd = "Table-maker -p \"$tb_file\"\nquit\n";
  open($tace_fh, "echo '$tb_cmd' | $tace $db |");
  while(<$tace_fh>) {
    print ;
    /^\"(\S+)\"/ and do {
      my ($rnai_id) = $1;
      if (/^\"\S+\"\s+\"(\S+)\"\s+\"(\S+)\"/) {
        my ($dna_text, $clone_id) = ($1, $2);
        $seen{$rnai_id}++;
        my $unique_rnai_id = sprintf("%s.DNA_text_%d", $rnai_id,  $seen{$rnai_id});
      
        push @rnai, [$unique_rnai_id, $clone_id, $dna_text];
      } elsif (/^\"\S+\"\s+\"(\S+)\"/) {
        my $dna_text = $1;
        if ($dna_text !~ /^[ACGTNacgtn]+$/) {
          $log->write_to("Will not attempt to align $rnai_id - dodgy chars in DNA_text\n");
        } else {
          $seen{$rnai_id}++;
          my $unique_rnai_id = sprintf("%s.DNA_text_%d", $rnai_id,  $seen{$rnai_id});
      
          push @rnai, [$unique_rnai_id, "NO_CLONE_GIVEN", $dna_text];
        }
      } else {
        $log->write_to("Will not attempt to align $rnai_id - dodgy DNA_text\n");
      }
    };
  }
  close($tace_fh);
  unlink $tb_file;

  #
  # Now get the objects the new Caltech pipeline, which have no DNA_text but a PCR_product instead
  #
  my (%rnai2pcr, %pcr, %pcr_seq);

  $tb_file = &query_with_pcr_product();
  $tb_cmd = "Table-maker -p \"$tb_file\"\nquit\n";

  open($tace_fh, "echo '$tb_cmd' | $tace $db |");
  while(<$tace_fh>) {
    /^\"(\S+)\"\s+\"(\S+)\"/ and do {
      my ($rnai_id, $prod) = ($1, $2);

      if (/^\"\S+\"\s+\"\S+\"\s+\"(\S+)\"\s+\"(\S+)\"\s+(\d+)\s+(\d+)/) {
        my ($seq, $prod2, $start, $end) = ($1, $2, $3, $4);

        next unless $prod eq $prod2;
        if (not exists $pcr{$prod}) {
          $pcr{$prod} = [$seq, $start, $end];
        }
        $rnai2pcr{$rnai_id}->{$prod} = 1;
      } else {
        $log->write_to("Will not attempt to align $rnai_id - No DNA_text, and attached PCR_product $prod is not mapped\n");
      }
    }
  }
  close($tace_fh);
  unlink $tb_file;
  
  # Fetch underlying sequences
  foreach my $prod (keys %pcr) {
    my ($seq, $start, $end) = @{$pcr{$prod}};
    ($start, $end) = ($end, $start) if $end < $start;

    my $seq_string = $coords->Sub_sequence($seq, $start - 1, $end - $start + 1);

    $pcr_seq{$prod} = $seq_string;
  }

  foreach my $rnai_id (keys %rnai2pcr) {
    foreach my $pcr_prod (keys %{$rnai2pcr{$rnai_id}}) {
      if (exists $pcr_seq{$pcr_prod}) {
        $seen{$rnai_id}++;
        my $unique_rnai_id = sprintf("%s.PCR_product_%d", $rnai_id,  $seen{$rnai_id});

        push @rnai, [$unique_rnai_id, $pcr_prod, $pcr_seq{$pcr_prod}];
      }
    }
  }

  return @rnai;
}

#########################
sub query_with_dna_text {

  my $tmdef = "/tmp/rnai_with_dnatext.$$.def";

  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $condition = ($only_unmapped) ? "Condition NOT Homol_homol" : "";

  my $query = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class RNAi 
From 1
$condition
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Text 
From 1 
Tag DNA_text 
 
Colonne 3 
Width 12 
Optional 
Visible 
Text 
Right_of 2 
Tag  HERE  

EOF

  print $qfh $query;
  close($qfh);
  return $tmdef;
}

#########################
sub query_with_pcr_product {

  my $tmdef = "/tmp/rnai_wit_pcrprod.$$.def";

  open my $qfh, ">$tmdef" or 
      $log->log_and_die("Could not open $tmdef for writing\n");  

  my $condition = ($only_unmapped) ? "Condition NOT Homol_homol" : "";

  my $query = <<"EOF";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class RNAi 
From 1 
$condition
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class PCR_product 
From 1 
Tag PCR_product 
 
Colonne 3 
Width 12 
Optional 
Visible 
Class 
Class Sequence 
From 2 
Tag Canonical_parent 
 
Colonne 4 
Width 12 
Optional 
Visible 
Class 
Class PCR_product 
From 3 
Tag PCR_product 
 
Colonne 5 
Width 12 
Optional 
Visible 
Integer 
Right_of 4 
Tag  HERE  
 
Colonne 6 
Width 12 
Optional 
Visible 
Integer 
Right_of 5 
Tag  HERE  
 
 
EOF

  print $qfh $query;
  close($qfh);
  return $tmdef;
}


##########################
sub filter_hits {
  my ($res, $qlengths, $check_block_size) = @_;

  my (@all_primary, @all_secondary);

  foreach my $qid (keys %$res) {

    my (@primary, @secondary);

    foreach my $quad (@{$res->{$qid}}) {
      my ($percent, $total_aligned, $max_block_size, $hit) = @$quad;

      if ( ($total_aligned eq $qlengths->{$qid} or $total_aligned >= $PRIMARY_LEN) and $percent >= $PRIMARY_QUAL) {
        push @primary, $quad;
      } elsif ($total_aligned >= $SECONDARY_LEN and $percent >= $SECONDARY_QUAL) {
        push @secondary, $quad;
      }
    }
    
    @primary = sort { $b->[1] <=> $a->[1] } @primary;
    @secondary = sort { $b->[1] <=> $a->[1] } @secondary;

    if (scalar(@primary) + scalar(@secondary) > $MAX_HITS) {
      if (scalar(@primary) >= $MAX_HITS) {
        @secondary = ();
        splice(@primary, $MAX_HITS);
      } else {
        splice(@secondary, $MAX_HITS - scalar(@primary));
      }
    }

    if (defined $check_block_size and scalar(@primary) + scalar(@secondary) > 1) {
      # always keep the best primary hit, regardless
      if (@primary) {
        my ($best, @others) = @primary;
        @primary = ($best, grep { $_->[2] >= $PRIMARY_LEN } @others);
      }
      @secondary = grep { $_->[2] >= $SECONDARY_LEN } @secondary;
    }

    push @all_primary, map { $_->[3] } @primary;
    push @all_secondary, map { $_->[3] } @secondary;
  }

  return (\@all_primary, \@all_secondary);
}

#############################
sub log_blat {
  my ($perc, $aligned, $max_bsize, $res) = @_;

  my @b = $res->get_SeqFeatures;
  $log->write_to(sprintf("RESULT: %s %d-%d %s %s %d %d\n", 
                         $b[0]->feature1->seq_id, 
                         $res->start, 
                         $res->end, 
                         $res->seq_id, 
                         $perc, 
                         $aligned, 
                         $max_bsize));
  foreach my $b (@b) {
    $log->write_to(sprintf(" SEG: %d-%d  %d-%d (%d)\n", 
                           $b->feature1->start, 
                           $b->feature1->end, 
                           $b->feature2->start, 
                           $b->feature2->end, 
                           ($b->feature1->end - $b->feature1->start + 1)));
  }
}
