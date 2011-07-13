#!/software/bin/perl -w

# Version: $Version: $
# Last updated by: $Author: klh $
# Last updated on: $Date: 2011-07-13 08:38:42 $

use strict;
use warnings;
use Getopt::Long;

use lib $ENV{'CVS_DIR'};

use Wormbase;
use Log_files;

use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;
use Ace;
use Coords_converter;
use Bio::Tools::Blat;
use Bio::SeqIO;

my $BLAT = "/software/worm/bin/blat/blat";
#my $BLAT = "/nfs/users/nfs_k/klh/bin/x86_64/blat";
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
    $query,
    $target, 
    $submit,
    $acefile,
    $pslfile,
    $noload,
    $outdir,
    $species, 
    $load,
    $keep_best,
    $verbose,
    );

GetOptions(
	   "debug=s"        => \$debug,
	   "test"           => \$test,
	   "species:s"      => \$species,
	   'store=s'        => \$store,
	   "database:s"     => \$db,
           "submit"         => \$submit,
	   "target:s"       => \$target,
	   "query:s"        => \$query,
           "acefile:s"      => \$acefile,
           "pslfile:s"      => \$pslfile,
           "load"           => \$load,
           "keepbest=s"     => \$keep_best,
           "verbose"        => \$verbose,
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
$outdir = $wormbase->blat;
$pslfile = "$outdir/RNAi_homols.psl" if not $pslfile;

my $coords = Coords_converter->invoke($db, 0, $wormbase);

if (not $query) {
  ###### Fetch the queries (if not supplied)
  my (@rnai, @query, @rnai_chunks, @out_ace_files, @out_lsf_files, @out_psl_files, %best_hit_only, %seen);
  
  my $tace = $wormbase->tace;
  my $tb_file = "$db/wquery/RNAi_dna.def";
  if (not -e $tb_file) {
    $log->log_and_die("Could not find Table-maker definition file $tb_file\n");
  }
  my $tb_cmd = "Table-maker -p \"$tb_file\"\nquit\n";
  open(my $tace_fh, "echo '$tb_cmd' | $tace $db |");
  while(<$tace_fh>) {
    if (/^\"(\S+)\"\s+\"(\S+)\"\s+\"(\S+)\"/) {
      my ($rnai_id, $dna_text, $clone_id) = ($1, $2, $3);
      $seen{$rnai_id}++;
      my $unique_rnai_id = sprintf("%s.DNA_text_%d", $rnai_id,  $seen{$rnai_id});
      
      if ($clone_id =~ /^mv_/ or
          $clone_id =~ /^yk\d+/) {
        $best_hit_only{$unique_rnai_id} = 1;
      }
      push @rnai, [$unique_rnai_id, $dna_text];
    }
  }
  close($tace_fh);

  my $chunk_size = 2000;      
  my $total_rnai = scalar(@rnai);
  my $chunk_count = int(scalar(@rnai) / $chunk_size);
  
  # randomise order for better distribution
  while(my $rnai = shift @rnai) {
    my $chunk_id = int(rand($chunk_count));
    
    push @{$rnai_chunks[$chunk_id]}, $rnai;
  }
  
  $log->write_to(sprintf("Fetched %d RNAi objects in %d chunks\n", $total_rnai, scalar(@rnai_chunks)));
  
  for(my $i=0; $i < @rnai_chunks; $i++) {
    my $file = $wormbase->blat . "/rnai_query.$i.fa";
    open my $fh, ">$file" or $log->log_and_die("Could not open $file for writing\n");
    foreach my $rnai (@{$rnai_chunks[$i]}) {
      printf $fh ">%s\n%s\n", @$rnai;
    }
    close($fh);
    push @query, $file;
  }
  
  if (not defined $keep_best) {
    $keep_best =  $wormbase->blat . "/rnai.best_hit_only.txt";
  }
  open my $pcf, ">$keep_best" 
      or $log->log_and_die("Could not open $keep_best for writing\n");
  foreach my $pc (keys %best_hit_only) {
    print $pcf "$pc\n";
  }
  close($pcf);

  my $lsf = LSF::JobManager->new();
  
  for(my $i=0; $i < @query; $i++) {
    my $query = $query[$i];
    my $jname = "worm_rna2genome.$i.$$";
    my $out_file = $wormbase->blat . "/rnai_out.$i.ace";
    my $psl = $wormbase->blat . "/rnai_out.$i.psl";
    my $lsf_out = $wormbase->blat . "/rnai_out.$i.lsf_out";

    my $command = ($store) 
        ? $wormbase->build_cmd_line("RNAi2Genome.pl -query $query -keepbest $keep_best -acefile $out_file -pslfile $psl", $store)
        : $wormbase->build_cmd("RNAi2Genome.pl -query $query -keepbest $keep_best -acefile $out_file -pslfile $psl");
    
    my @bsub_options = (-J => $jname, 
                        -o => $lsf_out,
                        -E => 'test -w ' . $wormbase->blat,
                        -M => 2600000,
                        -R => 'select[mem>=2600] rusage[mem=2600]'
                        );
    $lsf->submit(@bsub_options, $command);

    push @out_ace_files, $out_file;
    push @out_psl_files, $psl;
    push @out_lsf_files, $lsf_out;
  }     
       
  $lsf->wait_all_children( history => 1 );
  
  $log->write_to("All RNAi2Genome batch jobs have completed!\n");
  for my $job ( $lsf->jobs ) {    
    $log->error("$job exited non zero\n") if $job->history->exit_status != 0;
  }
  $lsf->clear;
  
  my %parent_seqs;
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
    unlink @out_ace_files, @out_psl_files, @out_lsf_files, @query;
  }

  if ($load) {
    $wormbase->load_to_database($db, $acefile, "RNAi_mapping", $log);
  }
} else {

  my (%keep_best_clones, %qlengths, $primary, $secondary, %results, %results_by_parent);
  
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
      
      my ($m_tid, $m_tstart, $m_tend) = $coords->LocateSpan($tid, $tstart, $tend);
      my $map_down = $coords->isa_clone($m_tid);
      my ($percent) = $result->get_tag_values('percent_id');
      
      my @blocks;
      foreach my $block (@b) {
        my $tstart = $block->feature1->start;
        my $tend   = $block->feature1->end;
        my $qstart = $block->feature2->start;
        my $qend   = $block->feature2->end;
        
        if ($map_down) {
          ($m_tid, $tstart, $tend) = $coords->LocateSpan($tid,$tstart,$tend);
        }
        
        if ($strand < 0) {
          ($tstart, $tend) = ($tend, $tstart);
        }
        push @blocks, [$tstart, $tend, $qstart, $qend];
      }
      
      my $parent = ($map_down) ? $m_tid : $tid;
      # strip off the suffix added earlier to make the ids unique
      $qid =~ s/\.DNA_text_\d+$//;

      push @{$results_by_parent{$parent}->{$qid}->{$method}}, [$percent, \@blocks];
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
}


$log->mail;
exit;



################################################################
#
################################################################

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
