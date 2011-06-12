#!/software/bin/perl -w

# Version: $Version: $
# Last updated by: $Author: klh $
# Last updated on: $Date: 2011-06-12 09:09:33 $

use strict;
use warnings;
use Getopt::Long;

use lib $ENV{'CVS_DIR'};
use lib '/software/worm/ensembl/bioperl-live';
use lib '/software/worm/lib';

use Wormbase;
use Log_files;

use LSF::JobManager;
use Ace;
use Coords_converter;
use Bio::Tools::Blat;

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
    $problem_file,
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
           "problemfile=s"  => \$problem_file,
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
                             -organsim => $species); 
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
  my (@rnai, @query, @rnai_chunks, @out_ace_files, @out_pslfiles, %bad_clones);
  
  my $ace = Ace->connect('-path' => $db, 
                         '-program' => $wormbase->tace ) or $log->log_and_die(Ace->error."\n");
  
  my $RNAis = $ace->fetch_many('-class' => "RNAi");
  while(my $rnai = $RNAis->next) {
    if ($rnai->DNA_text) {
      # we need to register yk clones and mv PCR products for a hack later
      if ($rnai->PCR_product and $rnai->PCR_product->name =~ /^mv_/ or
          $rnai->Sequence and $rnai->Sequence->name =~ /^yk\d+/) {
        $bad_clones{$rnai->name} = 1;
      }
      push @rnai, $rnai;
    }
  }
  
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
      printf $fh ">%s\n%s\n", $rnai->name, $rnai->DNA_text;
    }
    close($fh);
    push @query, $file;
  }
  
  my $prob_file = $wormbase->blat . "/rnai.problem_clones.txt";
  open my $pcf, ">$prob_file" 
      or $log->log_and_die("Could not open $prob_file for writing\n");
  foreach my $pc (keys %bad_clones) {
    print $pcf "$pc\n";
  }
  close($pcf);

  my $lsf = LSF::JobManager->new();
  
  for(my $i=0; $i < @query; $i++) {
    my $query = $query[$i];
    my $jname = "worm_rna2genome.$i.$$";
    my $out_file = $wormbase->blat . "/rnai_out.$i.ace";
    my $psl = $wormbase->blat . "/rnai_out.$i.psl";

    my $command = ($store) 
        ? $wormbase->build_cmd_line("RNAi2Genome.pl -query $query -problem $prob_file -acefile $out_file -pslfile $psl", $store)
        : $wormbase->build_cmd("RNAi2Genome.pl -query $query -problem $prob_file -acefile $out_file -pslfile $psl");
    
    my @bsub_options = (-J => $jname, 
                        #-o => "/nfs/users/nfs_k/klh/test." . $i . ".out",
                        #-e => "/nfs/users/nfs_k/klh/test." . $i . ".err",
                        -o => '/dev/null',
                        -e => '/dev/null',
                        -E => 'test -w ' . $wormbase->blat,
                        );
    $lsf->submit(@bsub_options, $command);

    push @out_ace_files, $out_file;
    push @out_pslfiles, $psl;
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
      
  #unlink @out_ace_file, @query, @out_pslfiles;

  if ($load) {
    $wormbase->load_to_database($db, $acefile, "RNAi_mapping", $log);
  }
} else {

  my (%problem_clones, $primary, $secondary, %results, %results_by_parent);
  
  if ($problem_file) {
    open my $pf, $problem_file or $log->log_and_die("Could not open $problem_file\n");
    while(<$pf>) {
      /^(\S+)/ and $problem_clones{$1} = 1;
    }
  }
  
  my $command = "$BLAT $target $query $pslfile $BLAT_OPTIONS";
  
  print STDERR "$command\n";
  system("$command") and $log->error("Command '$command' failed\n");
  
  open my $out_ace, ">$acefile";

  $log->write_to("Processing PSL file $pslfile...\n");
    
  open(my $fh, "<$pslfile") or $log->log_and_die("cant open blat output $pslfile : $!\n");
  my $parser = Bio::Tools::Blat->new('-fh' => $fh);
  
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

    &print_blat($percent, $total_aligned, $max_block_size, $result);
    
    # only need to keep the best primary hit for the problem clones
    if ($problem_clones{$result->seq_id}) {
      if ($total_aligned >= $PRIMARY_LEN and $percent >= $PRIMARY_QUAL and
          (not exists $results{$result->seq_id} or
           $total_aligned > $results{$result->seq_id}->[0]->[1])) {
        $results{$result->seq_id} = [[$percent, $total_aligned, $max_block_size, $result]];
      }
    } else {
      if ($total_aligned >= $PRIMARY_LEN and $percent >= $PRIMARY_QUAL or
          $total_aligned >= $SECONDARY_LEN and $percent > $SECONDARY_QUAL) {
        push @{$results{$result->seq_id}}, [$percent, $total_aligned, $max_block_size, $result];
      }
    }
  }
  
  ($primary, $secondary) = &filter_hits(\%results, $CHECK_BLOCK_SIZE);
  
  foreach my $type_pair ([$primary, 'RNAi_primary'], [$secondary, 'RNAi_secondary']) {
    my ($list, $method) = @$type_pair;
    
    foreach my $result (@$list) {
      
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
  my ($res, $check_block_size) = @_;

  my (@all_primary, @all_secondary);

  foreach my $qid (keys %$res) {

    my (@primary, @secondary);

    foreach my $quad (@{$res->{$qid}}) {
      my ($percent, $total_aligned, $max_block_size, $hit) = @$quad;

      if ($total_aligned >= $PRIMARY_LEN and $percent >= $PRIMARY_QUAL) {
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
      @primary = grep { $_->[2] >= $PRIMARY_LEN } @primary;
      @secondary = grep { $_->[2] >= $SECONDARY_LEN } @secondary;
    }

    #if (exists $problem_clones{$qid}) {
    #  push @all_primary, $primary[0]->[3]; 
    #} else {
    push @all_primary, map { $_->[3] } @primary;
    push @all_secondary, map { $_->[3] } @secondary;
    #}
  }

  return (\@all_primary, \@all_secondary);
}


sub print_blat {
  my ($perc, $aligned, $max_bsize, $res) = @_;

  my @b = $res->get_SeqFeatures;
  printf "RESULT: %s %d-%d %s %s %d %d\n", $b[0]->feature1->seq_id, $res->start, $res->end, $res->seq_id, $perc, $aligned, $max_bsize;;
  foreach my $b (@b) {
    printf " SEG: %d-%d  %d-%d (%d)\n", $b->feature1->start, $b->feature1->end, $b->feature2->start, $b->feature2->end, ($b->feature1->end - $b->feature1->start + 1);
  }
}
