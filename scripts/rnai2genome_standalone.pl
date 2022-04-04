#!/usr/bin/env perl
#
# rnai2genome_standalone.pl
#
# Take an fasta file with RNAi sequences, and use BLAT to perform in silico RNAi "experiment"
#


use strict;
use warnings;
use Getopt::Long;

use lib $ENV{CVS_DIR};

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

my (
  $workdir,
    $query,
    $target, 
    $pslfile,
    $outfile,
    $gffout,
    $keep_best,
    $verbose,
    $only_unmapped,
    );

GetOptions(
  "workdir=s"      => \$workdir,
  "target:s"       => \$target,
  "query:s"        => \$query,
  "outfile:s"      => \$outfile,
  "pslfile:s"      => \$pslfile,
  "verbose"        => \$verbose,
  "blatexe"        => \$BLAT,
  "keepbest=s"     => \$keep_best,
    );


############################
# recreate configuration   #
############################

die "Usage: rnai2genome_standalone.pl -query <query fasta> -target <genome fasta>\n"
    if not defined $query or not defined $target;


$workdir = "." if not defined $workdir;
$pslfile = "/tmp/rnai2genome.$$.psl" if not $pslfile;

my (%keep_best_clones, %qlengths, $primary, $secondary, %results, %results_by_parent, %has_hit);

if ($keep_best) {
  open my $pf, $keep_best or die("Could not open $keep_best\n");
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

system("$command") and die("Command '$command' failed\n");

my $out_fh;
if (defined $outfile) {
  open($out_fh, ">$outfile") or die("Could not open $outfile for writing\n");
} else {
  $out_fh = \*STDOUT;
}

open(my $fh, "<$pslfile") or die("cant open blat output $pslfile : $!\n");
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
    
    my (@blocks);
    
    foreach my $block (@b) {
      my $seg_tstart = $block->feature1->start;
      my $seg_tend   = $block->feature1->end;
      my $seg_qstart = $block->feature2->start;
      my $seg_qend   = $block->feature2->end;
      
      if ($strand < 0) {
        ($seg_tstart, $seg_tend) = ($seg_tend, $seg_tstart);
      }
      push @blocks, [$seg_tstart, $seg_tend, $seg_qstart, $seg_qend];
    }
    
    $qid =~ s/\.DNA_text_\d+$//;
    $qid =~ s/\.PCR_product_\d+$//;
    
    push @{$results_by_parent{$tid}->{$qid}->{$method}}, [$percent, $strand, \@blocks];
  }
}

foreach my $parent (keys %results_by_parent) {
  my $t_res = $results_by_parent{$parent};
  
  foreach my $query (keys %$t_res) {
    my $t_q_res = $t_res->{$query};
    
    foreach my $meth (keys %$t_q_res) {
      my @hits = @{$t_q_res->{$meth}};
      
      foreach my $hit (@hits) {
        my ($score, $strand, $blocks) = @$hit;
        foreach my $block (@$blocks) {
          printf($out_fh "%s\t%s\tRNAi_reagent\t%d\t%d\t\.\t%s\t%s\tTarget=%s %d %d +\n",
                 $parent,
                 $meth,
                 ($block->[0] < $block->[1]) ? $block->[0] : $block->[1],
                 ($block->[0] < $block->[1]) ? $block->[1] : $block->[0],
                 $score,
                 ($strand < 0) ? "-" : "+",
                 $query,
                 $block->[2], 
                 $block->[3]);
        }
      }  
    }
  }
}


exit(0);



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
