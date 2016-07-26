#!/usr/bin/env perl

use strict;
use Getopt::Long;
use List::Util qw(reduce);
use JSON;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

my (
  $dbname,
  $dbhost,
  $dbuser,
  $dbport,
  $dbpass,
  $outfile, $outfh,
  $do_assembly, $do_cegma, $do_busco,
);


&GetOptions(
  'dbname=s' => \$dbname,
  'user=s'   => \$dbuser,
  'host=s'   => \$dbhost,
  'port=s'   => \$dbport,
  'pass=s'   => \$dbpass,
  'assembly' => \$do_assembly,
  'cegma'    => \$do_cegma,
  'busco'    => \$do_busco,
  'outfile=s' => \$outfile,
);


if (defined $outfile) {
  open($outfh, ">$outfile") or die "Could not open $outfile for writing\n";
} else {
  $outfh = \*STDOUT;
}

my $db = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
  '-dbname' => $dbname,
  '-host'   => $dbhost,
  '-user'   => $dbuser,
  '-port'   => $dbport,
  '-pass'   => $dbpass
    );

my $output = {};

if ($do_assembly) {

  my $break = 10; # number of consecutive Ns to break scaffolds into contigs
  my $bins = 1000; # number of bins to place sequences into

  my (@ctgs, @scafs);
  
  {
    foreach my $slice (sort { $b->length <=> $a->length } @{$db->get_SliceAdaptor->fetch_all('toplevel')}) {
      push @scafs, uc($slice->seq);
    }
  }

  bin_seqs($output, \@scafs, $bins, 1);

  while (my $seq = shift @scafs) {
    push @ctgs,  split_scaf($seq, $break);
  }


  my $ctg_output = {};
  bin_seqs($ctg_output, \@ctgs, $bins);


  $output->{contigs}               = $ctg_output->{scaffolds};
  $output->{binned_contig_counts}  = $ctg_output->{binned_scaffold_counts};
  $output->{contig_count}          = $ctg_output->{scaffold_count};
  $output->{binned_contig_lengths} = $ctg_output->{binned_scaffold_lengths};
  $output->{contig_N50}            = $ctg_output->{scaffold_N50};
  $output->{contig_L50}            = $ctg_output->{scaffold_L50};
}

if ($do_cegma) {

  my $mc = $db->get_MetaContainerAdaptor;
  
  my $cegma_complete = $mc->single_value_by_key('assembly.cegma_complete') * 1;
  my $cegma_partial = $mc->single_value_by_key('assembly.cegma_partial') * 1;

  $output->{cegma_complete} = $cegma_complete;
  $output->{cegma_partial} = $cegma_partial;

}


if ($do_busco) {

  my $mc = $db->get_MetaContainerAdaptor;
  
  my $h = {
    C => $mc->single_value_by_key('assembly.busco_complete') * 1,
    D => $mc->single_value_by_key('assembly.busco_duplicated') * 1,
    F => $mc->single_value_by_key('assembly.busco_fragmented') * 1,
    M => $mc->single_value_by_key('assembly.busco_missing') * 1,
    n => $mc->single_value_by_key('assembly.busco_number') * 1,
  };

  $output->{busco} = $h;
}  

my $json = JSON->new;
$json->pretty(1);
print $outfh $json->encode($output),"\n";

##########################################################3

sub bin_seqs {
  my ($outh, $ref, $bins, $flag) = @_;

  # sort sequences (longest first)
  my @seqs = sort { length $b <=> length $a } @$ref;
  my $count = scalar @seqs;
  my @seq_lengths = map length, @seqs;
  my $span = reduce { $a + $b } @seq_lengths;
  my $fbin = $span / $bins;
  my $bin = int ($span / $bins);

  my ($sum, $nsum, $gcsum) = (0, 0, 0);
  my ($y, $z) = (0,0);
  my (@binned_lengths,@binned_counts, @binned_gcs, @binned_ns, @bin_seqs);
  my (@AoH_ns, @AoH_gcs);
  my $catted;  

  for (my $x = 0; $x < $count; $x++){
    $sum += $seq_lengths[$x];
    $catted .= $seqs[$x];
    push @bin_seqs, $seqs[$x];
    
    while ($sum >= ($y+1)*$fbin){

      $z += $bin;
      $binned_lengths[$y] = $seq_lengths[$x];
      $binned_counts[$y] = $x+1;
      if ($flag) {
        
        if (@bin_seqs){
          my ($min_ns, $max_ns, $sum_ns);
          my ($min_gcs, $max_gcs, $sum_gcs);
          
          foreach my $str (@bin_seqs) {

            my ($ns, $gcs); 

            for(my $i=0; $i < length($str); $i++) {
              my $char = substr($str, $i, 1);
              $ns++ if $char eq 'N';
              $gcs++ if $char eq 'G' or $char eq 'C';
            }

            $ns = sprintf("%.3f", 100 * ($ns / length($str)));
            # the following tells JSON that this is a number
            $ns *= 1;
            $sum_ns += $ns;
            
            $gcs = sprintf("%.3f", 100 * ($gcs / (length($str) - $ns)));
            $gcs *= 1;
            $sum_gcs += $gcs;
            
            $min_ns  = $ns  if not defined $min_ns or $ns < $min_ns;
            $max_ns  = $ns  if not defined $max_ns or $ns > $max_ns;
            $min_gcs = $gcs if not defined $min_gcs or $gcs < $min_gcs;
            $max_gcs = $gcs if not defined $max_gcs or $gcs > $max_gcs;
          }
          
          my $mean_ns = sprintf("%.3f", $sum_ns / @bin_seqs) * 1;
          my $mean_gcs = sprintf("%.3f", $sum_gcs / @bin_seqs) * 1;

          $AoH_ns[$y] = {
            'mean' => $mean_ns,
            'max'  => $max_ns,
            'min'  => $min_ns,
          };
          $AoH_gcs[$y] = {
            'mean' => $mean_gcs,
            'max'  => $max_gcs,
            'min'  => $min_gcs,
          };
          @bin_seqs = ();
        } else {
          $AoH_ns[$y]  = $AoH_ns[$y-1];
          $AoH_gcs[$y] = $AoH_gcs[$y-1]; 
        }

        # also bin gc and n content
        my $string = substr($catted,0,$bin,'');
        # apply a correction to accommodate non-integer bin-sizes
        my $correction = ($y+1)*$fbin - $z;
        my $extra = 0;
        if ($correction >= 0.5){
          $string .= substr($catted,0,int($correction+0.5),'');
          $extra = int($correction+0.5);
          $z += $correction;
        }
        $binned_ns[$y]  = () = $string =~ /n/gi;
        $binned_gcs[$y] = () = $string =~ /[gc]/gi;
        $nsum  += $binned_ns[$y];
        $gcsum += $binned_gcs[$y];
        #print STDERR "x=$x y=$y bin=$bin extra=$extra binned_lengths[y]=$binned_lengths[$y] binned_counts[y]=$binned_counts[$y] binned_ns[y]=$binned_ns[$y] binned_gc[y]=$binned_gcs[$y]\n";
        $binned_gcs[$y] = ($bin+$extra-$binned_ns[$y] > 0) 
            ?  ($binned_gcs[$y] / ($bin+$extra-$binned_ns[$y])) * 100
            : 0;
        $binned_ns[$y] /= ($bin+$extra) / 100;
      }
      $y++;
    }
  }
  $outh->{assembly} = $span;
  $outh->{ATGC} = $span - $nsum;
  $outh->{GC} = int($gcsum / ($span - $nsum) * 10000) / 100;
  $outh->{N} = $nsum;
  #$outh->{binned_GCs} = \@binned_gcs;
  #$outh->{binned_Ns} = \@binned_ns;
  $outh->{binned_GCs} = \@AoH_gcs;
  $outh->{binned_Ns} = \@AoH_ns;
  splice @seq_lengths, 1;
  $outh->{scaffolds} = \@seq_lengths;
  $outh->{scaffold_count} = $count;
  $outh->{binned_scaffold_counts} = \@binned_counts;
  $outh->{binned_scaffold_lengths} = \@binned_lengths;
  $outh->{scaffold_N50} = $binned_lengths[ int($bins / 2) - 1];
  $outh->{scaffold_L50} = $binned_counts[ int($bins / 2) - 1];

}



sub split_scaf {
  # break scaffolds into contigs on runs of >= $break Ns
  my ($seq,$break) = @_;
  my @ctgs = split /N{$break,}/i,$seq;
  return @ctgs;
}
