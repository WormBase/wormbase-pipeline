#!/usr/bin/env perl
#
# run_inverted.pl
#
# Usage : run_inverted.pl [-options] <sequence_file>

#################################################################################
# Initialise variables                                                          #
#################################################################################
 
use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Storable;
use Log_files;

use Bio::PrimarySeq;
use Bio::SeqIO;

##############################
# command-line options       #
##############################
                                                                                                                                     
my ($debug, $test, $store, $einverted, $sequence, $all, $species, $wormbase, $acefile); 

GetOptions (
  "debug:s"     => \$debug,
  "test"        => \$test,
  "store:s"     => \$store,
  "species:s"   => \$species,
  "einverted:s" => \$einverted,
  "all"         => \$all,
  "sequence:s"  => \$sequence,
  "acefile:s"   => \$acefile,
	    );

if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			     -organism => $species
			   );
}


$einverted = "einverted" if not defined $einverted;
$acefile = $wormbase->acefiles."/inverted_repeats.ace" if not defined $acefile;
my $log = Log_files->make_build_log($wormbase);

my (%clone2seq, @batches, %output);

if (defined $sequence) {
  push @batches, {
    size => length($sequence), 
    seqs => [Bio::PrimarySeq->new(-id => "anonymous", 
                                  -seq => $sequence)],
  };
  $clone2seq{"anonymous"} = $sequence;
} elsif ($all) {
  %clone2seq  = $wormbase->FetchData("clone2sequence_${\$wormbase->species}");
  foreach my $id (sort { length($clone2seq{$b}) <=> length($clone2seq{$a}) } keys %clone2seq) {
    my $seq = $clone2seq{$id};
    $seq =~ s/\s//g; 
    next if not $seq;
    
    if (not @batches or $batches[-1]->{size} > 1000000) {
      push @batches, {
        size => 0,
        seqs => [],
      };
    }
    $batches[-1]->{size} += length($seq);
    push @{$batches[-1]->{seqs}}, Bio::PrimarySeq->new(-id => $id,
                                                       -seq => $seq);
  }
} else {
  $log->log_and_die("You must supply either -sequence or -all\n");
}

# Loop through all clone which need to be dealt with

my ($clone, $score, $percent,$gaps);
my ($loop_1_start,$loop_1_stop);
my ($loop_2_start,$loop_2_stop);
my ($loop_len);
my $tag;

foreach my $batch (@batches) {
  $log->write_to(sprintf("Processing batch of %d seqs (total length %d)\n", scalar(@{$batch->{seqs}}), $batch->{size}));

  my $tmp_file = "/tmp/inverted_temp.$$.fa";
  my $seqio = Bio::SeqIO->new(-format => 'fasta',
                              -file => ">$tmp_file");
  foreach my $seq (@{$batch->{seqs}}) {
    $seqio->write_seq($seq);
  }
  $seqio->close;

  # inverted output format: 
  #
  # Score 52: 20/22 ( 90%) matches, 0 gaps
  #      3810 ataaaaactcgaattcaaaaaa 3831
  #                                 
  #      4164 tatttgtgagcttaattttttt 4143

  # Final ace output
  #
  # Feature Inverted        3810    3831    90      "loop 353"
  
  open (my $einv_fh, "$einverted $tmp_file -outfile stdout -gap 12 -threshold 50 -match 3 -mismatch -4 -outseq /dev/null -auto| ")
      or $log->log_and_die("Could not initiate einverted command\n");
  while (<$einv_fh>) {
    # Main data 
    if  (/^(\S+): Score (\d+)\: \S+ \(\s*(\d+)\%\) matches\, (\d+) gaps/) {
      ($clone, $score, $percent,$gaps) = ($1, $2, $3, $4);
      $tag = 1;
      next;
    }
    # start loop
    if (($tag == 1) && (/^\s+(\d+) \S+ (\d+)/)) {
      ($loop_1_start,$loop_1_stop) = ($1,$2);
      $tag++;
      next;
    } 
    
    # end loop
    if (($tag == 2) && (/^\s+(\d+) \S+ (\d+)/)) {
      ($loop_2_start,$loop_2_stop) = ($1,$2);
      
      $loop_len = $loop_2_start - $loop_1_start - 1;
      
      # output ace format for both stem structures at the same time
      
      if ($gaps > 1) {
        push @{$output{$clone}}, "Feature Inverted $loop_1_start $loop_1_stop $percent \"loop $loop_len, $gaps gaps\"\n";
        push @{$output{$clone}}, "Feature Inverted $loop_2_start $loop_2_stop $percent \"loop $loop_len, $gaps gaps\"\n";
      }
      elsif ($gaps == 1) {
        push @{$output{$clone}}, "Feature Inverted $loop_1_start $loop_1_stop $percent \"loop $loop_len, 1 gap\"\n";
        push @{$output{$clone}}, "Feature Inverted $loop_2_start $loop_2_stop $percent \"loop $loop_len, 1 gap\"\n";
      }
      elsif ($gaps == 0) {
        push @{$output{$clone}}, "Feature Inverted $loop_1_start $loop_1_stop $percent \"loop $loop_len\"\n";
        push @{$output{$clone}}, "Feature Inverted $loop_2_start $loop_2_stop $percent \"loop $loop_len\"\n";
      }
      
      $tag = 0;
    }
    
  }
  close($einv_fh);
  unlink $tmp_file;
}

open(my $outfh, ">$acefile") or $log->log_and_die("Failed to open output acefile: '$acefile'\n");

foreach my $clone (sort keys %output) {  
  print $outfh "\nSequence : \"$clone\"\n";
  print $outfh "Feature_data \"${clone}:inverted\" 1 ", length($clone2seq{$clone}), "\n\n";
                                                               
  print $outfh "Feature_data \"${clone}:inverted\"\n";
  foreach my $outline (@{$output{$clone}}) {
    print $outfh $outline;
  }
}

close($outfh);
$log->mail();
exit (0);


__END__
