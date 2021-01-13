#!/usr/bin/env perl
#
# find_upstream_seqs.pl
#
# Creates fasta file of upstream regions of genes (potential promoters)
##
# Last updated by: $Author: klh $                      
# Last updated on: $Date: 2012-06-20 08:43:54 $        

use strict;                                      
use Getopt::Long;
use Carp;
use Storable;

use Bio::SeqIO;
use Bio::PrimarySeq;

use lib $ENV{CVS_DIR};
use Wormbase;
use Log_files;


my ($help, $debug, $test, $verbose, $store, $wormbase, $species);
my ($output, %seqs, %genes); 

my $upstream_dist = 2500;

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "test"        => \$test,
	    "verbose"     => \$verbose,
	    "store:s"     => \$store,
	    "output=s"    => \$output,
            "upstream=s"  => \$upstream_dist,
	    "species:s"   => \$species,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species,
			     );
}


my $genome_seq_file = $wormbase->genome_seq;
my $seqin = Bio::SeqIO->new(-format => 'fasta',
                             -file   => $wormbase->genome_seq);
while( my $seq = $seqin->next_seq ) {
  $seqs{ $seq->id } = $seq;
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$output = $wormbase->sequences."/upstream_sequences.dna" unless $output;
my $seqout = Bio::SeqIO->new(-format => 'fasta',
                             -file   => ">$output");

my @gff_files;
if ($wormbase->assembly_type eq 'contig') {
  push @gff_files, sprintf("%s/gene.gff", $wormbase->gff_splits);
} else {
  foreach my $chr ($wormbase->get_chromosome_names(-mito => 1, -prefix => 1)) {
    push @gff_files, sprintf("%s/%s_gene.gff", $wormbase->gff_splits, $chr);
  }
}
my $seq_count = 0;
foreach my $gff_file (@gff_files) {
  open (my $fh, $gff_file) or $log->log_and_die("Could not open $gff_file for reading\n");
  while(<$fh>) {
    next if /^\#/;
    my @l = split(/\t/, $_);
    next if $l[2] ne "gene";

    my ($id) = $l[8] =~ /Gene \"(\S+)\"/; 

    push @{$genes{$l[0]}{$l[6]}}, {
      chr => $l[0],
      start => ($l[6] eq '-') ? $seqs{$l[0]}->length - $l[4] + 1 : $l[3],
      end   => ($l[6] eq '-') ? $seqs{$l[0]}->length - $l[3] + 1 : $l[4],
      strand => $l[6],
      id => $id,
    };
    $seq_count++;
  }
}

$log->write_to("Read $seq_count genes from GFF\n");

$seq_count = 0;
foreach my $chr (sort keys %genes) {
  my $chrseq = $seqs{ $chr} ; 

  foreach my $strand (keys %{$genes{$chr}}) {
    $chrseq = $chrseq->revcom if $strand eq '-';
    my $actual_seq = $chrseq->seq;

    my @genes = sort { $a->{start} <=> $b->{start} } @{$genes{$chr}->{$strand}};

    for(my $i=0; $i<@genes; $i++) {
      my $reg_end = $genes[$i]->{start} - 1;
      next if $reg_end < 1;

      my $reg_start = $reg_end - $upstream_dist + 1;
      $reg_start = 1 if $reg_start < 1;

      my ($subseq) = substr($actual_seq, $reg_start - 1, $reg_end - $reg_start + 1);

      my $disp_reg_start = $reg_start;
      my $disp_reg_end = $reg_end;
      if ($strand eq '-') {
        $disp_reg_start =  $seqs{$chr}->length - $reg_end + 1;
        $disp_reg_end =  $seqs{$chr}->length - $reg_start + 1;
      }

      my $outseq = Bio::PrimarySeq->new( -seq => uc($subseq),
                                         -description => sprintf("%s/%d-%d", $chr, $disp_reg_start, $disp_reg_end),
                                         -id => $genes[$i]->{id});
    
      $seqout->write_seq($outseq); 
      $seq_count++;
    }

  }
}
$seqout->close;

####################################
# print some statistics to the log #
####################################

$log->write_to("\nWrote $seq_count sequences to $output\n");
$log->mail();

exit(0);



########## 
