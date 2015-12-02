#!/usr/bin/env perl
#
# rnaseq_assign_strand.pl
#
# This script takes a BAM file containing aligned paired-end reads from an unstranded RNASeq experiment, 
# attempts to assign a strand to each pair (by looking at splice sites, and annotation), and writes a new
# SAM lines for each pair that have the FLAG field manipulated to appear like the reads came from an 
# stranded RNASeq protocol (specfically, the FR protocol)
#
# Usage: samtools view -h input.bam | perl rnaseq_assign_strand.pl -ignoreunassigned /path/to/genome.fa /path/to/wormbase/annotation.gff3 - | samtools sort -o input.strand_assigned.bam -O bam -T /tmp/bam_sort
#


use strict; 
use Getopt::Long;

my ($ignore_unassigned); 

&GetOptions("ignoreunassigned" => \$ignore_unassigned);

my $genome_file = shift @ARGV;
my $annot_file = shift @ARGV;

my $genome = &read_genome($genome_file);
my $exons  = &read_annotation($annot_file);

my (%pairs, $processed_pairs, $assigned_strand, $unassigned_strand);


while(<>) {
  if (/^\@/) {
    print;
    next;
  }

  my @l2 = split(/\t/, $_);

  # skip secondary alignments
  next if $l2[1] & 256;

  # skip reads that are not correctly oriented
  next unless $l2[1] & 2;

  my (@segs, @el);

  foreach my $el ($l2[5] =~ /\d+[NMS]/g) {
    push @el, $el
  };
  
  my $coord = $l2[3];
  foreach my $cig (@el) {
    my ($num, $let) = $cig =~ /^(\d+)(\S)/;
    
    if ($let eq 'M') {
      push @segs, [$coord, $coord + $num - 1];
    }
    
    $coord += $num;
  }

  my (%strand_votes);


  # try to match up with exons

  # remove exons with an end less than the start of the read. These will not match any read
  while(@{$exons->{$l2[2]}} and $exons->{$l2[2]}->[0]->[1] < $segs[0]->[0]) {
    #print "Discaring exon ", $l2[2],  " ", $exons->{$l2[2]}->[0]->[0], " ", $exons->{$l2[2]}->[0]->[1], "\n";
    shift @{$exons->{$l2[2]}};
  }
  my %matching;
  foreach my $seg (@segs) {
    #print "Checking @$seg\n"; 
    foreach my $ex (@{$exons->{$l2[2]}}) {
      #print "  checking againt exon @$ex\n";
      last if $ex->[0] > $seg->[1];

      if ($ex->[0] <= $seg->[1] and
          $ex->[1] >= $seg->[0]) {
        $matching{$ex} = $ex;
      }
    }
    last if keys %matching;
  }
  my @matching = values(%matching);
  foreach my $match (@matching) {
    $strand_votes{$match->[2]}++;
  }

  if (scalar(@segs) > 1) {
    my @intr;

    my $strand;

    for(my $i=1; $i < @segs; $i++) {
      my $lex = $segs[$i-1];
      my $rex = $segs[$i];

      my $int_st = $lex->[1] + 1;
      my $int_en = $rex->[0] - 1;

      #print "INTRON = $int_st $int_en\n";

      my $left = substr($genome->{$l2[2]}, $int_st - 1, 2);
      my $right = substr($genome->{$l2[2]}, $int_en - 2, 2);

      #print "$left .. $right\n";

      if ($left eq 'gt' and $right eq 'ag') {
        $strand_votes{"+"}++;
        last;
      } elsif ($left eq 'ct' and $right eq 'ac') {
        $strand_votes{"-"}++;
        last;
      }
    }
  }

  if (exists $pairs{$l2[0]}) {
    foreach my $strand (keys %{$pairs{$l2[0]}->{strand_votes}}) {
      $strand_votes{$strand} += $pairs{$l2[0]}->{strand_votes}->{strand};
    }
    my $plus = exists($strand_votes{"+"}) ? $strand_votes{"+"} : 0;
    my $minus = exists($strand_votes{"-"}) ? $strand_votes{"-"} : 0;

    my @l1 = @{$pairs{$l2[0]}->{line}};

    my ($l1_bits_on, $l1_bits_off, $l2_bits_on, $l2_bits_off);

    if ($plus != $minus) {
      if ($plus > $minus) {
        # @l1 should be read 1, and map to same strand as gene, i.e. forward strand
        # @l2 should be read 2, and map to opposite strand of gene, i.e. reverse strand
        
        # required: for @l1, switch on 0x020,0x040 (0x060) and switch off 0x010,Ox080 (0x090)
        #           for @l2, switch off 0x020,0x040 (0x060) and switch on 0x010,Ox080 (0x090)

        $l1_bits_on  = 0x060;
        $l1_bits_off = 0x090;
        $l2_bits_on  = 0x090;
        $l2_bits_off = 0x060;
      } else {
        # @l1 should be read 2, which should map to opposite strand to gene, i.e. forward strand
        # @l2 should be read 1, which should map to same strand as gene, i.e. reverse strand

        # required: for @l2, switch on 0x010,0x020 (0x030) and switch off 0x40,Ox80 (0x0C0)
        #           for @l1, switch off 0x10,0x20 (0x030) and switch on 0x40,Ox80 (0x0C0)

        $l1_bits_on  = 0x0C0;
        $l1_bits_off = 0x030;
        $l2_bits_on  = 0x030;
        $l2_bits_off = 0x0C0;
      }
        
      $l1[1] |= $l1_bits_on;
      $l2[1] |= $l2_bits_on;

      $l1[1] &= ~$l1_bits_off;
      $l2[1] &= ~$l2_bits_off;

      $assigned_strand++;

      print join("\t", @l1);
      print join("\t", @l2);
      
    } else {
      $unassigned_strand++;
      if (not $ignore_unassigned) {
        print join("\t", @l1);
        print join("\t", @l2);
      }
    }
    
    delete $pairs{$l2[0]};
    $processed_pairs++;

    if ($processed_pairs % 100000 == 0) {
      print STDERR "Done $processed_pairs pairs...\n";
    }
    
  } else {
    $pairs{$l2[0]} = {
      line => \@l2,
      strand_votes => \%strand_votes,
    };
  }

}

print STDERR "Assigned strand to $assigned_strand pairs\n";
print STDERR "Could not assigne strand to $unassigned_strand pairs\n";

##################################################

sub read_genome {
  my ($genome_file) = @_;

  my $open_cmd = ($genome_file =~ /\.gz$/) ? "gunzip -c $genome_file |" : "$genome_file";

  my (%genome, $id);

  open(my $fh, "$open_cmd");
  while(<$fh>) {
    /^\>(\S+)/ and do {
      $id = $1;
      $id =~ s/CHROMOSOME_//; 
      next;
    };
    
    /^(\S+)$/ and do {
      $genome{$id} .= $1;
    }
  }

  return \%genome;
}




sub read_annotation {
  my ($annot_file) = @_;

  my $open_cmd = ($genome_file =~ /\.gz$/) ? "gunzip -c $annot_file |" : "$annot_file";
  
  my (%cds_by_strand, %cds);
  open(my $fh, "$open_cmd");
  while(<$fh>) {
    next if /^\#/; 

    my @l = split(/\t/, $_);
    next if $l[2] ne 'CDS';

    push @{$cds_by_strand{$l[6]}->{$l[0]}}, [ $l[3], $l[4], $l[6] ];
  }

  foreach my $strand (keys %cds_by_strand) {
    foreach my $seqid (keys %{$cds_by_strand{$strand}}) {
      my @nr_segs;
      foreach my $ex (sort { $a->[0] <=> $b->[0] } @{$cds_by_strand{$strand}->{$seqid}}) {
        if (not @nr_segs or $nr_segs[-1]->[1] < $ex->[0] - 1) {
          push @nr_segs, $ex;
        } elsif ($ex->[1] > $nr_segs[-1]->[1]) {
          $nr_segs[-1]->[1] = $ex->[1];
        }
      }
      push @{$cds{$seqid}}, @nr_segs;
    }
  }

  foreach my $seqid (keys %cds) {
    $cds{$seqid} = [ sort { $a->[0] <=> $b->[0] } @{$cds{$seqid}} ];
  }

  return \%cds;
}
