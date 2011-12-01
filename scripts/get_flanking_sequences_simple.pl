#/usr/bin/env perl

use strict;
use Getopt::Long;

my ($flank_len,
    $between_bases,
    $gff3);

&GetOptions('flanklength=s' => \$flank_len,
            'betweenbases'  => \$between_bases,
            'gff3'          => \$gff3,
    );

$flank_len = 50 if not defined $flank_len;

my $genome_fasta = shift;

my $genome_h = &read_genome_fasta($genome_fasta);

while(<>) {
  /^#/ and do {
    print;
    next;
  };
  chomp;
  my @l = split(/\t+/, $_);

  my ($seq, $start, $end, $strand) = ($l[0], $l[3], $l[4], $l[6]);

  die "Sequence '$seq' was not found in fasta file\n" 
      if not exists $genome_h->{$seq};

  my $left_end = $start - 1;
  my $right_start = $end + 1;
  if ($between_bases) {
    if ($gff3) {
      die "For between-base features in GFF3, start should be equal to end, but $seq/$start-$end are not\n"
          if $start != $end;
      
      $left_end = $start;
      $right_start = $start + 1;
    } else {
      die "For between-base features in GFF2, start should be equal to end - 1, but $seq/$start-$end are not\n"
          if $start + 1 != $end;
      $left_end = $start;
      $right_start = $end;
    }
  }

  my $left_start = $left_end - $flank_len + 1;
  my $right_end = $right_start + $flank_len - 1;
  
  $left_start = 1 if $left_start < 0;
  $right_end = length($genome_h->{$seq}) if $right_end > length($genome_h->{$seq});

  print "$left_start $left_end $right_start $right_end\n";

  my $left_flank = substr($genome_h->{$seq}, $left_start - 1, ($left_end - $left_start + 1));
  my $right_flank = substr($genome_h->{$seq}, $right_start - 1, ($right_end - $right_start + 1));

  if ($strand eq '-') {
    ($left_flank, $right_flank) = ($right_flank, $left_flank);

    $left_flank =~ tr/ACGTN/TGCAN/;
    $right_flank =~ tr/ACGTN/TGCAN/;

    $left_flank = join("", reverse(split(//,$left_flank)));
    $right_flank = join("", reverse(split(//,$right_flank)));
  }



  if ($gff3) {
    $l[8] .= ";" if $l[8];
    $l[8] .= "Note=LeftFlank:$left_flank RightFlank:$right_flank\n";
  } else {
    $l[8] .= " ; " if $l[8];
    $l[8] .= "Note \"LeftFlank:$left_flank RightFlank:$right_flank\n";
  }

  print join("\t", @l);
}




#############################
sub read_genome_fasta {
  my ($genome_fa) = @_;

  my ($fah);
  if ($genome_fa =~ /\.gz$/) {
    open($fah, "gunzip -c $genome_fa |") or die "Could not open unzip stream to $genome_fa\n";
  } else {
    open($fah, $genome_fa) or die "Could not open $genome_fa for reading\n";
  }

  my ($current_id, %genome_ref);

  while(<$fah>) {
    /^>(\S+)/ and do {
      $current_id = $1;
      next;
    };

    /^(\S+)/ and $genome_ref{$current_id} .= uc($1);
  }

  return \%genome_ref;
}
