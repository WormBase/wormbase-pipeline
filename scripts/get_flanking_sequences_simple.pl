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

$flank_len = 40 if not defined $flank_len;

my $genome_fasta = shift;
my $genome_h = &read_genome_fasta($genome_fasta);

while(<>) {
  /^#/ and do {
    print;
    next;
  };
  chomp;
  my @l = split(/\t+/, $_);

  my ($chr, $gffstart, $gffend, $strand) = ($l[0], $l[3], $l[4], $l[6]);
  $chr =~ s/^CHROMOSOME_//;

  my ($lflank, $rflank);
  my ($reached_end, $reached_start);
  
  my ($lflank_end, $rflank_start) = ($gffstart - 1, $gffend + 1);
  if ($between_bases) {
    if ($gff3) {
      $lflank_end = $gffstart;
      $rflank_start = $lflank_end + 1;
    } else { 
      $lflank_end++;
      $rflank_start--;
    }
  }


  if ($rflank_start > length($genome_h->{$chr})) {
    $rflank = ">";
  }
  
  if ($lflank_end < 1) {
    $lflank = "<";
  }

  my $lflank_len = $flank_len;
  my $rflank_len = $flank_len;
  
  ITER: while(1) {
    my $lpos = $lflank_end - $lflank_len + 1;
    if ($lpos < 1) {
      $lpos = 1;
      $lflank_len = $lflank_end;
      $reached_start = 1;
    }
    my $rpos = $rflank_start;
    if ($rpos + $rflank_len - 1 > length($genome_h->{$chr})) {
      $rflank_len = length($genome_h->{$chr}) - $rpos;
      $reached_end = 1;
    }

    $lflank = substr($genome_h->{$chr}, $lpos - 1, $lflank_len);
    $rflank = substr($genome_h->{$chr}, $rpos - 1, $rflank_len);
    
    if (not $reached_start) {
      my $matches = scalar(&matches($genome_h->{$chr}, $lflank));
      if ($matches == 1) {
        $reached_start = 1;
      } else {
        $lflank_len *= 2;
        next ITER;
      }
    }
    if (not $reached_end) {
      my $matches = scalar(&matches($genome_h->{$chr}, $rflank));
      if ($matches == 1) {
        $reached_end = 1;
      } else {
        $rflank_len *= 2;
        next ITER;
      }
    }

    last;
  }

  if ($strand eq '-') {
    ($lflank, $rflank) = ($rflank, $lflank);
    
    $lflank =~ tr/ACGTN/TGCAN/;
    $rflank =~ tr/ACGTN/TGCAN/;
    
    $lflank = join("", reverse(split(//,$lflank)));
    $rflank = join("", reverse(split(//,$rflank)));
  } 
  
  if ($gff3) {
    $l[8] .= ";" if $l[8];
    $l[8] .= "leftflank=$lflank;rightflank=$rflank\n";
  } else {
    $l[8] .= " ; " if $l[8];
    $l[8] .= "LeftFlank \"$lflank\" ; RightFlank \"$rflank\"\n"
  }
  print STDERR "cannot get flanks for ($chr, $gffstart, $gffend, $strand), will use empty strings instead\n" unless ($lflank && $rflank);
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
      $current_id =~ s/^CHROMOSOME_//;
      next;
    };

    /^(\S+)/ and $genome_ref{$current_id} .= uc($1);
  }

  return \%genome_ref;
}


sub matches {
  my ($seq, $flank) = @_;

  my @matches;

  my $pos = -1;
  while (($pos = index($seq, $flank, $pos)) > -1) {
    push @matches, $pos;
    $pos++;
  }

  return @matches;
}
