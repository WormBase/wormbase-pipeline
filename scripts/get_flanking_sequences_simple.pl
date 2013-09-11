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

$flank_len = 30 if not defined $flank_len;

my $left_flank_len = $flank_len;
my $right_flank_len = $flank_len;

my $genome_fasta = shift;
my $genome_h = &read_genome_fasta($genome_fasta);
print STDERR "Read fasta\n";
while(<>) {
  /^#/ and do {
    print;
    next;
  };
  chomp;
  my @l = split(/\t+/, $_);

  my ($seq, $start, $end, $strand) = ($l[0], $l[3], $l[4], $l[6]);
  $seq =~ s/^CHROMOSOME_//;

  die "Sequence '$seq' was not found in fasta file\n" 
      if not exists $genome_h->{$seq};

  my $reference = substr($genome_h->{$seq}, $start - 1, $end - $start + 1);

  my $left_end = $start - 1;
  my $right_start = $end + 1;
  if ($between_bases) {
    $reference = "";

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
  
  my ($left_flank, $right_flank); 

  if ($right_start > length($genome_h->{$seq})) {
    $right_flank = ">";
  }
  
  if ($left_end < 1) {
    $left_flank = "<";
  }
  
  while(not defined $left_flank or not defined $right_flank) {
    my ($lcoord, $rcoord);

    $lcoord =  $left_end - $left_flank_len;
    if ($lcoord < 1) {
      $lcoord = 1;
      $left_flank_len = $left_end;
    }
    $rcoord =  $right_start - 1;
    if ($rcoord + $right_flank_len > length($genome_h->{$seq})) {
      $right_flank_len = length($genome_h->{$seq}) - $right_start + 1;
    }
    
    if (not defined $left_flank) {
      my $left = substr($genome_h->{$seq}, $lcoord, $left_flank_len);
      my $lidx = index($genome_h->{$seq}, $left);
      
      if ($lidx == $lcoord) {
        $left_flank = $left;
      } else {
        $left_flank_len += 10;
      }
    }
    
    if (not defined $right_flank) {
      my $right = substr($genome_h->{$seq}, $rcoord, $right_flank_len);
      my $ridx = rindex($genome_h->{$seq}, $right);
      
      if ($ridx == $rcoord) {
        $right_flank = $right;
      } else {
        $right_flank_len += 10;
      }
    }  
  }
  
  if ($strand eq '-') {
    ($left_flank, $right_flank) = ($right_flank, $left_flank);
    
    $left_flank =~ tr/ACGTN/TGCAN/;
    $right_flank =~ tr/ACGTN/TGCAN/;
    
    $left_flank = join("", reverse(split(//,$left_flank)));
    $right_flank = join("", reverse(split(//,$right_flank)));
  } 
  
  if ($gff3) {
    $l[8] .= ";" if $l[8];
    $l[8] .= "leftflank=$left_flank;rightflank=$right_flank;reference=$reference\n";
  } else {
    $l[8] .= " ; " if $l[8];
    $l[8] .= "LeftFlank \"$left_flank\" ; RightFlank \"$right_flank\" ; Reference \"$reference\"\n";
  }
  print STDERR "cannot get flanks for ($seq, $start, $end, $strand), will use empty strings instead\n" unless ($left_flank && $right_flank);
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
