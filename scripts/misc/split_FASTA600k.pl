#!/usr/local/ensembl/perl -w

my $file = shift;

die "please enter the filename to be split eg split_FASTA600k.pl slimtrembl.fasta\n" unless $file;
die "$file does not exist\n" unless (-e $file);

my $half_1 = $file."_1";
my $half_2 = $file."_2";

open (IN,"<$file") or die "$file\n";
open (FH1,">$half_1") or die "$half_1\n";
open (FH2,">$half_2") or die "$half_2\n";

my $print = *FH1;

while( <IN> ) {
  $x++ if />/;
  if( $x == 600000 ){
    $print = *FH2;
  }
  print $print $_;
}

close IN;
close FH1;
close FH2;

=pod

 Pass the name of a fasta file to be split.  The first 600,000 sequences will be in the first file the rest in another

  split_FASTA600k.pl slimtrembl.fasta

  will make slimtrembl.fasta_1 and slimtrembl.fasta_2

=cut
