#! /usr/bin/env perl 
# change mRNAs to nontranslating_transcript for transcripts of genes that fail the ProteinTranslation datacheck
# remove CDSs these had as a parent
# this is appropriate for cases where a handful of "protein coding" transcripts in a user submitted GFF fail the ProteinTranslation datacheck
# adapted from Wojtek's change_transcripts_to_noncoding.pl

use strict;
use warnings;

@ARGV or die "Usage: $0 in.gff transcripts.txt";

my $f = shift;
open (my $fh, '<', $f) or die;

my $g = shift;
open (my $gh, '<', $g) or die;

my %ts;
while (<$gh>) {
  chomp;
  $ts{$_}++;
}
sub extract_attribute {
   my ($attr,$str) = @_;
   $str =~ /$attr=(.*?);/;
   $str =~ /$attr=(.*)$/ unless $1;
   chomp $1;
   return $1;
}
my %ids;
my @cds;
while (<$fh>) {
    print and next if /^#/;
    my @c = split "\t";
    if ($c[2] eq "CDS"){ 
	push @cds, \@c; 
	next;
    }
    if ( $c[2] eq "mRNA" ) {
        my $name = &extract_attribute("Name", $c[8]);
	my $id = &extract_attribute("ID", $c[8]);
        warn "$c[8] without a name? " unless $name;
	if ($ts{$name}) {
         $c[2] = "nontranslating_transcript";
         $ids{$id}++;
        }
    }
   print join "\t", @c;
}

foreach my $cds (@cds){
	my $parent = &extract_attribute("Parent", $cds->[8]);
	next if $ids{$parent};
	print join "\t", @$cds;
}
