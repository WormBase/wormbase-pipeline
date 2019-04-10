#!/usr/bin/env perl
use strict;
use warnings;
use feature 'say';

for my $taxon (@ARGV){
  say for get_assemblies_for_taxon($taxon);
}

sub get_assemblies_for_taxon {
  my ($taxon) = @_;
  my @result;
  local $/ = "--";
  open(my $fh, "curl -s 'https://www.ncbi.nlm.nih.gov/assembly/organism/browse-organism/$taxon/all/?p\$l=Ajax&page=1&pageSize=1000&sortCol=4&sortDir=-1' | grep -B 2 href |");
  while(<$fh>){
    my ($accession) = /<a href="https:\/\/www.ncbi.nlm.nih.gov\/assembly\/(.*?)\/">.*?</;
    push @result, $accession;
  }
  return @result;
}
