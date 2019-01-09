#!/usr/bin/env perl
use strict;
use warnings;
use YAML;
use XML::Simple;
while(<>){
  chomp;
  my $doc = `curl -s "https://www.ncbi.nlm.nih.gov/assembly/$_?report=xml&format=text" `;
  print Dump XMLin(XMLin($doc));
}
