#!/usr/bin/env perl
use strict;
use warnings;
use YAML;
use XML::Simple;
my @assemblies;
@assemblies = split "\n", do {local $/; <>} unless -t STDIN;
@assemblies = @ARGV if @ARGV;
die unless @assemblies;
for(@assemblies){
  my $doc = `curl -s "https://www.ncbi.nlm.nih.gov/assembly/$_?report=xml&format=text" `;
  print Dump XMLin(XMLin($doc));
}
