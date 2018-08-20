#!/usr/bin/env perl
# Given an input stream and a config, save lines in target folder
# Duplicate the lines when there are multiple species per taxon
use File::Spec;
use YAML;
die "Usage: $0 <input file> species_per_taxon.yaml <target folder>" unless @ARGV == 3 ;

die "Expected a folder: $ARGV[2]" unless -d  $ARGV[2];

my $species_per_taxon = YAML::LoadFile($ARGV[1]) or die "could not read $ARGV[1]"; 
die "Seems empty: $ARGV[1]" unless %$species_per_taxon;
local $gaf_file;
if ($ARGV[0] eq "-") {
  $gaf_file = STDIN;
} else {
  open ($gaf_file, "<", $ARGV[0]) or die "could not read $ARGV[0]";
} 
my %fhs;
while(<$gaf_file>){
  my ($pre, $taxon, $post) = split /tgt_taxon=([0-9]*)/;
  for my $species (@{$species_per_taxon->{$taxon} // [] }){
    open($fhs{$species}, ">", File::Spec->catfile($ARGV[2], "annotations_ensembl-$species.gpa")) unless $fhs{$species};
    print { $fhs{$species} } "${pre}tgt_species=${species}${post}";
  }
}

for my $k (keys  %fhs) {
  close $fhs{$k};
}
