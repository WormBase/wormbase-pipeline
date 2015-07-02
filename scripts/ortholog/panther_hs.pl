#!/bin/env perl
# this script will 
open INF, 'blah.txt';

my %u2e;

while(<INF>){
  chomp;
  my @F=split;
  $u2e{$F[1]}=$F[0];
}

my %sp2lin= (
  MOUSE => 'Mus musculus',
  HUMAN => 'Homo sapiens',
  DANRE => 'Danio rerio',
  DROME => 'Drosophila melanogaster',
  YEAST => 'Saccharomyces cerevisiae',
);

while(<>){
    next unless /WBGene/;
    my ($from)=$_=~/(WBGene\d+)/;
    my @F=split;
    my ($to) = grep{$_=~!/WBGene/}($F[0],$F[1]);
    next unless $to;

    my ($thing) = $to=~/^(\w+)\|/;
    my ($uniprot)=$to=~/UniProtKB=(\w+)/;
    
    next unless $u2e{$uniprot};

    print "Gene : \"$from\"\n";
    print 'Ortholog "',$u2e{$uniprot},"\" \"$sp2lin{$thing}\" From_analysis Panther\n\n";

    print "Gene : \"",$u2e{$uniprot},"\"\n";
    print 'Ortholog "',$from,'" "Caenorhabditis elegans" from_analysis Panther',"\n\n"
}
