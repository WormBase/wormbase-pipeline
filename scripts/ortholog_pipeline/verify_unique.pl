#!/usr/bin/perl

# Verify that the ortholog set only contains unique one-to-one mappings
#  Todd Harris (harris@cshl.org)
#  DATE : 05 May 2003

use strict;

my $counter;
my %orthos;
while (<>) {
  $counter++;
  chomp;
  my ($cb,$ce,$eval,$method) = split("\t");
  $orthos{$ce} = $cb;
}	


# Invert the hash.
# Don't want multiple elegans or briggsae

my %inverted = map { $orthos{$_} => $_ } keys %orthos;

print "Total orthos in list : $counter\n";
print "Briggsae to elegans mappings: ",scalar keys %inverted,"\n";
print "Elegans to briggsae mappings: ",scalar keys %orthos,"\n";
print "WARNING: DUPLICATES FOUND IN LIST!!\n" if (scalar keys %orthos != scalar keys %inverted);
