#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Ace;

my (@sources, $target, $database, $file, @merges);


&GetOptions(
  "source=s@" => \@sources,
  "target=s"  => \$target,
  "file=s"    => \$file,
  "database=s" => \$database,
    );

die "You must supply a valid database path to (e.g.) geneace\n"
    if not $database or not -d $database;

if ($file) {
  open(my $fh, $file) or die "Could not open $file for reading\n";
  while(<$fh>) {
    my ($tgt, @src) = split(/\s+/, $_);

    push @merges, [$tgt, @src];
  }
} elsif (@sources and $target) {
  push @merges, [$target, @sources];
} else {
  die "You must supply a target and at least one source, or a file of such\n";
}

my @lt = localtime(time);
my $date = sprintf("%4d%2d%2d", $lt[5] + 1900, $lt[4] +1, $lt[3]); 

my $db = Ace->connect(-path => $database ) or die "Could not connect to Acedb database $database\n"; 

foreach my $merge (@merges) {
  my ($tgt, @src) = @$merge;
  
  my @tgt_lines;

  foreach my $src (@src) {
    my $src_obj = $db->fetch( Variation => "$src" );
    if (not defined $src_obj) {
      die "Could not fetch object for $src\n";
    }

    my @strain_lines = $src_obj->at("Origin.Strain");
    my @analysis = $src_obj->at("Origin.Analysis");
    my @other_name = $src_obj->at("Name.Other_name");
    my @refs = $src_obj->at("Reference");

    foreach my $nm (@other_name) {
      push @tgt_lines, sprintf("Other_name %s", $nm);
    }
    foreach my $sl (@strain_lines) {
      my $txt = sprintf("Strain %s %s %s", $sl->name, $sl->right, $sl->right->right);
      push @tgt_lines, $txt;
    }
    foreach my $ana (@analysis) {
      push @tgt_lines, sprintf("Analysis %s", $ana);
    }
    push @tgt_lines, "Other_name $src";
    foreach my $ref (@refs) {
      push @tgt_lines, sprintf("Reference %s", $ref);
    }

    &kill_var($src, $tgt);
  }

  print "\nVariation : \"$tgt\"\n";
  foreach my $ln (@tgt_lines) {
    print "$ln\n";
  }
}



sub kill_var {
  my ($var, $merged_into) = @_;

  print "\nVariation : \"$var\"\n";
  print "-D Evidence\n";
  print "-D Public_name\n";
  print "-D Sequence_details\n";
  print "-D Variation_type\n";
  print "-D Strain\n";
  print "-D Species\n";
  print "-D Laboratory\n";
  print "-D Author\n";
  print "-D Person\n";
  print "-D Affects\n";
  print "-D Reference\n";
  print "-D Analysis\n";

  print "\nVariation : \"$var\"\n";
  print "Dead Remark \"[$date] Merged into $merged_into\"\n";
  print "Merged_into $merged_into\n\n";

}
