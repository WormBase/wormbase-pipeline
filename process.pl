#!/usr/bin/perl

use strict;
use warnings;

my @species_list = ('echinococcus_multilocularis_prjeb122', 'echinococcus_granulosus_prjeb121', 'hymenolepis_microstoma_prjeb124');
my $counter = 0;

foreach my $species (@species_list) {
 
  warn "Species: $species";
 
  my $file = "../in/$species.txt";
  (my $out = $file) =~ s/in/out/;
  open(INFILE, $file);
  open(OUTFILE, ">$out");
  
  print OUTFILE "[ENSEMBL_INTERNAL_BIGWIG_SOURCES]\n";
  
  my $groups;
  my $files;
  
  foreach(<INFILE>) {
    chomp;
    my @parts = split("\t", $_);
    next unless $parts[0];
  
    # Remove all the Excel crap
    foreach(@parts) {
      $_ =~ s/^"//g;
      $_ =~ s/"$//g;
    }
  
    # Process colour
    $parts[4] =~ s/ //g;
    my @rgb = split(",", $parts[4]);
    my $hex = '#';
    foreach my $col (@rgb) {
      $hex .= sprintf("%02X", $col);
    }
  
    # Create the description
    my $desc = sprintf(
                '<br />ENA Sample ID: <a href="http://www.ebi.ac.uk/ena/data/view/%s">%s</a><br />
                 ENA Project ID: <a href="http://www.ebi.ac.uk/ena/data/view/%s">%s</a><br />
                 ArrayExpress ID: <a href="http://www.ebi.ac.uk/arrayexpress/experiments/%s">%s</a>',
              $parts[5], $parts[5], $parts[6], $parts[6], $parts[7], $parts[7]);
    if($parts[8]) {
      $desc .= sprintf('<br />PubMed ID: <a href="http://europepmc.org/abstract/MED/%s">%s</a>', $parts[8], $parts[8]);
    }
    $desc =~ s/\n//g;
  
    # Generate the URL
    my $url = "http://ngs.sanger.ac.uk/production/parasites/wormbase/RNASeq_alignments/myHub_$species/$species/$parts[5].bw";

    # Create the unique track ID
    my $track_id = sprintf("%03d", $counter) . "_" . $parts[5];
  
    # Create the config file text
    $groups .= sprintf("%s=%s\n", $track_id, $parts[6]);
    $files .= sprintf("[%s]\nsource_name=%s\ncaption=%s\ndescription=%s\nsource_url=%s\nsource_type=rnaseq\ndisplay=tiling\ncolour=%s\n\n",
      $track_id, $parts[1], $parts[0], $desc, $url, $hex);
  
    $counter++;

  }
  
  print OUTFILE $groups, "\n";
  print OUTFILE $files;
  
  close(INFILE);
  close(OUTFILE);
 
}
 
