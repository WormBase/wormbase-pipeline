#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use LWP::UserAgent;
use JSON;

my @species_list = ('echinococcus_granulosus_prjeb121', 'echinococcus_multilocularis_prjeb122', 'hymenolepis_microstoma_prjeb124', 'strongyloides_stercoralis_prjeb528');
my $counter = 0;

# Create the hub.txt file
mkdir 'myHub' unless -d 'myHub';
open(OUTFILE, '>myHub/hub.txt');
print OUTFILE "hub WBPS-RNASeq\nshortLabel RNA-Seq Alignments\nlongLabel RNA-Seq Alignments\ngenomesFile genomes.txt\nemail parasite-help\@sanger.ac.uk\n";
close(OUTFILE);

open(OUTFILE, '>myHub/genomes.txt');
close(OUTFILE);

foreach my $species (@species_list) {
 
  warn "Species: $species";
  $counter = 0;

  # Write to genomes.txt
  open(OUTFILE, '>>myHub/genomes.txt');
  ## Use the REST API to lookup the assembly name
  my $url = "http://parasite.wormbase.org/api/info/assembly/$species?content-type=application/json";
  my $ua = LWP::UserAgent->new();
  my $response = $ua->get($url);
  if ($response->is_success) {
    my $output = from_json($response->decoded_content);
    my $assembly = $output->{'assembly_name'};
    print OUTFILE "genome $assembly\ntrackDb $species/trackDb.txt\n\n";
  }
  close(OUTFILE);
 
  my $file = "./in/$species.txt";
  open(INFILE, $file);
  mkdir "myHub/$species" unless -d "myHub/$species";
  open(OUTFILE, ">myHub/$species/trackDb.txt");

  my $groups;
  my $files;
  my %names;
  
  foreach(<INFILE>) {
    chomp;
    my @parts = split("\t", $_);
    next unless $parts[0];
  
    # Remove all the Excel crap
    foreach(@parts) {
      chomp;
      $_ =~ s/^"//g;
      $_ =~ s/"$//g;
    }
 
    # Create the unique track ID
    my $track_id = sprintf("%03d", $counter) . "_" . $parts[5];
 
    # Process colour
    $parts[4] =~ s/ //g;
    my @rgb = split(",", $parts[4]);
    my $hex = '#';
    foreach my $col (@rgb) {
      $hex .= sprintf("%02X", $col);
    }
  
    # Create the project description
    my $desc = sprintf('
                 ENA Project ID: <a href="http://www.ebi.ac.uk/ena/data/view/%s">%s</a><br />
                 ArrayExpress ID: <a href="http://www.ebi.ac.uk/arrayexpress/experiments/%s">%s</a>',
              $parts[6], $parts[6], $parts[7], $parts[7]);
    if($parts[8]) {
      $desc .= sprintf('<br />PubMed ID: <a href="http://europepmc.org/abstract/MED/%s">%s</a>', $parts[8], $parts[8]);
    }
    $desc =~ s/\n//g;
    mkdir "myHub/$species/doc" unless -d "myHub/$species/doc";
    open(HTMLOUT, ">myHub/$species/doc/$parts[6].html");
    print HTMLOUT $desc;
    close(HTMLOUT);

    # Create the sample description
    $desc = sprintf(
                'ENA Sample ID: <a href="http://www.ebi.ac.uk/ena/data/view/%s">%s</a><br />
                 ENA Project ID: <a href="http://www.ebi.ac.uk/ena/data/view/%s">%s</a><br />
                 ArrayExpress ID: <a href="http://www.ebi.ac.uk/arrayexpress/experiments/%s">%s</a>',
              $parts[5], $parts[5], $parts[6], $parts[6], $parts[7], $parts[7]);
    if($parts[8]) {
      $desc .= sprintf('<br />PubMed ID: <a href="http://europepmc.org/abstract/MED/%s">%s</a>', $parts[8], $parts[8]);
    }
    $desc =~ s/\n//g;
    mkdir "myHub/$species/doc" unless -d "myHub/$species/doc";
    open(HTMLOUT, ">myHub/$species/doc/$track_id.html");
    print HTMLOUT $desc;
    close(HTMLOUT);
 
    # Generate the URL
    my $url = "http://ngs.sanger.ac.uk/production/parasites/wormbase/RNASeq_alignments/myHub_$species/$species/$parts[5].bw";

    # Create the trackDB text
    $groups .= sprintf("track %s\nsuperTrack on\ngroup %s\nshortLabel %s\nlongLabel %s\nhtml doc/%s\n\n", $parts[6], $parts[6], $parts[6], $parts[6], $parts[6]) unless $names{$parts[6]};
    $files .= sprintf("track %s\nparent %s\nbigDataUrl %s\nshortLabel %s\nlongLabel %s\ncolor %s\nhtml doc/%s\n\n", $track_id, $parts[6], $url, $parts[0], $parts[1], $parts[4], $track_id);
    $counter++;
    $names{$parts[6]} = 1;

  }
  
  print OUTFILE $groups, "\n";
  print OUTFILE $files;
  
  close(INFILE);
  close(OUTFILE);
 
}
 
