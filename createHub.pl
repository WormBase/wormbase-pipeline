#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use LWP::UserAgent;
use HTML::Entities;
use XML::Simple;
use JSON;

my @species_list = `ls ./in/*.txt | xargs -n1 basename`;
foreach(@species_list) {
  chomp;
  $_ =~ s/\.txt//;
}

my $counter = 0;

# Create the hub.txt file
mkdir 'myHub' unless -d 'myHub';
open(OUTFILE, '>myHub/hub.txt');
print OUTFILE "hub WBPS-RNASeq\nshortLabel WormBase ParaSite\nlongLabel RNA-Seq Alignments for WormBase ParaSite\ngenomesFile genomes.txt\nemail parasite-help\@sanger.ac.uk\n";
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

    # Deal with study/project being wrong way round
    if($parts[6] =~ /^(ERS|SRS)/ && $parts[5] =~ /^(ERP|SRP)/) {
      my $temp = $parts[5];
      $parts[5] = $parts[6];
      $parts[6] = $temp;
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
    my $proj_desc = sprintf('
                 ENA Project ID: <a href="http://www.ebi.ac.uk/ena/data/view/%s">%s</a><br />
                 ArrayExpress ID: <a href="http://www.ebi.ac.uk/arrayexpress/experiments/%s">%s</a>',
              $parts[6], $parts[6], $parts[7], $parts[7]);
    if($parts[6]) {
      my $ref = get_ENA_project($parts[6]);
      $proj_desc .= "<br /><br />$ref<br />" if $ref;
    }
    if($parts[8]) {
      my $ref = get_reference($parts[8]);
      if($ref) {
        $proj_desc .= "<br />Reference<br />$ref"; 
      } else {
        $proj_desc .= sprintf('<br />PubMed ID: <a href="http://europepmc.org/abstract/MED/%s">%s</a>', $parts[8], $parts[8]);
      }
    }
    $proj_desc =~ s/\n//g;
    mkdir "myHub/$species/doc" unless -d "myHub/$species/doc";
    open(HTMLOUT, ">myHub/$species/doc/$parts[6].html");
    print HTMLOUT $proj_desc;
    close(HTMLOUT);

    # Generate the URL
    my $url = "http://ngs.sanger.ac.uk/production/parasites/wormbase/RNASeq_alignments/$species/$parts[5].bw";

    # Re-format the WB life stage
    $parts[10] =~ s/^WBls_/WBls:/ if $parts[10];

    # Create the sample description
    my $desc = sprintf(
                'ENA Sample ID: <a href="http://www.ebi.ac.uk/ena/data/view/%s">%s</a><br />
                 %s%s%s',
              $parts[5], $parts[5], $parts[9] ? "ChEBI: <a href=\"https://www.ebi.ac.uk/chebi/searchId.do?chebiId=$parts[9]\">$parts[9]</a><br />" : '',
              $parts[10] ? "WormBase Life Stage: <a href=\"http://www.wormbase.org/species/all/life_stage/$parts[10]\">$parts[10]</a><br />" : '',  $proj_desc);
    $desc =~ s/\n//g;
    $desc .= "<br /><br />This data comes from URL: $url";
    mkdir "myHub/$species/doc" unless -d "myHub/$species/doc";
    open(HTMLOUT, ">myHub/$species/doc/$track_id.html");
    print HTMLOUT $desc;
    close(HTMLOUT);
 
    # Create the trackDB text
    $groups .= sprintf("track %s\nsuperTrack on\ngroup %s\nshortLabel %s\nlongLabel %s\nhtml doc/%s\n\n", $parts[6], $parts[6], $parts[6], $parts[6], $parts[6]) unless $names{$parts[6]};
    $files .= sprintf("track %s\nparent %s\ntype bigWig\nbigDataUrl %s\nshortLabel %s\nlongLabel %s\ncolor %s\nhtml doc/%s\n\n", $track_id, $parts[6], $url, $parts[0], $parts[1], $parts[4], $track_id);
    $counter++;
    $names{$parts[6]} = 1;

  }
  
  print OUTFILE $groups, "\n";
  print OUTFILE $files;
  
  close(INFILE);
  close(OUTFILE);
 
}

sub get_reference {
  # Form a reference from a PubMed ID
  my ($pmid) = @_;
  my $url = "http://www.ebi.ac.uk/europepmc/webservices/rest/search/query=$pmid";
  my $ua = LWP::UserAgent->new();
  my $response = $ua->get($url);
  if ($response->is_success) {
    my $result = XMLin($response->decoded_content);
    my $text = encode_entities("$result->{resultList}->{result}->{authorString} ") . "<a href=\"http://europepmc.org/abstract/MED/$pmid\">" . encode_entities("$result->{resultList}->{result}->{title}") . "</a> <em>" . encode_entities($result->{resultList}->{result}->{journalTitle}) . "</em>" . encode_entities(", $result->{resultList}->{result}->{pubYear}");      # encode_entities will encode any symbolic characters (such as ligatures in author names) into the correct HTML
    $text .= encode_entities(";$result->{resultList}->{result}->{journalVolume}($result->{resultList}->{result}->{issue}):$result->{resultList}->{result}->{pageInfo}") if $result->{resultList}->{result}->{journalVolume} && $result->{resultList}->{result}->{issue} && $result->{resultList}->{result}->{pageInfo};
    return $text;
  }
}

sub get_ENA_project {
  my ($id) = @_;
  my $url = "http://www.ebi.ac.uk/ena/data/view/$id&display=xml";
  my $ua = LWP::UserAgent->new();
  my $response = $ua->get($url);
  if ($response->is_success) {
    my $result = XMLin($response->decoded_content);
    my $text = $result->{STUDY}->{DESCRIPTOR}->{STUDY_DESCRIPTION} ? "Project Description<br />$result->{STUDY}->{DESCRIPTOR}->{STUDY_DESCRIPTION}" : '';
    return $text;
  }
}
 
