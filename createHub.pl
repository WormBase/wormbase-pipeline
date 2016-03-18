#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use LWP::UserAgent;
use HTML::Entities;
use XML::Simple;
use JSON;

my @species_list = `ls ./in/*.ini | xargs -n1 basename`;
foreach(@species_list) {
  chomp;
  $_ =~ s/\.ini//;
}

my $counter = 0;

# Create the hub.txt file
mkdir 'myHub' unless -d 'myHub';
open(OUTFILE, '>myHub/hub.txt');
print OUTFILE "hub WBPS-RNASeq\nshortLabel RNA-Seq Alignments\nlongLabel RNA-Seq Alignments for WormBase ParaSite\ngenomesFile genomes.txt\nemail parasite-help\@sanger.ac.uk\n";
close(OUTFILE);

open(OUTFILE, '>myHub/genomes.txt');
close(OUTFILE);

foreach my $in_file (@species_list) {
 
  warn "Species: $in_file";
  $counter = 0;

  my $file = "./in/$in_file.ini";
  open(INFILE, $file);

  # Get the BioProject
  my $bioproject;
  foreach(<INFILE>) {
    chomp;
    $_ =~ /^general_bioproject=(.*?)$/;
    $bioproject = $1 if $1;
    last if $bioproject;
  }
  warn "Skipping as no BioProject" && next unless $bioproject;
  my $species = $in_file . "_" . lc($bioproject);
  warn "BioProject: $bioproject";

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
 
  mkdir "myHub/$species" unless -d "myHub/$species";
  open(OUTFILE, ">myHub/$species/trackDb.txt");

  my $groups;
  my $files;
  my %names;

  # Get each study
  my %studies;
  my $current = 'general';
  seek(INFILE, 0, 0);
  foreach(<INFILE>) {
    chomp;
    next if $_ =~ /^$/;
    if($_ =~ /^\[(.*)\]$/) {
      $current = $1;
    } else {
      $studies{$current} .= $_ . "\n";
    }
  }
  
  # Get the details for each study
  foreach my $study (keys %studies) {
    next if $study eq 'general';
    warn "-- Study: $study";
    my %ini;
    # Put the key value pairs into a hash
    foreach(split("\n", $studies{$study})) {
      chomp;
      $_ =~ /^(.*?)=(.*)$/;
      my $key = $1;
      my $value = $2;
      $ini{$key} = $value;
    }
    # Should we include this study?
    next unless $ini{'use'} eq 'YES';
    # Create the trackDb text
    $groups .= sprintf("track %s\nsuperTrack on\ngroup %s\nshortLabel %s\nlongLabel %s\nhtml doc/%s\n\n", $study, $study, $study, $ini{'study_title'} || $study, $study);
    # Create the project description
    my $proj_desc = sprintf('ENA Project ID: <a href="http://www.ebi.ac.uk/ena/data/view/%s">%s</a><br />', $study, $study);
    $proj_desc .= sprintf('ArrayExpress ID: <a href="http://www.ebi.ac.uk/arrayexpress/experiments/%s">%s</a><br />', $ini{'ArrayExpress_ID'}, $ini{'ArrayExpress_ID'}) if $ini{'ArrayExpress_ID'};
    if($ini{'study_accession'}) {
      my $ref = get_ENA_project($ini{'study_accession'});
      $proj_desc .= "<br />$ref<br />" if $ref;
      unless($ref) {
        my $ref = get_ENA_project($study);
        $proj_desc .= "<br />$ref<br />" if $ref;
      }
    }
    if($ini{'pubmed'}) {
      my $ref = get_reference($ini{'pubmed'});
      if($ref) {
        $proj_desc .= "<br />Reference<br />$ref";
      } else {
        $proj_desc .= sprintf('<br />PubMed ID: <a href="http://europepmc.org/abstract/MED/%s">%s</a>', $ini{'pubmed'}, $ini{'pubmed'});
      }
    }
    $proj_desc =~ s/\n//g;
    mkdir "myHub/$species/doc" unless -d "myHub/$species/doc";
    open(HTMLOUT, ">myHub/$species/doc/$study.html");
    print HTMLOUT $proj_desc;
    close(HTMLOUT);
    # Get the unique sample IDs
    my @samples;
    foreach my $key (keys %ini) {
      if($key =~ /^library_sample/) {
        push(@samples, $ini{$key}) unless grep{$_ eq $ini{$key}} @samples;
      }
    }
    # Loop through each sample
    foreach my $sample (@samples) {
      warn "  -- Sample: $sample";
      $counter++;
      # Create the unique track ID
      my $track_id = sprintf("%03d", $counter) . "_" . $sample;
      my $url = sprintf("http://ngs.sanger.ac.uk/production/parasites/wormbase/RNASeq_alignments/%s/%s.bw", lc($species), $sample);
      # Create the sample description
      my $desc = sprintf(
                  'ENA Sample ID: <a href="http://www.ebi.ac.uk/ena/data/view/%s">%s</a><br />
                   %s%s%s<br />',
                $sample, $sample,
                $ini{"sample_ChEBI_ID_$sample"} ? sprintf(
                    'ChEBI: <a href="https://www.ebi.ac.uk/chebi/searchId.do?chebiId=%s">%s</a><br />',
                    $ini{"sample_ChEBI_ID_$sample"}, $ini{"sample_ChEBI_ID_$sample"}) : '',
                $ini{"sample_WormBaseLifeStage_$sample"} ? sprintf(
                    'WormBase Life Stage: <a href="http://www.wormbase.org/species/all/life_stage/%s">%s</a><br />',
                    $ini{"sample_WormBaseLifeStage_$sample"}, $ini{"sample_WormBaseLifeStage_$sample"}) : '',
                $proj_desc);
      $desc =~ s/\n//g;
      $desc .= "<br /><br />This data comes from URL: $url";
      mkdir "myHub/$species/doc" unless -d "myHub/$species/doc";
      open(HTMLOUT, ">myHub/$species/doc/$track_id.html");
      print HTMLOUT $desc;
      close(HTMLOUT);
      # Create the trackDb text
      $files .= sprintf("track %s\nparent %s\ntype bigWig\nbigDataUrl %s\nshortLabel %s\nlongLabel %s\ncolor %s\nhtml doc/%s\n\n", $track_id, $study, $url, $ini{"sample_shortLabel_$sample"}, $ini{"sample_longLabel_$sample"}, $ini{'Colour'}, $track_id);
    }
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
    my $text;
    if($result->{PROJECT}->{DESCRIPTION}) {
      my $formatted = encode_entities($result->{PROJECT}->{DESCRIPTION});
      $text = "Project Description<br />$formatted";
    } elsif($result->{STUDY}->{DESCRIPTOR}->{STUDY_DESCRIPTION}) {
      my $formatted = encode_entities($result->{STUDY}->{DESCRIPTOR}->{STUDY_DESCRIPTION});
      $text = "Study Description<br />$formatted";
    }
    return $text;
  }
}
 
