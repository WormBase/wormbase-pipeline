#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long qw(GetOptions);
use HTML::Entities;
use JSON;
use Log::Log4perl qw(:easy);
use LWP::UserAgent;
use Pod::Usage qw(pod2usage);
use XML::Simple;

Log::Log4perl->easy_init($INFO);
my $logger = get_logger();

my $in;
my $out;
my $track_hub = 0;
my $jbrowse = 0;

GetOptions('in=s'      => \$in,
           'out=s'     => \$out,
           'track_hub' => \$track_hub,
           'jbrowse'   => \$jbrowse
          );
if(!$in || !$out || (!$track_hub && !$jbrowse)) {
  pod2usage(1);
  exit;
}
          
die "Can only output one of track hub or JBrowse at any time" if $track_hub && $jbrowse;

my @species_list = `ls $in/*.ini | xargs -n1 basename`;
foreach(@species_list) {
  chomp;
  $_ =~ s/\.ini//;
}

my $counter = 0;

# Create the hub.txt file
mkdir $out unless -d $out;
if($track_hub) {
  open(OUTFILE, ">$out/hub.txt");
  print OUTFILE "hub WBPS-RNASeq\nshortLabel RNA-Seq Alignments\nlongLabel RNA-Seq Alignments for WormBase ParaSite\ngenomesFile genomes.txt\nemail parasite-help\@sanger.ac.uk\n";
  close(OUTFILE);

  open(OUTFILE, ">$out/genomes.txt");
  close(OUTFILE);
}

foreach my $in_file (@species_list) {
 
  $logger->info("Species: $in_file");
  $counter = 0;

  my $file = "$in/$in_file.ini";
  open(INFILE, $file);

  # Get the BioProject
  my $bioproject;
  foreach(<INFILE>) {
    chomp;
    $_ =~ /^general_bioproject=(.*?)$/;
    $bioproject = $1 if $1;
    last if $bioproject;
  }
  $logger->info("Skipping as no BioProject") && next unless $bioproject;
  my $species = $in_file . "_" . lc($bioproject);
  $logger->info("BioProject: $bioproject");

  if($track_hub) {
    # Write to genomes.txt
    open(OUTFILE, ">>$out/genomes.txt");
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
 
    mkdir "$out/$species" unless -d "$out/$species";
    open(OUTFILE, ">$out/$species/trackDb.txt");
  }
  
  if($jbrowse) {
    my $species_lc = lc $species;
    mkdir "$out/$species_lc" unless -d "$out/$species_lc";
    mkdir "$out/$species_lc/data" unless -d "$out/$species_lc/data";
    open(OUTFILE, ">>$out/$species_lc/data/tracks.conf");
  }

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
    next if $study =~ /^general$/i;
    $logger->info("-- Study: $study");
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
    if($track_hub) {
      mkdir "$out/$species/doc" unless -d "$out/$species/doc";
      open(HTMLOUT, ">$out/$species/doc/$study.html");
      print HTMLOUT $proj_desc;
      close(HTMLOUT);
    }
    # Get the unique sample IDs
    my @samples;
    foreach my $key (keys %ini) {
      if($key =~ /^library_sample/) {
        push(@samples, $ini{$key}) unless grep{$_ eq $ini{$key}} @samples;
      }
    }
    # Put the tracks into some logical order
    my %descriptions;
    foreach my $sample (@samples) {
      $descriptions{$sample} = $ini{"sample_shortLabel_$sample"};
    }
    my @samples_ordered;
    for my $k (sort {$descriptions{$a} cmp $descriptions{$b}} keys %descriptions) {
      push(@samples_ordered, $k);
    }
    # Loop through each sample
    foreach my $sample (@samples_ordered) {
      $logger->info("  -- Sample: $sample");
      # Skip samples that do not have a description
      next unless $ini{"sample_shortLabel_$sample"} && $ini{"sample_longLabel_$sample"};
      $counter++;
      # Should the track be enabled?
      $ini{"display"} = 'off' unless $ini{"display"};
      my $display = $ini{"display"} =~ /^on$/i ? 'full' : 'hide';
      # Create the unique track ID
      my $track_id = sprintf("%03d", $counter) . "_" . $sample;
      my $url = sprintf("http://ngs.sanger.ac.uk/production/parasites/wormbase/RNASeq_alignments/%s/%s.bw", lc($species), $sample);
      # Get the life stage URLs
      my $life_stage = '';
      if($ini{"sample_WormBaseLifeStage_$sample"}) {
        my @stages = split(',', $ini{"sample_WormBaseLifeStage_$sample"});
        my @urls;
        foreach my $stage (@stages) {
          $stage =~ s/ //;
          push @urls, sprintf('<a href="%s">%s</a>', get_life_stage_url($stage), $stage);
        }
        $life_stage = sprintf("WormBase Life Stage: %s<br />", join('; ', @urls));
      }
      # Create the sample description
      my $desc = sprintf(
                  'ENA Sample ID: <a href="http://www.ebi.ac.uk/ena/data/view/%s">%s</a><br />
                   %s%s%s<br />',
                $sample, $sample,
                $ini{"sample_ChEBI_ID_$sample"} ? sprintf(
                    'ChEBI: <a href="https://www.ebi.ac.uk/chebi/searchId.do?chebiId=%s">%s</a><br />',
                    $ini{"sample_ChEBI_ID_$sample"}, $ini{"sample_ChEBI_ID_$sample"}) : '',
                $life_stage,
                $proj_desc);
      $desc =~ s/\n//g;
      $desc =~ s/<br \/>/\n<br \/>\n/g;
      my $ftp = sprintf("ftp://ngs.sanger.ac.uk/production/parasites/wormbase/RNASeq_alignments/%s", lc($species));
      $desc .= sprintf('<br /><br />This data comes from URL: <a href="%s">%s</a><br />Download data (BAM and BigWig): <a href="%s">%s</a>', $url, $url, $ftp, $ftp);
      if($track_hub) {
        mkdir "$out/$species/doc" unless -d "$out/$species/doc";
        open(HTMLOUT, ">$out/$species/doc/$track_id.html");
        print HTMLOUT $desc;
        close(HTMLOUT);
        # Create the trackDb text
        $files .= sprintf("track %s\nparent %s\ntype bigWig\nbigDataUrl %s\nshortLabel %s\nlongLabel %s\ncolor %s\nhtml doc/%s\nvisibility %s\n\n",
                          $track_id,
                          $study,
                          $url,
                          $ini{"sample_shortLabel_$sample"},
                          $ini{"sample_longLabel_$sample"},
                          $ini{'Colour'} || "0,0,0",
                          $track_id,
                          $display
                        );
      }
      if($jbrowse) {
        $files .= sprintf("[tracks.%s]\nstoreClass = JBrowse/Store/SeqFeature/BigWig\ntype = JBrowse/View/Track/Wiggle/XYPlot\nurlTemplate = %s\nkey = %s\ncategory = RNA-Seq/%s\nautoscale = local\nyScalePosition = right\nstyle.pos_color = rgb(%s)\n\n",
                          $track_id,
                          $url,
                          $ini{"sample_longLabel_$sample"},
                          $study,
                          $ini{'Colour'} || "0,0,0"
                        );                     
      }
    }
  }

  if($track_hub) {
    print OUTFILE $groups, "\n";
    print OUTFILE $files;
    close(OUTFILE);
  }
  if($jbrowse) {
    print OUTFILE $files;
    close(OUTFILE);
  }

  close(INFILE);

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

sub get_life_stage_url {
  my ($id) = @_;
  my $url = "http://www.wormbase.org/rest/widget/life_stage/$id/overview?content-type=application%2Fjson";
  my $ua = LWP::UserAgent->new();
  my $response = $ua->get($url);
  if ($response->is_success) {
    return "http://www.wormbase.org/species/all/life_stage/$id";
  } else {
    $id =~ s/:/_/;
    return "http://www.ebi.ac.uk/ols/beta/ontologies/wbls/terms?iri=http://purl.obolibrary.org/obo/$id";
  }
}
 

__END__

=head1 NAME

create_rnaseq_track_conf.pl - takes the RNA-seq ini configuration files to produce either a complete track hub or JBrowse conf files

=head1 USAGE

  create_rnaseq_track_conf.pl           \
      --in <path to ini files>          \
      --out <path to track hub output>  \
      [ --jbrowse ]                     \
      [ --track_hub ]
