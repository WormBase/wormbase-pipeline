#!/usr/bin/perl

use strict;

use File::Basename;
use File::Copy qw(copy);
use IO::Uncompress::Gunzip qw(gunzip);
use Log::Log4perl qw(:easy);

## TODO: Take arguments from the command line
my $ftp_dir = "/ebi/ftp/pub/databases/wormbase/parasite/releases/current";
my $jbrowse_path = "/homes/bbolt/dev/jbrowse/JBrowse-1.12.1";
my $out_dir = "/homes/bbolt/dev/jbrowse/JBrowse-1.12.1/parasite";
my $gff3_config_file = "/homes/bbolt/dev/wormbase-pipeline/parasite/scripts/production/jbrowse/gff3_tracks.tsv";

Log::Log4perl->easy_init($INFO);
my $logger = get_logger();

## Load the GFF3 config into memory - we only want to do this once
my @gff3_config;
open(FILE, $gff3_config_file) or die "Cannot open GFF3 track configuration file";
foreach(<FILE>) {
  chomp;
  my @parts = split("\t", $_);
  push(@gff3_config, {
    'type' => $parts[0],
    'category' => $parts[1],
    'trackLabel' => $parts[2],
    'urltemplate' => $parts[3]
  });
}
close(FILE);

## Get all the species directories from the FTP site
opendir(DIR, "$ftp_dir/species/") or die "Cannot open species directory";
my @species = readdir(DIR);
closedir(DIR);
### TODO: For testing, we're just limiting to one species for now
@species = qw(schistosoma_mansoni);
###
for my $species (@species) {
  next if $species =~ /^\.*$/;
  my $species_name = basename($species);
  
  ## Get all the genomes for a species
  opendir(DIR, "$ftp_dir/species/$species") or die "Cannot open species directory";
  my @genomes = readdir(DIR);
  closedir(DIR);
  foreach my $genome (@genomes) {
    next if $genome =~ /^\.*$/;
    my $bioproject = basename($genome);
    my $prod_name = lc sprintf("%s_%s", $species_name, $bioproject);
    $logger->info("Working on $prod_name");
    mkdir "$out_dir/$prod_name";
    mkdir "$out_dir/$prod_name/data_files";
    
    ## Get a list of the available files for this genome
    opendir(DIR, "$ftp_dir/species/$species/$bioproject") or die "Cannot open genome directory";
    my @files = readdir(DIR);
    closedir(DIR);

    ## Find the name of the files we are interested in
    my ($fasta_file, $gff3_file);
    foreach(@files) {
      $_ =~ s/.gz$//;
      if($_ =~ /genomic\.fa$/) {
        $fasta_file = $_;
      }
      if($_ =~ /annotations\.gff3$/) {
        $gff3_file = $_;
      }
    }

    ## Process the FASTA
    $logger->info("Converting $fasta_file to JSON");
    copy("$ftp_dir/species/$species/$bioproject/$fasta_file.gz", "$out_dir/$prod_name/data_files/$fasta_file.gz");
    gunzip("$out_dir/$prod_name/data_files/$fasta_file.gz", "$out_dir/$prod_name/data_files/$fasta_file");
    `perl $jbrowse_path/bin/prepare-refseqs.pl --fasta $out_dir/$prod_name/data_files/$fasta_file --out $out_dir/$prod_name/data`;
    
    ## Process the GFF3 - note there are many feature types here which are extracted into single tracks, so we load them from an external file to make them configurable
    $logger->info("Converting $gff3_file to JSON");
    copy("$ftp_dir/species/$species/$bioproject/$gff3_file.gz", "$out_dir/$prod_name/data_files/$gff3_file.gz");
    gunzip("$out_dir/$prod_name/data_files/$gff3_file.gz", "$out_dir/$prod_name/data_files/$gff3_file");
    foreach my $gff3_track (@gff3_config) {
      warn sprintf("-- feature type %s", $gff3_track->{'type'});
      my $sys_cmd = sprintf(
        qq(perl %s/bin/flatfile-to-json.pl --gff "%s/%s/data_files/%s" --type %s --trackLabel "%s" --urltemplate "%s" --metadata '{"category": "%s"}' --out "%s/%s/data"),
        $jbrowse_path,
        $out_dir,
        $prod_name,
        $gff3_file,
        $gff3_track->{'type'},
        $gff3_track->{'trackLabel'},
        $gff3_track->{'urltemplate'},
        $gff3_track->{'category'},
        $out_dir,
        $prod_name
      );
      `$sys_cmd`;
    }
    
  }
}


## Viewable at http://<path_to_jbrowse_server>/?data=parasite/schistosoma_mansoni_prjea36577/data/
