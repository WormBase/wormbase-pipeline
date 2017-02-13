#!/usr/bin/perl

use strict;

use Data::Dumper;
use File::Basename;
use File::Copy qw(copy);
use Getopt::Long;
use IO::Uncompress::Gunzip qw(gunzip);
use Log::Log4perl qw(:easy);
use Pod::Usage qw(pod2usage);

Log::Log4perl->easy_init($INFO);
my $logger = get_logger();

## Get the command line options
if(!@ARGV) {
  pod2usage(1);
  exit;
}
my ($ftp_dir, $jbrowse_path, $out_dir, $gff3_config_file, $species_name);
GetOptions(
    'ftp_dir=s'          => \$ftp_dir,
    'jbrowse_path=s'     => \$jbrowse_path,
    'out_dir=s'          => \$out_dir,
    'gff3_config_file=s' => \$gff3_config_file,
    'species=s'          => \$species_name
  );

## Load the GFF3 config into memory - we only want to do this once
my @gff3_config;
open(FILE, $gff3_config_file) or $logger->logdie("Cannot open GFF3 track configuration file: $!");
foreach(<FILE>) {
  chomp;
  my @parts = split("\t", $_);
  push(@gff3_config, {
    'type' => $parts[0],
    'category' => $parts[1],
    'trackLabel' => $parts[2],
    'trackType' => $parts[3]
  });
}
close(FILE);

## Get all the species directories from the FTP site
opendir(DIR, "$ftp_dir/species/") or $logger->logdie("Cannot open species directory: $!");
my @species = readdir(DIR);
closedir(DIR);
@species = qq($species_name) if $species_name;
for my $species (@species) {
  next if $species =~ /^\.*$/;
  my $species_name = basename($species);
  
  ## Get all the genomes for a species
  opendir(DIR, "$ftp_dir/species/$species") or $logger->logdie("Cannot open species directory: $!");
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
    opendir(DIR, "$ftp_dir/species/$species/$bioproject") or $logger->logdie("Cannot open genome directory: $!");
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

    ## Create trackList.json
    $logger->info("Creating trackList.json");
    open(TRACKJSON, ">$out_dir/$prod_name/data/trackList.json");
    print TRACKJSON '{ "include" : ["functions.conf"] }';
    close(TRACKJSON);

    ## Copy the functions file
    copy("./includes/functions.conf",  "$out_dir/$prod_name/data/functions.conf");

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
      $logger->info(sprintf("-- feature type %s", $gff3_track->{'type'}));
      (my $track_label = $gff3_track->{'trackLabel'}) =~ s/\s/_/g;
      my $sys_cmd = sprintf(
        qq(perl %s/bin/flatfile-to-json.pl --gff "%s/%s/data_files/%s" --type %s --key "%s" --trackLabel "%s" --trackType %s --metadata '{ "category": "%s", "menuTemplate" : [{ "label" : "View gene at WormBase ParaSite", "action" : "newWindow", "url" : "/Gene/Summary?g={name}" }] }' %s --out "%s/%s/data"),
        $jbrowse_path,
        $out_dir,
        $prod_name,
        $gff3_file,
        $gff3_track->{'type'},
        $gff3_track->{'trackLabel'},
        $track_label,
        $gff3_track->{'trackType'},
        $gff3_track->{'category'},
	$gff3_track->{'type'} =~ /^gene/ ? qq(--clientConfig '{ "color" : "geneColor", "label" : "geneLabel" }') : '';
        $out_dir,
        $prod_name
      );
      `$sys_cmd`;
    }
    
  }
}


## Viewable at http://<path_to_jbrowse_server>/?data=parasite/<production_name>/data/


__END__

=head1 NAME

create_species_config.pl - produce a JBrowse configuration from a WormBase ParaSite FTP directory

=head1 USAGE

  create_species_config.pl                                \
      --ftp_dir <path to filesystem location of FTP site> \
      --jbrowse_path <path to JBrowse check out>          \
      --out_dir <output directory>                        \
      --gff3_config_file <path to GFF3 mapping file>
