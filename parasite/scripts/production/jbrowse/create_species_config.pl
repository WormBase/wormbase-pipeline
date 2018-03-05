#!/usr/bin/perl

use strict;

use Data::Dumper;
use File::Basename;
use File::Copy qw(copy);
use Getopt::Long;
use IO::Uncompress::Gunzip qw(gunzip);
use Log::Log4perl qw(:easy);
use Pod::Usage qw(pod2usage);
use FindBin qw($Bin);
Log::Log4perl->easy_init($INFO);
my $logger = get_logger();

## Get the command line options
if(!@ARGV) {
  pod2usage(1);
  exit;
}
my ($ftp_dir, $jbrowse_path, $out_dir, $species_name);
GetOptions(
    'ftp_dir=s'           => \$ftp_dir,
    'jbrowse_path=s'      => \$jbrowse_path,
    'out_dir=s'           => \$out_dir,
    'species=s'           => \$species_name
  );

## Load the GFF3 config into memory - we only want to do this once
my @gff3_config;
open(FILE, "$Bin/gff3_tracks.tsv") or $logger->logdie("Cannot open GFF3 track configuration file: $!");
foreach(<FILE>) {
  chomp;
  my @parts = split("\t", $_);
  push(@gff3_config, {
    'feature' => $parts[0],
    'type' => $parts[1],
    'category' => $parts[2],
    'trackLabel' => $parts[3],
    'trackType' => $parts[4]
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
    mkdir "$out_dir/$prod_name/data";
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
    print TRACKJSON '{ "include" : ["functions.conf"], "tracks" : [ ] }';
    close(TRACKJSON);
    
    ## Copy the functions file
    copy("$Bin/includes/functions.conf",  "$out_dir/$prod_name/data/functions.conf");
    
    ## Get the FASTA
    $logger->info("Retrieving and unzipping $fasta_file.gz");
    copy("$ftp_dir/species/$species/$bioproject/$fasta_file.gz", "$out_dir/$prod_name/data_files/$fasta_file.gz");
    gunzip("$out_dir/$prod_name/data_files/$fasta_file.gz", "$out_dir/$prod_name/data_files/$fasta_file");
    
    ## Process the FASTA
    $logger->info("Converting $fasta_file to JSON");
    `perl $jbrowse_path/bin/prepare-refseqs.pl --fasta $out_dir/$prod_name/data_files/$fasta_file --out $out_dir/$prod_name/data --compress`;
    
    ## Get the GFF3
    $logger->info("Retrieving and unzipping $gff3_file.gz");
    copy("$ftp_dir/species/$species/$bioproject/$gff3_file.gz", "$out_dir/$prod_name/data_files/$gff3_file.gz");
    gunzip("$out_dir/$prod_name/data_files/$gff3_file.gz", "$out_dir/$prod_name/data_files/$gff3_file");
    
    ## Create a stripped down version of the GFF3 to speed up processing
    $logger->info("Creating stripped down GFF3 $gff3_file.short");
    my @all_gff3_features = map { split(",", $_->{'feature'}) } @gff3_config;
    open(GFF3IN, "$out_dir/$prod_name/data_files/$gff3_file") or $logger->die("Could not open GFF3\n");;
    open(GFF3OUT, ">$out_dir/$prod_name/data_files/$gff3_file.short") or $logger->die ("Could not open outfile for writing\n");
    while(<GFF3IN>) {
      my @gff_elements = split("\t", $_);
      print GFF3OUT $_ if($gff_elements[1] ~~ \@all_gff3_features);
    }
    close(GFF3IN);
    close(GFF3OUT);

    ## Process the GFF3 by converting to JSON
    $logger->info("Converting $gff3_file.short to JSON");
    foreach my $gff3_track (@gff3_config) {
      $logger->info(sprintf("-- feature type %s", $gff3_track->{'type'}));
      (my $track_label = $gff3_track->{'trackLabel'}) =~ s/\s/_/g;
      my $sys_cmd = sprintf(
        qq(perl %s/bin/flatfile-to-json.pl --gff "%s/%s/data_files/%s.short" --type %s --key "%s" --trackLabel "%s" --trackType %s --nameAttributes "name,id,alias,locus" --metadata '{ "category": "%s", "menuTemplate" : [{ "label" : "View gene at WormBase ParaSite", "action" : "newWindow", "url" : "/Gene/Summary?g={name}" }] }' %s --out "%s/%s/data" --compress),
        $jbrowse_path,
        $out_dir,
        $prod_name,
        $gff3_file,
        $gff3_track->{'type'},
        $gff3_track->{'trackLabel'},
        $track_label,
        $gff3_track->{'trackType'},
        $gff3_track->{'category'},
	$gff3_track->{'type'} =~ /^gene/ ? qq(--clientConfig '{ "color" : "{geneColor}", "label" : "{geneLabel}" }') : '',
        $out_dir,
        $prod_name
      );
      `$sys_cmd`;
    }
    
    ## Create the search index
    $logger->info("Running generate-names.pl to index feature names for search");
    my $sys_cmd = sprintf('perl %s/bin/generate-names.pl --out %s/%s/data --compress --completionLimit 0', $jbrowse_path, $out_dir, $prod_name);
    `$sys_cmd`;

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
