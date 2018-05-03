#!/usr/bin/env perl
use strict;
use YAML;

my @wormbase_core_species = qw(
  caenorhabditis_elegans_prjna13758
  caenorhabditis_briggsae_prjna10731
  caenorhabditis_brenneri_prjna20035
  caenorhabditis_remanei_prjna53967
  caenorhabditis_japonica_prjna12591
  pristionchus_pacificus_prjna12644
  brugia_malayi_prjna10729
  onchocerca_volvulus_prjeb513
  strongyloides_ratti_prjeb125
  trichuris_muris_prjeb126
);
#Split by comma if there are more selenoproteins.
my $seleno_proteins = {
  caenorhabditis_elegans_prjna13758 => "WBGene00015553",
  caenorhabditis_briggsae_prjna10731 => "WBGene00028139",
  caenorhabditis_brenneri_prjna20035 => "WBGene00158831",
  caenorhabditis_remanei_prjna53967 => "WBGene00068657",
  caenorhabditis_japonica_prjna12591 => "WBGene00122465",
  brugia_malayi_prjna10729 => "WBGene00222286",
  onchocerca_volvulus_prjeb513 => "WBGene00241445",
};
# Could be incomplete! This is per assembly, we rarely change assemblies, and usually not through this script.
my $mitochondrial_seq_regions = {
  trichuris_muris_prjeb126 => "TMUE_MITO"
}; 

sub wormbase_ftp_dir {
  my ($species, $wormbase_version) = @_;
  my ($spe, $cies, $bioproject) = split /_/, $species;
  return join("/", "/nfs/ftp/pub/databases/wormbase/releases", "WS$wormbase_version","species", lc((substr $spe, 0, 1 ) . "_" . $cies) , uc($bioproject));
}
sub wormbase_gff3 {
  my ($species, $wormbase_version) = @_;
  my ($spe, $cies, $bioproject) = split /_/, $species;
  return join ("/", &wormbase_ftp_dir(@_), join(".", lc((substr $spe, 0, 1 ) . "_" . $cies), uc($bioproject), "WS$wormbase_version", "annotations.gff3.gz" ));
}
sub wormbase_fasta {
  my ($species, $wormbase_version) = @_;
  my ($spe, $cies, $bioproject) = split /_/, $species;
  return join ("/", &wormbase_ftp_dir(@_), join(".", lc((substr $spe, 0, 1 ) . "_" . $cies), uc($bioproject), "WS$wormbase_version", "genomic.fa.gz" ));
}
sub wormbase_core_db_name {
  my ($species, $wormbase_version, $parasite_version, $ensembl_version) = @_;
  my ($spe, $cies, $bioproject) = split /_/, $species;
  return join("_", $spe, $cies, $bioproject, "core", $parasite_version, $ensembl_version, $wormbase_version);
}

sub config {
  my ($species, @__) = @_;
  my $result = {
    gff3 => &wormbase_gff3(@_),
    fasta => &wormbase_fasta(@_),
    core_database => {
      dbname => &wormbase_core_db_name(@_),
    },
  };
  my $s = $seleno_proteins->{$species};
  $result->{"seleno"} = $s if $s;
  
  my $m = $mitochondrial_seq_regions->{$species};
  $result->{"mitochondrial"} = $m if $m;
  return $result;
}

@ARGV = @wormbase_core_species if not @ARGV;
my %species; 
for my $arg (@ARGV){
  for (@wormbase_core_species) {
    $species{$_}++ if /$arg/; 
  }
}
my $wormbase_version = $ENV{WORMBASE_VERSION};
my $parasite_version = $ENV{PARASITE_VERSION};
my $ensembl_version = $ENV{ENSEMBL_VERSION};
die "Usage: WORMBASE_VERSION=265 PARASITE_VERSION=11 ENSEMBL_VERSION=92 $0 <arg patterns, empty for all species>" 
  unless %species and $wormbase_version and $parasite_version and $ensembl_version;
for my $species (keys %species) {
  $species{$species} = &config($species, $wormbase_version, $parasite_version, $ensembl_version);
}

print Dump \%species;