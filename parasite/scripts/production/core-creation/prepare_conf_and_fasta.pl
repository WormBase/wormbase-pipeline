#!/usr/bin/env perl
# Sorts out FASTA split-or-not and conf template
# Creates the conf in a given directory
# Alternative conf template can be passed on STDIN
use strict;
use warnings;
use YAML;
use CoreCreation::Fasta;
use File::Spec;
use File::Basename;
use File::Slurp;

my $conf = Load(do {local $/; <>} || "--- {}");
die "usage: PARASITE_DATA=... $0 < conf" unless $conf;

my ($data_dir_name, @others) = keys %{$conf};
die "conf format: {species => {...}}" unless $data_dir_name and not @others and $data_dir_name =~ /[a-z]+_[a-z]+_[a-z0-9]+/;
my $data_dir_path = join ("/", $ENV{PARASITE_DATA}, $data_dir_name);
my $conf_path = File::Spec->catfile($data_dir_path, "$data_dir_name.conf");

if ( not -s $conf_path or $ENV{REDO_FASTA}) {
  $conf = $conf->{$data_dir_name};
  $conf->{gff3} //=  File::Spec->catfile($data_dir_path, "$data_dir_name.gff3");
  die "Expected gff at: ".$conf->{gff3} unless -f $conf->{gff3} and File::Spec->file_name_is_absolute($conf->{gff3});
  my $check_sources_column = "grep -c $conf->{gff_sources} $conf->{gff3}";
  die "Failed: $check_sources_column" unless 0 < `$check_sources_column`;
  my $fasta_path = File::Spec->catfile($data_dir_path,"$data_dir_name.fa");
  die "Expected fasta at: $fasta_path" unless -f $fasta_path;
  $conf->{toplevel} = "scaffold";
  my $fasta = CoreCreation::Fasta->new($fasta_path);
  if($fasta->needs_contig_structure and not $ENV{SKIP_SPLIT_FASTA}){
    (my $split_fasta_path = $fasta_path) =~ s/.fa/.split.fa/;
    (my $agp_path = $fasta_path) =~ s/.fa/.toplevel.agp/;
    $fasta->split(fasta => $split_fasta_path, agp => $agp_path);
    $conf->{fasta} = $split_fasta_path;
    $conf->{agp} = $agp_path;
    $conf->{seqlevel} = "contig";
  } else {
    $conf->{fasta} = $fasta_path;
    $conf->{seqlevel} = "scaffold";
  }
  my $mito = $fasta->mito;
  $conf->{mitochondrial} = $mito if $mito;
  
  
  open(FH, '>', $conf_path) or die $!;
  print FH Dump({$data_dir_name => $conf});
  close(FH);
}

my $text = File::Slurp::read_file($conf_path);
if ($text =~ /\?/) {
  die "$conf_path: complete the config file and remove all the ?s!";
} 
print Dump (YAML::LoadFile($conf_path)); 
