#!/usr/bin/env perl

use YAML;
use CoreCreation::Fasta;
use File::Spec;
use File::Basename;
use File::Slurp;
my $data_dir_path = $ARGV[0];
die "Usage: $0 <data dir>" unless -d $data_dir_path;
$data_dir_path =  File::Spec->rel2abs($data_dir_path) unless File::Spec->file_name_is_absolute($data_dir_path);
my $data_dir_name = File::Basename::basename($data_dir_path); 
my $conf_path = File::Spec->catfile($data_dir_path, "$data_dir_name.conf");

while (true) {
  if ( not -f $conf_path) {
    my $conf = Load(do {local $/; <DATA>});
    my $gff_path = File::Spec->catfile($data_dir_path, "$data_dir_name.gff3");
    die "Expected gff at: $gff_path" unless -f $gff_path;
    my $fasta_path = File::Spec->catfile($data_dir_path,"$data_dir_name.fa");
    die "Expected fasta at: $fasta_path" unless -f $fasta_path;
    $conf->{gff3} = $gff_path;
    $conf->{toplevel} = "scaffold";
    my $fasta = CoreCreation::Fasta->new($fasta_path);
    if($fasta->needs_contig_structure){
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
    
    
    open(FH, '>', $conf_path) or die $!;
    print FH Dump({ $data_dir_name => $conf});
    close(FH);
  }
  
  my $text = File::Slurp::read_file($conf_path);
  if ($text =~ /\?/) {
    system("vi $conf_path") and die; #:cq in vim to escape this
  } else {
    last;
  }  
}
print Dump (YAML::LoadFile($conf_path)); 

__DATA__
---
taxon_id: ?
assembly_version: ?
#look through the GFF's source column. worm_lite.pl will ignore things not mentioned here.
gff_sources: ?
meta:
 assembly.accession: ?
 provider.name: ?
 provider.url: ? 
 species.biosample: ?
# Consult WormBook if needed. Remove for Platyhelminthes.
 species.nematode_clade: ?
