#!/usr/bin/env perl
use strict;
use warnings;
use SpeciesFtp;
use File::Basename;

die "Usage: PARASITE_VERSION=... PARASITE_STAGING_MYSQL=... $0 <pattern>" unless $ENV{PARASITE_VERSION} and $ENV{PARASITE_STAGING_MYSQL};
my ($species_pattern) = @ARGV;
my $sql_cmd = $ENV{PARASITE_STAGING_MYSQL};
my $dataset_names_sql = "select name, sql_name from parasite_mart_$ENV{PARASITE_VERSION}.dataset_names where src_db like \\\"%_core_$ENV{PARASITE_VERSION}_%\\\"";
open (my $dataset_names_fh, "$sql_cmd -Ne \" $dataset_names_sql \" | " ) or die $!;

my %species_by_biomart_name;
while(<$dataset_names_fh>){
  my ($biomart_name, $species) = split;
  $species_by_biomart_name{$biomart_name} = $species;
}
print STDERR sprintf("Present in BioMart: %s species\n", scalar keys  %species_by_biomart_name);
die unless %species_by_biomart_name;

my %species_with_paralogs_biomart_names;
open (my $baby_mart_tables_fh, "$sql_cmd parasite_mart_$ENV{PARASITE_VERSION} -Ne \"show tables\" | ");
while(<$baby_mart_tables_fh>){
  my ($biomart_name, $biomart_name_2) = /(.*)_gene__paralog_(.*)__dm/;
  if($biomart_name and $species_by_biomart_name{$biomart_name}){
     die $_ unless $biomart_name eq $biomart_name_2;
     $species_with_paralogs_biomart_names{$species_by_biomart_name{$biomart_name}} = $biomart_name;
  }
}
print STDERR sprintf("Has paralog tables in BioMart: %s species\n", scalar keys  %species_with_paralogs_biomart_names);

my %species_with_orthologs_biomart_names;

open (my $merged_mart_tables_fh, "$sql_cmd parasite_mart_$ENV{PARASITE_VERSION}_merged -Ne \"show tables\" | ");
while(<$merged_mart_tables_fh>){
  my ($biomart_name) = /wbps_gene__homolog_(.*)__dm/;
  if($biomart_name and $species_by_biomart_name{$biomart_name}){
     $species_with_orthologs_biomart_names{$species_by_biomart_name{$biomart_name}} = $biomart_name;
  }
}
print STDERR sprintf("Has ortholog tables in merged BioMart: %s species\n", scalar keys  %species_with_orthologs_biomart_names);
die unless %species_with_orthologs_biomart_names;

my $GET_PARALOGS_TEMPLATE = <<'EOF';
SELECT stable_id_1023 AS gene_id,
       stable_id_4016_r2 AS paralog_gene_id,
       perc_id_4015_r1 AS target_identity,
       perc_id_4015 AS query_identity
FROM parasite_mart_{{PARASITE_VERSION}}.{{SPECIES_BIOMART_NAME}}_gene__paralog_{{SPECIES_BIOMART_NAME}}__dm
JOIN parasite_mart_{{PARASITE_VERSION}}.{{SPECIES_BIOMART_NAME}}_gene__gene__main USING (gene_id_1020_key)
WHERE stable_id_4016_r2 IS NOT NULL
ORDER BY gene_id, paralog_gene_id
EOF

my $GET_ORTHOLOGS_TEMPLATE = <<'EOF';
SELECT h.stable_id_4016_r2 AS gene_id,
       d.species_name AS ortholog_species_name,
       g.stable_id_1023 AS ortholog_gene_id,
       h.perc_id_4015 AS target_identity,
       h.perc_id_4015_r1 AS query_identity,
       IFNULL(g.description_1020, "") AS ortholog_description
FROM parasite_mart_{{PARASITE_VERSION}}_merged.wbps_gene__homolog_{{SPECIES_BIOMART_NAME}}__dm h
JOIN parasite_mart_{{PARASITE_VERSION}}_merged.wbps_gene__gene__main g ON (h.gene_id_1020_key = g.gene_id_1020_key)
JOIN parasite_mart_{{PARASITE_VERSION}}.dataset_names d ON (g.species_id_1010_key = d.name)
WHERE stable_id_4016_r2 IS NOT NULL
ORDER BY gene_id,
         ortholog_species_name,
         ortholog_gene_id
EOF

sub sql {
  my ($sql_template, $species_biomart_name) = @_;
  $sql_template =~ s/{{PARASITE_VERSION}}/$ENV{PARASITE_VERSION}/g;
  $sql_template =~ s/{{SPECIES_BIOMART_NAME}}/$species_biomart_name/g;
  $sql_template =~  s/"/\\"/g;
  chomp $sql_template;
  return $sql_template;
}
sub dump_sql_to_path {
  my ($sql, $path) = @_;
  unlink $path if -f $path;
  unlink "$path.gz" if -f "$path.gz";
  die "No parent directory: $path" unless -d (dirname $path);
  my $cmd = "$sql_cmd -e \" $sql \" > $path ";
  my $msg = `$cmd`;
  die $msg if $msg =~ /error/i;
  system("gzip $path") and die "Could not gzip: $path";
}
for my $species (sort keys %species_with_paralogs_biomart_names){
  next if $species_pattern and $species !~ /$species_pattern/;
  my $species_biomart_name = $species_with_paralogs_biomart_names{$species};
  my $path = SpeciesFtp->current_staging->path_to($species, "paralogs.tsv", 0);
  print STDERR localtime . " $species\t$path.gz\n";
  if (!-e "$path.gz") {
    dump_sql_to_path(
        sql($GET_PARALOGS_TEMPLATE, $species_biomart_name),
        $path,
    );
  }
}
for my $species (sort keys %species_with_orthologs_biomart_names){
  next if $species_pattern and $species !~ /$species_pattern/;
  my $species_biomart_name = $species_with_orthologs_biomart_names{$species};
  my $path = SpeciesFtp->current_staging->path_to($species, "orthologs.tsv", 0);
  print STDERR localtime ." $species\t$path.gz\n";
  if (!-e "$path.gz") {
    dump_sql_to_path(
      sql($GET_ORTHOLOGS_TEMPLATE, $species_biomart_name),
      $path,
    );
  }
}
