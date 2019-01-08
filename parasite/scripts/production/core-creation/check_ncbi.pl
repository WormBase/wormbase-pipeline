#!/usr/bin/env perl
use strict;
use warnings;
use YAML;
use ProductionMysql;
use Data::Compare;
use Carp;
my $production_mysql = ProductionMysql->staging;
# Get these by BioProject
# Doesn't work when BioProject changed


my %assemblies = ( %{get_assemblies_for_taxon (6157)}, %{get_assemblies_for_taxon (6231)}); 

my %report = (
  PARASITE_ONLY => [],
  NCBI_ONLY => {},
  MATCHING_METADATA => {},
  MATCHING_SCAFFOLD_NAMES => {},
  MATCHING_SCAFFOLD_LENGTHS => {},
  MISMATCHED => {},
);

my @core_dbs = $production_mysql->core_databases;
my @core_dbs_parasite_only;
for my $core_db ($production_mysql->core_databases ){
   my ($spe, $cies, $bioproject) = split "_", $core_db;
   my $species = "${spe}_${cies}";
   if ($assemblies{$species}{$bioproject}){
      my $accession = $production_mysql->meta_value($core_db, "assembly.accession") // "";
      my $name = $production_mysql->meta_value($core_db, "assembly.name") // "";
      $assemblies{$species}{$bioproject}{core_db} =  {
         name => $core_db,
         assembly_accession => $accession,
         assembly_name => $name,
     };
   } else {
      push @{$report{PARASITE_ONLY}}, $core_db;
   }
}

for my $species (keys %assemblies) {
  BIOPROJECT:
  for my $bioproject (keys %{$assemblies{$species}}){
    my @ncbi_assemblies = @{$assemblies{$species}{$bioproject}{ncbi}};
    my $core_db = $assemblies{$species}{$bioproject}{core_db};
    unless ($core_db) {
       $report{NCBI_ONLY}{$species}{$bioproject} = \@ncbi_assemblies;
       next BIOPROJECT;
    }
    my ($ncbi_assembly, @others) = grep {same_assembly_names($core_db, $_)} @ncbi_assemblies;
    confess @ncbi_assemblies, $core_db if @others;
    if ($ncbi_assembly) {
       my @other_assemblies = grep {not same_assembly_names($ncbi_assembly, $_ )} @ncbi_assemblies;
       $report{MATCHING_METADATA}{$species}{$bioproject} = {%{$ncbi_assembly}, core_db => $core_db->{name}, other_assemblies => \@other_assemblies };
       next BIOPROJECT;
    }
    my $scaffold_lengths_core = get_scaffold_lengths_for_core_db($core_db->{name});
    for my $ncbi_assembly (@ncbi_assemblies){
       my @other_assemblies = grep {not same_assembly_names($ncbi_assembly, $_ )} @ncbi_assemblies;
       my $scaffold_lengths_ncbi = get_scaffold_lengths_from_report_url($ncbi_assembly->{report_url});
       if(Data::Compare::Compare ($scaffold_lengths_ncbi, $scaffold_lengths_core)){
           $report{MATCHING_SCAFFOLD_NAMES}{$species}{$bioproject} = {%{$ncbi_assembly}, core_db => $core_db->{name}, other_assemblies => \@other_assemblies };
           next BIOPROJECT;
       } elsif (Data::Compare::Compare([sort values %{$scaffold_lengths_ncbi}], [sort values %{$scaffold_lengths_core}])){
           $report{MATCHING_SCAFFOLD_LENGTHS}{$species}{$bioproject} = {%{$ncbi_assembly}, core_db => $core_db->{name}, other_assemblies => \@other_assemblies };
           next BIOPROJECT;
       }
    }
    $report{MISMATCHED}{$species}{$bioproject} = { %{$core_db}, ncbi => \@ncbi_assemblies};
  }
}
print Dump \%report;
sub same_assembly_names {
  my ($o1, $o2) = @_;
  return $o1->{assembly_accession} eq $o2->{assembly_accession} && $o1->{assembly_name} eq $o2->{assembly_name};
}

sub get_scaffold_lengths_from_report_url {
  my ($report_url) = @_;
  my %result;
  open(my $fh, "curl -s $report_url | grep -v '^#' | cut -f 5,9 | ");
  while (<$fh>){
    chomp;
    my ($scaffold, $length) = split "\t";
    $result{$scaffold} = $length;
  }
  return \%result;
}
sub get_scaffold_lengths_for_core_db {
  my ($core_db) = @_;
  my %result;
  for my $slice (@{ $production_mysql->adaptor($core_db, "Slice")->fetch_all("toplevel") }){
    my ($x, @xs) = map { $_->name} grep {$_->dbname eq "INSDC"} @{ $slice->get_all_synonyms };
    $result{$x || ($slice->seq_region_name) } = $slice->seq_region_length ; 
  } 
  return \%result;
}

sub get_assemblies_for_taxon {
  my ($taxon) = @_;
  my %result;
  local $/ = "--";
  open(my $fh, "curl -s 'https://www.ncbi.nlm.nih.gov/assembly/organism/browse-organism/$taxon/all/?p\$l=Ajax&page=1&pageSize=1000&sortCol=4&sortDir=-1' | grep -B 2 href |");
  while(<$fh>){
    my ($species) = /<td>(.*?) \(/;
    $species = lc $species;
    $species =~ tr/ /_/;

    my ($accession, $name) = /<a href="https:\/\/www.ncbi.nlm.nih.gov\/assembly\/(.*?)\/">(.*?)</;
    my $doc = `curl -s "https://www.ncbi.nlm.nih.gov/assembly/$accession?report=xml&format=text" `;
    my ($bioproject) = $doc =~ m{BioprojectAccn.*(PR\w+\d+).*\/BioprojectAccn};
    $bioproject = lc $bioproject;
    my ($report_url) = $doc =~ m{(ftp://.*assembly_report.txt)};
    push @{$result{$species}{$bioproject}{ncbi}}, {
      assembly_accession => $accession,
      assembly_name => $name,
      bioproject => $bioproject,
      report_url => $report_url,
    };
  }
  return \%result;
}
