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

my $stdin;
{
  local $/;
  $stdin = <>;
}

my %assemblies;
for my $doc (Load($stdin)){

    my $bioproject = $doc->{GB_BioProjects}{Bioproj}{BioprojectAccn} // "";
    
    $bioproject = lc $bioproject;

    my $species = $doc->{SpeciesName};
    $species =~ s/sp. //;
    $species =~ tr/ /_/;
    $species = lc $species;

    push @{$assemblies{$species}{$bioproject}{ncbi}}, {
      assembly_accession => $doc->{AssemblyAccession},
      assembly_name => $doc->{AssemblyName},
      assembly_url => "https://www.ncbi.nlm.nih.gov/assembly/".$doc->{AssemblyAccession},
      report_url =>  $doc->{FtpPath_Assembly_rpt},
      submission_date => $doc->{SubmissionDate},
      provider_name => $doc->{SubmitterOrganization},
    };
}

for my $core_db ($production_mysql->core_databases ){
   my ($spe, $cies, $bioproject) = split "_", $core_db;
   next if $bioproject eq "core";
   my $species = "${spe}_${cies}";
   $assemblies{$species}{$bioproject}{core_db}{name} = $core_db;
}

my %report = (
  PARASITE_ONLY => {},
  NCBI_ONLY => {},
  MATCHING_METADATA => {},
  MATCHING_SCAFFOLD_NAMES => {},
  MATCHING_SCAFFOLD_LENGTHS => {},
  MISMATCHED => {},
  MATCH_NOT_INTENDED => {},
);

my %match_not_intended = (
  ancylostoma_caninum_prjna72585 => "WashU 50 helminths",
  ascaris_suum_prjna80881 => "Gasser lab A. suum Australian isolate",
  bursaphelenchus_xylophilus_prjea64437 => "Old Sanger genome",
  caenorhabditis_angaria_prjna51225 => "Caltech old caenorhabditis",
  caenorhabditis_elegans_prjna13758 => "Follow WormBase on C. elegans, not archives",
  fasciola_hepatica_prjna179522 => "WashU 50 helminths",
  haemonchus_contortus_prjna205202 => "Gasser lab alternative Haemonchus",
  oesophagostomum_dentatum_prjna72579 => "WashU 50 helminths",
  steinernema_carpocapsae_prjna202318 => "Caltech Steinernema",
  steinernema_feltiae_prjna204661 => "Caltech Steinernema",
  steinernema_glaseri_prjna204943 => "Caltech Steinernema",
  steinernema_monticolum_prjna205067 => "Caltech Steinernema",
  steinernema_scapterisci_prjna204942 =>  "Caltech Steinernema",
  taenia_solium_prjna170813 => "Mexico T. solium, an up-tinkered version of what we have",
  teladorsagia_circumcincta_prjna72569 => "WashU 50 helminths",
  trichinella_nativa_prjna179527 => "WashU 50 helminths",
);
for my $species (sort keys %assemblies) {
  BIOPROJECT:
  for my $bioproject (sort keys %{$assemblies{$species}}){
    my @ncbi_assemblies = @{$assemblies{$species}{$bioproject}{ncbi} // [] };
    my $core_db = $assemblies{$species}{$bioproject}{core_db};
    unless ($core_db) {
       print Dump { NCBI_ONLY => { $species => { $bioproject => \@ncbi_assemblies }}}; 
       next BIOPROJECT;
    }
    $assemblies{$species}{$bioproject}{core_db}{genebuild_start_date} = $production_mysql->meta_value($core_db->{name}, "genebuild.start_date") // "";
    $assemblies{$species}{$bioproject}{core_db}{provider_name} = $production_mysql->meta_value($core_db->{name}, "provider.name") // "";
    unless(@ncbi_assemblies){
       print Dump { PARASITE_ONLY => {$core_db->{name} => 1} };
       next BIOPROJECT;
    }
    if ($match_not_intended{"${species}_${bioproject}"}){
       print Dump { MATCH_NOT_INTENDED => {$core_db->{name} => $match_not_intended{"${species}_${bioproject}"} } };
       next BIOPROJECT;
    }
    $assemblies{$species}{$bioproject}{core_db}{assembly_accession} = $production_mysql->meta_value($core_db->{name}, "assembly.accession") // "";
    $assemblies{$species}{$bioproject}{core_db}{assembly_name} = $production_mysql->meta_value($core_db->{name}, "assembly.name") // "";

    my ($ncbi_assembly, @others) = grep {same_assembly_names($core_db, $_)} @ncbi_assemblies;
    croak @ncbi_assemblies, $core_db if @others;
    if ($ncbi_assembly) {
       my @other_assemblies = grep {not same_assembly_names($ncbi_assembly, $_ )} @ncbi_assemblies;
       delete $ncbi_assembly->{report_url};
       delete $_->{report_url} for @other_assemblies;
       print Dump { MATCHING_METADATA => { $species => { $bioproject => {%{$ncbi_assembly}, core_db_name => $core_db->{name}, other_assemblies => \@other_assemblies } }}};
       next BIOPROJECT;
    }
    my $scaffold_lengths_core = @ncbi_assemblies ? get_scaffold_lengths_for_core_db($core_db->{name}): {};
    my $scaffold_stats_core = scaffold_stats($scaffold_lengths_core);
    $core_db->{stats} = $scaffold_stats_core;
    for my $ncbi_assembly (@ncbi_assemblies){
       my @other_assemblies = grep {not same_assembly_names($ncbi_assembly, $_ )} @ncbi_assemblies;
       my $scaffold_lengths_ncbi = get_scaffold_lengths_from_report_url($ncbi_assembly->{report_url});
       if(Data::Compare::Compare ($scaffold_lengths_ncbi, $scaffold_lengths_core)){
           delete $core_db->{stats};
           delete $ncbi_assembly->{report_url};
           delete $_->{report_url} for @other_assemblies;
           print Dump { MATCHING_SCAFFOLD_NAMES => { $species => { $bioproject => {%{$ncbi_assembly}, core_db => $core_db, other_assemblies => \@other_assemblies } }}};
           next BIOPROJECT;
       } 
       my $scaffold_stats_ncbi = scaffold_stats($scaffold_lengths_ncbi);
       $ncbi_assembly->{stats} = $scaffold_stats_ncbi;
       if (Data::Compare::Compare([sort values %{$scaffold_lengths_ncbi}], [sort values %{$scaffold_lengths_core}])){
           delete $ncbi_assembly->{report_url};
           delete $_->{report_url} for @other_assemblies;
           delete $_->{stats} for @other_assemblies;
           print Dump { MATCHING_SCAFFOLD_LENGTHS  => { $species => { $bioproject => {
              %{$ncbi_assembly},
              core_db => $core_db,
              other_assemblies => \@other_assemblies,
            } }}};
           next BIOPROJECT;
       }
    }
    print Dump { MISMATCHED => { $species => { $bioproject => { %{$core_db}, ncbi => \@ncbi_assemblies} }}};
  }
}

sub scaffold_stats {
   my($h) = @_;
   my $longest_scaffold = "";
   my $num_scaffolds = 0;
   my $longest_scaffold_length = 0;
   while( my( $k, $v) = each %{$h}){
      $num_scaffolds++;
      if ($v > $longest_scaffold_length) {
         $longest_scaffold_length = $v;
         $longest_scaffold = $k;
      }
   }
   return {longest_scaffold =>  $longest_scaffold, num_scaffolds => $num_scaffolds, longest_scaffold_length => $longest_scaffold_length };
}

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

