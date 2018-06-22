#!/usr/bin/perl

use LWP::UserAgent;
use JSON;
use ProductionMysql;
use YAML;
my $core_db = shift or die "Usage: $0 <core_db>";

my ($spe, $cies, $bioproject) = split "_", $core_db;

my $assembly = ProductionMysql->staging->meta_value($core_db, "assembly.name" );

my %result;
$result{species} = join "_", $spe, $cies, $bioproject;
$result{assembly} = $assembly;

for my $study (@{&rnaseqerGetStudiesByOrganism("${spe}_${cies}", $assembly)}) {
  my %study;
  
  my $runs = &rnaseqerGetSampleAttributesPerRunByStudy($study);
  my %runs;
  for my $run (keys %$runs) {
    $runs{$run}{attributes} = $runs{$run};
    $runs{$run}{bigwig} = &find_path_to_rnaseqer_alignment_file("${spe}_${cies}", $run);
  }
 
  $study{title}="TODO fetch from ENA with Bruce's code"; 

  $result{studies}{$study}=\%study;
}
print Dump(\%result);


sub get_json_from_rnaseqer {
  my $response = LWP::UserAgent->new()->get(
    join ("/", "https://www.ebi.ac.uk/fg/rnaseq/api/json", @_ )
  );
  die "Rnaseqer @_ error : " .$response->status_line unless $response->is_success;
  return from_json($response->decoded_content);
}

sub rnaseqerGetStudiesByOrganism {
  my ($species, $assembly) = @_;
  my @result;
  for my $study (@{&get_json_from_rnaseqer("getStudiesByOrganism", $species)}){
    push @result, $study->study_id if $study->assembly_used eq $assembly;
  }
  return @result;
}

sub rnaseqerGetSampleAttributesPerRunByStudy {
  return &get_json_from_rnaseqer("getSampleAttributesPerRunByStudy", @_);
}
sub find_path_to_rnaseqer_alignment_file {

("${spe}_${cies}", $run);
