
package GenomeBrowser::RnaseqerMetadata;
use List::MoreUtils qw(uniq);
use parent GenomeBrowser::LocallyCachedResource;

# Assembly -> study -> run -> type -> value

sub access {
  my $h = shift;
  while(@_ and ref $h){
    $h = $h->{shift @_};
  }
  return [sort keys %$h ] if ref $h;
  return [] unless $h;
  return $h;
}
sub _fetch {
  my ($class, $species) = @_;
  $species= lc($species);
  $species =~ s/[^a-zA-Z0-9]+/_/g;

  my $run_records = $class->_get_rnaseqer_runs_for_organism($species);

  my %run_attributes;
  for my $study_id (uniq (map {$_->{STUDY_ID}} @$run_records)){
     for my $attribute_record (@{ 
       $class->_get_rnaseqer_sample_attributes_per_run_for_study($study_id)
     }){
       $run_attributes
         {$attribute_record->{RUN_ID}}
         {$attribute_record->{TYPE}} 
         = $attribute_record->{VALUE};
     }
  }
  my %data;
  for my $run_record (@{$class->_get_rnaseqer_runs_for_organism($species)}){
      $data
        {$run_record->{ASSEMBLY_USED}}
        {$run_record->{STUDY_ID}}
        {$run_record->{RUN_IDS}}
        = $run_attributes{$run_record->{RUN_IDS}}; 
  }
  return \%data;
}
sub _get_rnaseqer_json {
  my $class = shift;
  return $class->get_json(
    join ("/",
      "https://www.ebi.ac.uk/fg/rnaseq/api/json",
      @_
    )
  );
}
# Ask for runs that had at least 30% reads mapped
# We want to exclude failures and queued entries
# RNASeq-er doesn't have complete records anyway!
sub _get_rnaseqer_runs_for_organism {
  my ($class, $species) = @_;
  return $class->_get_rnaseqer_json(
    "30", "getRunsByOrganism", $species
  );
}
sub _get_rnaseqer_sample_attributes_per_run_for_study {
  my ($class, $study) = @_;
  return $class->_get_rnaseqer_json(
    "getSampleAttributesPerRunByStudy", $study
  );
}
sub _url {
  $species = shift;
  $species= lc($species);
  $species =~ s/[^a-zA-Z0-9]+/%20/g;
  return "https://www.ebi.ac.uk/arrayexpress/json/v3/experiments?species=$species";
}

sub _create_from_payload {
  $payload = shift;
  return {
    primary_accession_to_factor_type => &_get_factor_types_from_ae_response($payload),
    secondary_to_primary_accession => &_get_secondary_to_primary_accession_from_ae_response($payload),
  };
}

sub _get_factor_types_from_ae_response{
  my $payload = shift;
  my %result;
  for my $experiment (@{$payload->{experiments}->{experiment}}){
    my @factor_types = map { $_->{name}} @{ $experiment->{experimentalvariable}};
    $result{$experiment->{accession}} = \@factor_types;
  }
  return \%result;
}

sub _get_secondary_to_primary_accession_from_ae_response{
  my $payload = shift;
  my %result; 
  for my $experiment (@{$payload->{experiments}->{experiment}}){
    for my $secondary_accession (@{$experiment->{secondaryaccession}}){
       $result{$secondary_accession}=$experiment->{accession}
    }
  }
  return \%result;
}

1;
