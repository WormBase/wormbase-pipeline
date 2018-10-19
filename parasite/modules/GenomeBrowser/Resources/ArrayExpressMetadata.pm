
package GenomeBrowser::Resources::ArrayExpressMetadata;

use parent GenomeBrowser::Resources::LocallyCachedResource;

sub factor_types {
   my ($self,$study_accession) = @_;

   my $key = $self->{secondary_to_primary_accession}->{$study_accession};
   $key //= $study_accession;
   return $self->{primary_accession_to_factor_type}->{$key};

}
sub _fetch {
  my ($class, $species) = @_;
  return &_create_from_payload(
    $class->get_json( &_url($species))
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
# We use these values to pick up factors, so make sure they're the same format!
sub _normalise_characteristic {
  my $type = shift;
  $type = lc($type);
  $type =~ s/\W+/_/g;
  return $type;
}
sub _get_factor_types_from_ae_response{
  my $payload = shift;
  my %result;
  for my $experiment (@{$payload->{experiments}->{experiment}}){
    my @factor_types = map { _normalise_characteristic($_->{name})} @{ $experiment->{experimentalvariable}};
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
