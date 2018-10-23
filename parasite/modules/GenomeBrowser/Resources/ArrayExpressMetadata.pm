
package GenomeBrowser::Resources::ArrayExpressMetadata;

use parent GenomeBrowser::Resources::LocallyCachedResource;

sub _val {
   my ($self, $study_accession, $val) = @_;
   my $key = $self->{secondary_to_primary_accession}->{$study_accession};
   $key //= $study_accession;
   return $self->{$val}{$key};
}

sub factor_types { return &_val(@_, "primary_accession_to_factor_type");}
sub pubmed { return &_val(@_, "primary_accession_to_pubmed");}

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
    primary_accession_to_pubmed => _get_pubmed_from_payload($payload),
    primary_accession_to_factor_type => &_get_factor_types_from_payload($payload),
    secondary_to_primary_accession => &_get_secondary_to_primary_accession_from_payload($payload),
  };
}
# We use these values to pick up factors, so make sure they're the same format!
sub _normalise_characteristic {
  my $type = shift;
  $type = lc($type);
  $type =~ s/\W+/_/g;
  return $type;
}
sub _get_factor_types_from_payload{
  my $payload = shift;
  my %result;
  for my $experiment (@{$payload->{experiments}->{experiment}}){
    my @factor_types = map { _normalise_characteristic($_->{name})} @{ $experiment->{experimentalvariable}};
    $result{$experiment->{accession}} = \@factor_types;
  }
  return \%result;
}

sub _get_secondary_to_primary_accession_from_payload{
  my $payload = shift;
  my %result; 
  for my $experiment (@{$payload->{experiments}->{experiment}}){
    for my $secondary_accession (@{$experiment->{secondaryaccession}}){
       $result{$secondary_accession}=$experiment->{accession}
    }
  }
  return \%result;
}
sub _get_pubmed_from_payload {
  my $payload = shift;
  my %result;
  for my $experiment (@{$payload->{experiments}{experiment}}){
     next unless $experiment->{bibliography};
     $result{$experiment->{accession}} = [map {$_->{accession} || ()} @{$experiment->{bibliography}}];
  }
  return \%result;
}
1;
