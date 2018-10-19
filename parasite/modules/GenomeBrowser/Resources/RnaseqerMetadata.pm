
package GenomeBrowser::Resources::RnaseqerMetadata;
use List::MoreUtils qw(uniq);
use File::Basename;
use parent GenomeBrowser::Resources::LocallyCachedResource;

# Assembly -> study -> run -> type -> value

sub access {
  my $self = shift;
  my $h = $self->{metadata};
  while(@_ and ref $h){
    $h = $h->{shift @_};
  }
  return [sort keys %$h ] if ref $h;
  return [] unless $h;
  return $h;
}

sub data_location {
  my ($self, $run_id) = @_; 
  my $bigwig_location = $self->{location_per_run_id}{$run_id};
  return dirname $bigwig_location; 
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
       my ($type, $value) = _normalise_type_and_value($attribute_record->{TYPE}, $attribute_record->{VALUE});
       $run_attributes{$attribute_record->{RUN_ID}}{$type} = $value if $type and $value;
     }
     $run_attributes{$attribute_record->{RUN_ID}} = _normalise_characteristics_hash($run_attributes{$attribute_record->{RUN_ID}});
  }
  my %data;
  for my $run_record (@$run_records){
      $data
        {$run_record->{ASSEMBLY_USED}}
        {$run_record->{STUDY_ID}}
        {$run_record->{RUN_IDS}}
        = $run_attributes{$run_record->{RUN_IDS}}; 
  }

  my %location_per_run_id;
  for my $run_record (@$run_records){
  (my $bigwig_we_want = $run_record->{BIGWIG_LOCATION}) =~ s/.bw$/.nospliced.bw/;
      $location_per_run_id
         {$run_record->{RUN_IDS}}
         = $bigwig_we_want;
  }
  return {metadata => \%data, location_per_run_id => \%location_per_run_id};
}
sub _normalise_type_and_value {
  my ($type, $value) = @_;

  $type = lc($type);
  $type =~ s/\W+/_/g;

#Sometimes there's curation like: age+time unit
  return $type, $value if $type eq "age" and $value =~s/^\W+$//;
# But sometimes age means developmental stage
  $type =~ s/^age$/developmental_stage/;

# Try give developmental stages first-class support
  $type =~ s/^stage$/developmental_stage/;
  $type =~ s/life_cycle_stage/developmental_stage/;
  $type =~ s/dev_stage/developmental_stage/;
  $type =~ s/development_stage/developmental_stage/;
  if( $type eq "developmental_stage"){
     $value =~ s/^ADULT$/adult/;
     $value =~ s/^Adult$/adult/;
     $value =~ s/^adult worms?$/adult/;
  }
  $type =~ s/^tissue$/organism_part/;
  if($type eq "organism_part" ){
     $value = lc($value);
#http://purl.obolibrary.org/obo/UBERON_0000468
     $value =~ s/^whole$/whole organism/;
     $value =~ s/^whole body$/whole organism/;
     $value =~ s/^whole ?-?_?worms? ?-?_?(tissue)?$/whole organism/;
     $value =~ s/^whole ?-?_?organisi?m$/whole organism/;
     $value =~ s/^whole ?-?_?animals?$/whole organism/;
     $value =~ s/^intact ?-?_?animals?$/whole organism/;
  }
  if($type eq "host_disease"){
     $value = lc($value);
     $value =~ s/^schistosomiase$/schistosomiasis/;
  }
  $type =~ s/time_point/timepoint/;

  $value =~s/^not applicable.*$//i;
  $value =~s/^unknown$//i;
  $value =~s/^N\\?A$//i;
  $value =~s/^\W+$//;

  return $type, $value;
}

# Pretends to be referentially transparent
sub _normalise_characteristics_hash {
  my $o = shift;
  if(not $o->{sex} and $o->{developmental_stage} =~ /^adult ?-?_?males?$/){
     $o->{sex} = "male";
     $o->{developmental_stage} = "adult";
  }
  if(not $o->{sex} and $o->{developmental_stage} =~ /^adult ?-?_?females?$/){
     $o->{sex} = "female";
     $o->{developmental_stage} = "adult";
  }
  if(not $o->{sex} and $o->{organism_part} =~ /^whole ?-?_?males?$/){
     $o->{sex} = "male";
     $o->{organism_part} = "whole organism";
  }
  if(not $o->{sex} and $o->{organism_part} =~ /^whole ?-?_?females?$/){
     $o->{sex} = "female";
     $o->{organism_part} = "whole organism";
  }
  return $o;
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
  my $payload = $class->_get_rnaseqer_json(
    "30", "getRunsByOrganism", $species
  );
  my @a = grep {$_->{BIGWIG_LOCATION} ne "NA"} @$payload;
  return \@a;
}
sub _get_rnaseqer_sample_attributes_per_run_for_study {
  my ($class, $study) = @_;
  return $class->_get_rnaseqer_json(
    "getSampleAttributesPerRunByStudy", $study
  );
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
