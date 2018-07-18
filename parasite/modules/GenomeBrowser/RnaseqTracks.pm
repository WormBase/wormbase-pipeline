
package GenomeBrowser::RnaseqTracks;

use GenomeBrowser::ArrayExpressMetadata;
use GenomeBrowser::RnaseqerMetadata;
use GenomeBrowser::Studies;
use GenomeBrowser::Factors;
# This represents track metadata
# Creating the object fetches info and saves it on disk
# You can then hand-edit the files, especially the factors file
# TODO: Be able to create hub.txt from this

# in: 
#   - data production directory where we can cache resources
#   - current core database

sub new {
  my ($class, $root_dir) = @_;
  bless {root_dir => $root_dir}, $class; 
}
sub get {
  my ($self, $species, $assembly) = @_;
  my $root_dir = $self->{root_dir};
  my ($spe, $cies, $bioproject) = split "_", $species;
  my $rnaseqer_metadata = GenomeBrowser::RnaseqerMetadata->new($root_dir, "${spe}_${cies}");
  my $array_express_metadata = GenomeBrowser::ArrayExpressMetadata->new($root_dir, "${spe}_${cies}");
  my $studies = GenomeBrowser::Studies->new($root_dir, "${spe}_${cies}", $assembly, $rnaseqer_metadata); 
  my $factors = GenomeBrowser::Factors->new($root_dir, "${spe}_${cies}", $assembly, $rnaseqer_metadata, $array_express_metadata);
  my @tracks;
  for my $study_id (@{$rnaseqer_metadata->access($assembly)}){
    my $study_attributes = $studies->{$study_id};
    for my $run_id (@{$rnaseqer_metadata->access($assembly, $study_id)}){
       my $run_attributes = $rnaseqer_metadata->{$assembly}{$study_id}{$run_id};
       my %attributes = %{{ %$study_attributes, %$run_attributes }};

       push @tracks, {
          run_id => $run_id,
          attributes => \%attributes,
          label => 
            &_label($run_id, map {exists $attributes{$_} ? $attributes{$_}: ()} @$factors),
       };
    }
  }
  return $factors, @tracks;
}
# You can move this out if convenient!
sub _label {
  my ($run_id, @factor_values) = @_;
  return $run_id unless @factor_values;
  return "$run_id: ". join (", ",
    grep !/[A-Z]+[0-9]+/, @factor_values
  ); 
}
1;
