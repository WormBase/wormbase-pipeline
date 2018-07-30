
package GenomeBrowser::RnaseqTracks;

use GenomeBrowser::ArrayExpressMetadata;
use GenomeBrowser::RnaseqerMetadata;
use GenomeBrowser::Studies;
use GenomeBrowser::Links;
use GenomeBrowser::Factors;
use GenomeBrowser::RnaseqerStats;
# TODO: Be able to create hub.txt from this
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
  my $links = GenomeBrowser::Links->new($root_dir, "${spe}_${cies}", $assembly, $rnaseqer_metadata);
  my $rnaseqer_stats = GenomeBrowser::RnaseqerStats->new($root_dir, "${spe}_${cies}", $assembly, $rnaseqer_metadata); 
  my $factors = GenomeBrowser::Factors->new($root_dir, "${spe}_${cies}", $assembly, $rnaseqer_metadata, $array_express_metadata);
  my @tracks;
  for my $study_id (@{$rnaseqer_metadata->access($assembly)}){
    my $study = $studies->{$study_id};
    for my $run_id (@{$rnaseqer_metadata->access($assembly, $study_id)}){
       my $stats = $rnaseqer_stats->get_formatted_stats($run_id);
       my $links_to_this_run = $links->{$run_id};
       my %attributes;
       for my $characteristic_type (@{$rnaseqer_metadata->access($assembly, $study_id, $run_id)}){
         $attributes{$characteristic_type} = $rnaseqer_metadata->access($assembly, $study_id, $run_id, $characteristic_type);
       }
       push @tracks, {
          run_id => $run_id,
          attributes => {%$study, %$stats, %$links_to_this_run, %attributes},
          label => 
            &_label($run_id, &_restrict_sample_name($cies, $attributes{sample_name}),  map {exists $attributes{$_} ? $attributes{$_}: ()} @$factors),
       };
    }
  }
  return $factors, $rnaseqer_metadata->{location_per_run_id}, @tracks;
}
sub _restrict_sample_name {
  my ($cies, $sample_name) = @_;
  return "" unless $sample_name;
  return "" if $sample_name =~ /$cies/i and $sample_name =~ /sample from/i;
  return $sample_name;
}

#"Sample Name" is a special attribute in ENA, and therefore usually present
#Sometimes it's not informative enough
sub _sample_name_seems_good_enough {
  my $sample_name = shift;
  my @ws = split /\W+/, $sample_name;
  return ($sample_name and length ($sample_name) > 10 and @ws > 1);
}
sub _label {
  my ($run_id,$sample_name, @factor_values) = @_;

  my $description;
  if(_sample_name_seems_good_enough($sample_name)){
    @factor_values = grep (!/^[A-Z]+[0-9]+$/, @factor_values);
    $description = $sample_name;
  } else {
    @factor_values = grep (/\w+/, @factor_values);
    $description = join(", ", @factor_values);
  }
  return $description ? "$run_id: $description" : $run_id;
}
1;
