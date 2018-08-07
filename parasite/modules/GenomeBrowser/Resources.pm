
package GenomeBrowser::Resources;

use GenomeBrowser::ArrayExpressMetadata;
use GenomeBrowser::RnaseqerMetadata;
use GenomeBrowser::StudyAttributes;
use GenomeBrowser::Links;
use GenomeBrowser::Factors;
use GenomeBrowser::RnaseqerStats;
# TODO: Be able to create hub.txt from this
# Structure the returned data as study->run
# The client code can squash the study and run code if it likes to
# Let them have short/long label
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
  my $study_attributes = GenomeBrowser::StudyAttributes->new($root_dir, "${spe}_${cies}", $rnaseqer_metadata); 
  my $links = GenomeBrowser::Links->new($root_dir, "${spe}_${cies}", $rnaseqer_metadata);
  my $rnaseqer_stats = GenomeBrowser::RnaseqerStats->new($root_dir, "${spe}_${cies}", $rnaseqer_metadata); 
  my $factors = GenomeBrowser::Factors->new($root_dir, "${spe}_${cies}", $rnaseqer_metadata, $array_express_metadata);
  my @studies;
  for my $study_id (@{$rnaseqer_metadata->access($assembly)}){
    unless ($study_attributes->{$assembly}{$study_id}){
       print STDERR "Study $study_id not in ENA, skipping\n";
       next;
    }
    my @runs;
    for my $run_id (@{$rnaseqer_metadata->access($assembly, $study_id)}){
       my $stats = $rnaseqer_stats->get_formatted_stats($run_id);
       my $links_to_this_run = $links->{$run_id};
       my %attributes;
       for my $characteristic_type (@{$rnaseqer_metadata->access($assembly, $study_id, $run_id)}){
         $attributes{$characteristic_type} = $rnaseqer_metadata->access($assembly, $study_id, $run_id, $characteristic_type);
       }
       push @runs, {
          run_id => $run_id,
          attributes => {%$stats, %$links_to_this_run, %attributes},
          run_description => 
            &_run_description($run_id, &_restrict_sample_name($cies, $attributes{sample_name}),  map {exists $attributes{$_} ? $attributes{$_}: ()} @$factors),
       };
    }
    push @studies, {
      study_id => $study_id,
      runs => \@runs,
      study_description => &_study_description(%$study_attributes),
      attributes => $study_attributes->{$assembly}{$study_id},
    };
  }
  return $factors, $rnaseqer_metadata->{location_per_run_id}, @studies;
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
sub _run_description {
  my ($run_id,$sample_name, @factor_values) = @_;

  my $description;
  if(&_sample_name_seems_good_enough($sample_name)){
    @factor_values = grep (!/^[A-Z]+[0-9]+$/, @factor_values);
    $description = $sample_name;
  } else {
    @factor_values = grep (/\w+/, @factor_values);
    $description = join(", ", @factor_values);
  }
  return $description ? "$run_id: $description" : $run_id;
}
sub _study_description {
  my (%study_attributes) = @_;
  return ( $study_attributes{"Study description"} || $study_attributes{study} || $study_attributes{study_id}); 
}
1;
