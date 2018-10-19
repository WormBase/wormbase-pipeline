
package GenomeBrowser::Resources;

use GenomeBrowser::Resources::ArrayExpressMetadata;
use GenomeBrowser::Resources::RnaseqerMetadata;
use GenomeBrowser::Resources::StudyAttributes;
use GenomeBrowser::Resources::Links;
use GenomeBrowser::Resources::Factors;
use GenomeBrowser::Resources::RnaseqerStats;
use GenomeBrowser::Descriptions;
sub new {
  my ($class, $root_dir) = @_;
  bless {root_dir => $root_dir, descriptions => GenomeBrowser::Descriptions->new}, $class; 
}
sub get {
  my ($self, $species, $assembly) = @_;
  my $root_dir = $self->{root_dir};
  $species = lc($species);
  $species ~= s/([a-z]_[a-z]).*/$1/g;
  my $rnaseqer_metadata = GenomeBrowser::Resources::RnaseqerMetadata->new($root_dir, $species);
  my $array_express_metadata = GenomeBrowser::Resources::ArrayExpressMetadata->new($root_dir, $species);
  my $study_attributes = GenomeBrowser::Resources::StudyAttributes->new($root_dir, $species, $rnaseqer_metadata); 
  my $links = GenomeBrowser::Resources::Links->new($root_dir, $species, $rnaseqer_metadata);
  my $rnaseqer_stats = GenomeBrowser::Resources::RnaseqerStats->new($root_dir, $species, $rnaseqer_metadata); 
  my $factors = GenomeBrowser::Resources::Factors->new($root_dir, $species, $rnaseqer_metadata, $array_express_metadata);
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
       my ($run_description_short, $run_description_full) = 
          $self->{descriptions}->run_description($species, $study_id, $run_id, $factors, \%attributes);
       push @runs, {
          run_id => $run_id,
          attributes => {%$stats, %$links_to_this_run, %attributes},
          run_description_short => $run_description_short,
          run_description_full => $run_description_full 
       };
    }
    my ($study_description_short, $study_description_full) =
         $self->{descriptions}->study_description($species, $study_id, $study_attributes->{$assembly}{$study_id});
    push @studies, {
      study_id => $study_id,
      runs => \@runs,
      study_description_short => $study_description_short,
      study_description_full => $study_description_full, 
      attributes => $study_attributes->{$assembly}{$study_id}{properties},
    };
  }
  return $factors, $rnaseqer_metadata->{location_per_run_id}, @studies;
}
1;
