
package GenomeBrowser::Resources;

use GenomeBrowser::Resources::ArrayExpressMetadata;
use GenomeBrowser::Resources::RnaseqerMetadata;
use GenomeBrowser::Resources::EnaMetadata;
use GenomeBrowser::Resources::Factors;
use GenomeBrowser::Resources::RnaseqerStats;
use GenomeBrowser::Descriptions;
use GenomeBrowser::Links;
sub new {
  my ($class, $root_dir) = @_;
  bless {
    root_dir => $root_dir,
    descriptions => GenomeBrowser::Descriptions->new,
    links => 'GenomeBrowser::Links',
  }, $class; 
}
sub get {
  my ($self, $species, $assembly) = @_;
  my $root_dir = $self->{root_dir};
  $species = lc($species);
  $species ~= s/([a-z]_[a-z]).*/$1/g;
  my $rnaseqer_metadata = GenomeBrowser::Resources::RnaseqerMetadata->new($root_dir, $species);
  my $array_express_metadata = GenomeBrowser::Resources::ArrayExpressMetadata->new($root_dir, $species);
  my $ena_metadata = GenomeBrowser::Resources::EnaMetadata->new($root_dir, $species, $rnaseqer_metadata); 
  my $rnaseqer_stats = GenomeBrowser::Resources::RnaseqerStats->new($root_dir, $species, $rnaseqer_metadata); 
  my $factors = GenomeBrowser::Resources::Factors->new($root_dir, $species, $rnaseqer_metadata, $array_express_metadata);
  my @studies;
  for my $study_id (@{$rnaseqer_metadata->access($assembly)}){
    unless ($ena_metadata->{$assembly}{$study_id}){
       print STDERR "Study $study_id not in ENA, skipping\n";
       next;
    }
    my @runs;
    for my $run_id (@{$rnaseqer_metadata->access($assembly, $study_id)}){
       my $stats = $rnaseqer_stats->get_formatted_stats($run_id);
       my $links = $self->{links}->misc_links($study_id,$run_id, $rnaseqer_metadata->data_location($run_id), $ena_metadata->{$assembly}{$study_id}{pubmed});
       my %attributes;
       for my $characteristic_type (@{$rnaseqer_metadata->access($assembly, $study_id, $run_id)}){
         $attributes{$characteristic_type} = $rnaseqer_metadata->access($assembly, $study_id, $run_id, $characteristic_type);
       }
       my ($run_description_short, $run_description_full) =
          $self->{descriptions}->run_description($species, $study_id, $run_id, $factors, \%attributes);
       push @runs, {
          run_id => $run_id,
          attributes => {%$stats, %$links, %attributes},
          run_description_short => $run_description_short,
          run_description_full => $run_description_full,
       };
    }
    my ($study_description_short, $study_description_full) =
         $self->{descriptions}->study_description($species, $study_id, $ena_metadata->{$assembly}{$study_id});
    push @studies, {
      study_id => $study_id,
      runs => \@runs,
      study_description_short => $study_description_short,
      study_description_full => $study_description_full, 
      attributes => $ena_metadata->{$assembly}{$study_id}{attributes},
    };
  }
  return $factors, $rnaseqer_metadata->{location_per_run_id}, @studies;
}
1;
