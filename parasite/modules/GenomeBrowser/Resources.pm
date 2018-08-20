
package GenomeBrowser::Resources;

use GenomeBrowser::ArrayExpressMetadata;
use GenomeBrowser::RnaseqerMetadata;
use GenomeBrowser::StudyAttributes;
use GenomeBrowser::Links;
use GenomeBrowser::Factors;
use GenomeBrowser::RnaseqerStats;
use GenomeBrowser::Descriptions;
sub new {
  my ($class, $root_dir) = @_;
  bless {root_dir => $root_dir}, $class; 
}
sub get {
  my ($self, $species, $good_assembly) = @_;
  my $root_dir = $self->{root_dir};
  my ($spe, $cies, $bioproject) = split "_", $species;
  my $rnaseqer_metadata = GenomeBrowser::RnaseqerMetadata->new($root_dir, "${spe}_${cies}");
  my $array_express_metadata = GenomeBrowser::ArrayExpressMetadata->new($root_dir, "${spe}_${cies}");
  my $study_attributes = GenomeBrowser::StudyAttributes->new($root_dir, "${spe}_${cies}", $rnaseqer_metadata); 
  my $links = GenomeBrowser::Links->new($root_dir, "${spe}_${cies}", $rnaseqer_metadata);
  my $rnaseqer_stats = GenomeBrowser::RnaseqerStats->new($root_dir, "${spe}_${cies}", $rnaseqer_metadata); 
  my $factors = GenomeBrowser::Factors->new($root_dir, "${spe}_${cies}", $rnaseqer_metadata, $array_express_metadata);
  my $descriptions = GenomeBrowser::Descriptions->new("${spe}_${cies}");
  my @studies;
  return unless grep(/$good_assembly/, @{$rnaseqer_metadata->access}); # Patch for ISL bug. TODO remove!
  for my $assembly (($good_assembly , '')){ # Patch for ISL bug. TODO remove!
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
          $descriptions->run_description($study_id, $run_id, $factors, \%attributes);
       push @runs, {
          run_id => $run_id,
          attributes => {%$stats, %$links_to_this_run, %attributes},
          run_description_short => $run_description_short,
          run_description_full => $run_description_full 
       };
    }
    push @studies, {
      study_id => $study_id,
      runs => \@runs,
      study_description => $descriptions->study_description($study_id, $study_attributes->{$assembly}{$study_id}),
      attributes => $study_attributes->{$assembly}{$study_id},
    };
  }
} #ENDPATCH
  return $factors, $rnaseqer_metadata->{location_per_run_id}, @studies;
}
1;
