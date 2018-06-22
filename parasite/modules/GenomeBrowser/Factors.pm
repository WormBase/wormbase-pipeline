
package GenomeBrowser::Factors;

use List::MoreUtils qw(uniq);
use parent GenomeBrowser::LocallyCachedResource;

# Factors for species are a union across studies of:
# - either what ArrayExpress says the factor types are for the study
# - or RNASeq-er characteristic types that vary across runs for the study
# Compare in ArrayExpress: a factor is an "important" sample characteristic

sub _fetch {
  my ($class, $species, $assembly, $rnaseqer_metadata, $array_express_metadata) = @_;
  my %data;
  for my $study_id ($rnaseqer_metadata->access($assembly)){
    my %rnaseqer_characteristics;
    for my $run ($rnaseqer_metadata->access($assembly, $study_id)){
       for my $characteristic_type ($rnaseqer_metadata->access($assembly, $study_id, $run)){
         $rnaseqer_characteristics{$characteristic_type}++;
       } 
    }
    my @factors = @{$array_express_metadata->factor_types($study_id) // []};
    @factors = keys %rnaseqer_characteristics unless @factors;
    @factors = map {$_ if $rnaseqer_characteristics{$_} > 1 } @factors;
    $data{$study_id}=\@factors;
  }
  return \{sort (uniq( map {@{$_}} (values %data)))};
}
