
package GenomeBrowser::Factors;
use List::Util qw(sum);
use List::MoreUtils qw(uniq);
use parent GenomeBrowser::LocallyCachedResource;

# Factors for species are a union across studies of:
# - either what ArrayExpress says the factor types are for the study
# - or RNASeq-er characteristic types that vary across runs for the study
# Compare in ArrayExpress: a factor is an "important" sample characteristic

my @blacklist = (
 "synonym",
 "bioproject_id",
  "species",
  "organism",
  "replicate"
);
sub not_in_blacklist {
  my $item = shift;
  $item = lc($item);
  $item =~ s/\W+/_/;
  return not($item ~~ @blacklist);
}
sub _fetch {
  my ( $class, $species, $assembly, $rnaseqer_metadata,
    $array_express_metadata ) = @_;
  my %data;
  for my $study_id ( @{ $rnaseqer_metadata->access($assembly) } ) {
    my %rnaseqer_characteristics;
    my @runs = @{ $rnaseqer_metadata->access( $assembly, $study_id ) };
    for my $run (@runs) {
      for my $characteristic_type (
        @{ $rnaseqer_metadata->access( $assembly, $study_id, $run ) } )
      {
        $rnaseqer_characteristics{$characteristic_type}{
          $rnaseqer_metadata->access( $assembly, $study_id, $run,
            $characteristic_type )
        }++;
      }
    }
    my @factors;
    my @ae_factors =
      @{ $array_express_metadata->factor_types($study_id) // [] };
    if (@ae_factors) {
      @factors =
        grep { exists $rnaseqer_characteristics{$_} } @ae_factors;
    }
    else {
      for ( keys %rnaseqer_characteristics ) {
        my %d      = %{ $rnaseqer_characteristics{$_} };
        my @values = keys %d;
        my @counts = values %d;
        push @factors, $_ if @values > 1 or sum(@counts) == 1;
      }
    }
    $data{$study_id} = \@factors;
  }
  my @result = map { @{$_} } ( values %data );
  @result = grep { not_in_blacklist($_) } @result;
  @result = uniq @result;
  @result = sort @result;
  return \@result;
}


1;