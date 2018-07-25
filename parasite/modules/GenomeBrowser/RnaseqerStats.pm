
package GenomeBrowser::RnaseqerStats;
use File::Basename;
use parent GenomeBrowser::LocallyCachedResource;


sub _fetch {
    my ( $class, $species, $assembly, $rnaseqer_metadata ) = @_;
    my %data;
    for my $study_id ( @{ $rnaseqer_metadata->access($assembly) } ) {
      for my $run_id (@{ $rnaseqer_metadata->access($assembly, $study_id) } ) {
        my $bigwig_location = $rnaseqer_metadata->{location_per_run_id}{$run_id};
        my $location = dirname $bigwig_location; 
        my $stats = &_get_pairs(
            $class->get_csv("$location/$run_id.pe.hits.bam.stats.csv")
        );
        my $a = $stats->{All_entries};
        my $u = $stats->{UniquelyMappedReads};
        $data{$run_id}{library_size} = $a || 0 ;
        $data{$run_id}{fraction_reads_uniquely_mapped} = $a ? ($u || 0) / $a : 0;
      }
    }
    return \%data;
}
sub get_formatted_stats {
  my ($self, $run_id) = @_;
  return $self->{$run_id};
}

sub _get_pairs {
  my ($lines) = @_;
  for (@$lines){
    my ($k, $number, @others) = @$_;
    die "Could not parse:", @$_ if @others;
    $result{$k} = 0+$number;
  }
  return \%result;
}
1;
