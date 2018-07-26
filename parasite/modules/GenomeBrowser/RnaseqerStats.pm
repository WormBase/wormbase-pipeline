
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
            $class->get_csv("$location/$run_id.pe.hits.bam.stats.csv", "$location/$run_id.se.hits.bam.stats.csv")
        );
        my $a = $stats->{All_entries};
        my $u = $stats->{UniquelyMappedReads};
        $data{$run_id}{library_size} = $a || 0 ;
        $data{$run_id}{fraction_reads_uniquely_mapped} = $a ? ($u || 0) / $a : 0;
      }
    }
    return \%data;
}
sub round_library_size {
   my $n = shift;
   return "?" unless $n;
   return "75mln+" if $n   > 75000000;
   return "50-75mln" if $n > 50000000;
   return "25-50mln" if $n > 25000000;
   return  "5-25mln" if $n > 5000000;
   return   "1-5mln" if $n > 1000000;
   return "200k-1mln" if $n >200000;
   return "under 200k"; 
}

sub round_fraction_reads_uniquely_mapped {
  my $f = shift;
  return "?" unless $f;
  return "90+%" if $f > 0.9;
  return "80-90%" if $f > 0.8;
  return "70-80%" if $f > 0.7;
  return "60-70%" if $f > 0.6;
  return "50-60%" if $f > 0.5;
  return "under 50%";
}

sub get_formatted_stats {
  my ($self, $run_id) = @_;
  return {
    library_total_amount_of_reads=>$self->{$run_id}{library_size},
    library_size_approx => round_library_size($self->{$run_id}{library_size}),
    mapping_fraction_of_uniquely_mapped_reads=>$self->{$run_id}{fraction_reads_uniquely_mapped},
    mapping_quality_approx => round_fraction_reads_uniquely_mapped($self->{$run_id}{fraction_reads_uniquely_mapped}),
  };
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
