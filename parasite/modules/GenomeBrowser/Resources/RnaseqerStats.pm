
package GenomeBrowser::Resources::RnaseqerStats;
use parent GenomeBrowser::Resources::LocallyCachedResource;


sub _fetch {
    my ( $class, $species, $rnaseqer_metadata ) = @_;
    my %data;
    for my $assembly( @{$rnaseqer_metadata->access}){
      for my $study_id ( @{ $rnaseqer_metadata->access($assembly) } ) {
        for my $run_id (@{ $rnaseqer_metadata->access($assembly, $study_id) } ) {
          my $location = $rnaseqer_metadata->data_location($run_id); 
          my $stats = &_get_pairs(
              $class->get_csv("$location/$run_id.se.hits.bam.stats.csv", "$location/$run_id.pe.hits.bam.stats.csv")
          );
          my $a = $stats->{All_entries};
          my $u = $stats->{UniquelyMappedReads};
          $data{$run_id}{library_size} = $a || 0 ;
          $data{$run_id}{fraction_reads_uniquely_mapped} = $a ? sprintf("%.3f", ($u || 0) / $a) : 0;
        }
      }
    }
    return \%data;
}
sub round_library_size {
   my $n = shift;
   return "?" unless $n;
   return "75M+" if $n   > 75000000;
   return "50-75M" if $n > 50000000;
   return "25-50M" if $n > 25000000;
   return  "5-25M" if $n > 5000000;
   return   "1-5M" if $n > 1000000;
   return "200k-1M" if $n >200000;
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
    library_size_reads=> $self->{$run_id}{library_size},
    library_size_reads_approximate => round_library_size($self->{$run_id}{library_size}),
    fraction_of_reads_mapping_uniquely_approximate => round_fraction_reads_uniquely_mapped($self->{$run_id}{fraction_reads_uniquely_mapped}),
    fraction_of_reads_mapping_uniquely=> $self->{$run_id}{fraction_reads_uniquely_mapped}
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
