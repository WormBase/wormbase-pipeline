
package GenomeBrowser::Resources::GeoMetadata;
use parent GenomeBrowser::Resources::LocallyCachedResource;
use Data::Dumper;
my $EUTILS_URL = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils';

sub _fetch {
    my ( $class, $species, $rnaseqer_metadata ) = @_;
    my %data;
    for my $assembly( @{$rnaseqer_metadata->access}){
      STUDY:
      for my $study_id ( @{ $rnaseqer_metadata->access($assembly) } ) {
         my ($web, $key) = &_session_bits_from_esearch_payload($class->get_xml(
            "$EUTILS_URL/esearch.fcgi?db=gds&term=${study_id}\[accn]&usehistory=y"
         ));
         next STUDY unless $web and $key;
         $data{$assembly}{$study_id} = &_data_from_esummary_payload($class->get_xml(
            "$EUTILS_URL/esummary.fcgi?db=gds&query_key=$key&WebEnv=$web"
         ));
      }
   }
   return \%data;
}

sub _session_bits_from_esearch_payload {
   my $payload = shift;
   return if $payload->{WarningList}{OutputMessage} eq 'No items found.';

   my $web = $payload->{WebEnv};
   my $key = $payload->{QueryKey};
   
   die Data::Dumper::Dumper( $payload) unless $web and $key;
   return $web, $key;
}

sub _data_from_esummary_payload {
   my $payload = shift;
   my @pubmed_ids =  map {$_->{content}||()} map {my $o = $_->{Item}; ref $o eq 'ARRAY' ? @$o : $o} grep ({$_->{Name} eq 'PubMedIds'} @{$payload->{DocSum}{Item}});
   return @pubmed_ids ? {pubmed => \@pubmed_ids}: {};
}
