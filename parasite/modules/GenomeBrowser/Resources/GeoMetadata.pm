
package GenomeBrowser::Resources::GeoMetadata;
use parent GenomeBrowser::Resources::LocallyCachedResource;

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
   use Test::More;
   die explain @_;
   return "web", "key";
}

sub _data_from_esummary_payload {
   use Test::More;
   die explain @_;
   return {pubmed => []};
}
