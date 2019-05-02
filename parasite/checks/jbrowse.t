use strict;
use warnings;
use JSON;
use LWP;
use URI;
use Test::More;
use ProductionMysql;

test_tracks_for_species("http://test.parasite.wormbase.org", $_)
  for ProductionMysql->staging->species(@ARGV);
done_testing;

sub test_tracks_for_species {
  my ($url_base, $species) = @_;
  my $track_list_url = join("/", $url_base, "jbrowse-data", $species, "data", "trackList.json");
  my $response = LWP::UserAgent->new->get($track_list_url);
  ok ($response->is_success, "$species config live: $track_list_url") or return;
  my $track_list = decode_json($response->decoded_content);
  my @tracks = @{$track_list->{tracks}};
  subtest "$species tracks" => sub {
     test_track($track_list_url, $_) for @tracks;
  }
}
sub test_track {
  my ($track_list_url, $track_config) = @_;
  my $template = $track_config->{urlTemplate};
  $template =~ s/{refseq_dirpath}.*//;
  $template =~ s/{refseq}.*//;
  my $url = URI->new_abs($template, $track_list_url);
  return if $url =~ /ngs/;
  my $response = LWP::UserAgent->new->head($url);
  ok($response->is_success, "HEAD $url ok");
}
