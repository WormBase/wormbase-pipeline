#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use WWW::Mechanize;
my $mech = WWW::Mechanize->new();

my $base_url = "http://test.parasite.wormbase.org";


$mech->get("$base_url/rest/");
my @documentation_page_links = $mech->find_all_links(url_regex => qr{documentation});
ok(($mech->success and @documentation_page_links), "Links to documentation on $mech->uri");
for my $documentation_url (map {join ("/", $base_url , "rest", $_->url) } grep {$_->url !~ /post/} @documentation_page_links){
  my ($section) = $documentation_url =~ m{/(\w+)$};
  subtest "REST $section" => sub {
    $mech->get($documentation_url);
    my @example_links = $mech->find_all_links(url_regex => qr{content-type});
    if ($documentation_url =~ /lookup_genome/){
      # https://www.ebi.ac.uk/panda/jira/projects/PARASITE/issues/PARASITE-354      
      ok($mech->success, "TODO $documentation_url examples");
      return;
    } 
    ok(($mech->success and @example_links), "$documentation_url has links to examples");
    for my $example_url (map {join ("/", $base_url , "rest", $_->url) } @example_links){
       next if $example_url =~ m{/genomes/taxonomy/}; # https://www.ebi.ac.uk/panda/jira/projects/PARASITE/issues/PARASITE-353
       $mech->get($example_url);
       ok($mech->success, $example_url);
    }
  } 
}
subtest "REST POST doc pages" => sub {
  for my $documentation_url (map {join ("/", $base_url , "rest", $_->url) } grep {$_->url =~ /post/} @documentation_page_links){
     $mech->get($documentation_url);
     ok($mech->success, "$documentation_url loads fine");
  }
};

done_testing;
