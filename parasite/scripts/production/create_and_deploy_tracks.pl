#!/usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

use Log::Any::Adapter;
use Log::Log4perl qw(:easy);
Log::Log4perl->easy_init($DEBUG);
Log::Any::Adapter->set('Log4perl');

use lib "/nfs/panda/ensemblgenomes/wormbase/parasite/wbps-expression/lib";
use lib "/nfs/panda/ensemblgenomes/wormbase/parasite/wbps-expression/local/lib/perl5";
use WbpsExpression;

use GenomeBrowser::JBrowseDisplay;
my @species;
my $do_sync = 0;
my $skip_flatfile_to_json = 0;

GetOptions(
  'species=s@'    => \@species,
  'do_sync' => \$do_sync,
  'skip_flatfile_to_json' => \$skip_flatfile_to_json,
) || exit_with_usage();




for my $species (@species){
  GenomeBrowser::JBrowseDisplay->new->make_tracks(
    $species,
    do_sync=>$do_sync,
    skip_flatfile_to_json => $skip_flatfile_to_json,
  );
}

sub exit_with_usage {
  print STDERR "Usage: $0 <species do_sync do_rnaseq do_jbrowse >";
  exit 1;
}

