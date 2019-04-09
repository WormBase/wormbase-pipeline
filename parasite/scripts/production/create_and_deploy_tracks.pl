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
my $jbrowse_tools_skip = 0;
my $jbrowse_tools_nofatal = 0;

GetOptions(
  'species=s@'    => \@species,
  'do_sync' => \$do_sync,
  'jbrowse_tools_skip' => \$jbrowse_tools_skip,
  'jbrowse_tools_nofatal' => \$jbrowse_tools_nofatal,
) || exit_with_usage();




for my $species (@species){
  GenomeBrowser::JBrowseDisplay->new->make_tracks(
    $species,
    do_sync=>$do_sync,
    jbrowse_tools_skip => $jbrowse_tools_skip,
    jbrowse_tools_nofatal=> $jbrowse_tools_nofatal,
  );
}

sub exit_with_usage {
  print STDERR "Usage: $0 <params>";
  exit 1;
}

