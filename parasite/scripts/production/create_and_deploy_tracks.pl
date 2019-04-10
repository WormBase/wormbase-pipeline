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
use GenomeBrowser::RnaseqHub;

my @species;
my $jbrowse_skip = 0;
my $hub_skip = 0;
my $sanger_deployment_skip = 0;
my $jbrowse_tools_skip = 0;
my $jbrowse_tools_nofatal = 0;

GetOptions(
  'species=s@'    => \@species,
  'jbrowse_skip' => \$jbrowse_skip,
  'hub_skip' => \$hub_skip,
  'sanger_deployment_skip' => \$sanger_deployment_skip,
  'jbrowse_tools_skip' => \$jbrowse_tools_skip,
  'jbrowse_tools_nofatal' => \$jbrowse_tools_nofatal,
) || exit_with_usage();

for my $species (@species){
  next if $jbrowse_skip;
  GenomeBrowser::JBrowseDisplay->new->make_displays(
    $species,
    sanger_deployment_skip=>$sanger_deployment_skip,
    jbrowse_tools_skip => $jbrowse_tools_skip,
    jbrowse_tools_nofatal=> $jbrowse_tools_nofatal,
  );
}
for my $species (@species){
  next if $hub_skip;
  GenomeBrowser::RnaseqHub->new->make_hubs(
    $species,
    sanger_deployment_skip=>$sanger_deployment_skip,
  );
}

# TODO copy to FTP, too
sub exit_with_usage {
  print STDERR "Usage: $0 <params>";
  exit 1;
}

