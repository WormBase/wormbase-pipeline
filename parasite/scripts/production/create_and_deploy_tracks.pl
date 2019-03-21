# /usr/bin/env perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

sub exit_with_usage {
  print STDERR "Usage: EXPRESSION_CODE=... $0 <options>";
  exit 1;
}
exit_with_usage unless $ENV{EXPRESSION_CODE};

my ( $debug, @species, @notspecies, $allspecies, $setup, $dna, $genes, $rules, $inputids, $meta, $pipeline_setup, $test, $yfile );

GetOptions(
  'species=s@'    => \@species,
  'notspecies=s@' => \@notspecies,
  'allspecies'    => \$allspecies,
  'setup'         => \$setup,
  'load_meta'     => \$meta,
  'load_dna'      => \$dna,
  'load_genes'    => \$genes,
  'load_pipeline' => \$pipeline_setup,
  'load_rules'    => \$rules,
  'load_iids'     => \$inputids,
  'debug'         => \$debug,
  'test'          => \$test,
  'yfile=s'       => \$yfile,

) || exit_with_usage;



use lib "$ENV{EXPRESSION_CODE}/lib";
use lib "$ENV{EXPRESSION_CODE}/local/lib/perl5";
use WbpsExpression;

print "ok\n"; 

