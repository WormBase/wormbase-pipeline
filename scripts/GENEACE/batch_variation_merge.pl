#!/usr/bin/env perl
#
# script to merge variations in the form:
#   batch_variation_merge.pl -user mt3 -file FILE.TXT [-debug mt3]
#
#   the File is of format:
#   WBVar1	WBVar2	WBVar3
#   variation_to_stay	variations_to_merge


use Getopt::Long;
use IO::File;

use lib $ENV{CVS_DIR};
use lib "$ENV{CVS_DIR}/NAMEDB/lib";
#use lib '/nfs/users/nfs_g/gw3/Nameserver-API';

use NameDB_handler;
use Log_files;
use Wormbase;
use strict;

my ($test, $help, $debug, $verbose, $store, $wormbase, $species);
my ($file, $debug);
GetOptions(
	   "test"       => \$test,
           'file:s'     => \$file,
           'debug:s'    => \$debug,
	   "store:s"    => \$store,
	   "species:s"  => \$species,
)||die(@!);

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -organism => $species,
			     -test => $test,
                             );
}

my $log;
if (defined $debug) {$log = Log_files->make_log("NAMEDB:$file", $debug);}
else {$log = Log_files->make_log("NAMEDB:$file");}

$log->write_to("Contacting NameServer .....\n");
my $db = NameDB_handler->new($wormbase, $test);

my (@merges, %all_ids);

my $inf = IO::File->new($file,'r')||die(@!);
while (<$inf>){
  chomp;
  my ($target, @to_merge) = split;

  #
  # check that none of the ids have been seen previously in the file
  #
  if (grep { exists $all_ids{$_} } ($target, @to_merge)) {
    $db->dienice("VALIDATION FAILED: Duplicate ids in list");
  }
  map { $all_ids{$_} = 1 } ($target, @to_merge);
  
  #
  # Check that to-merge list does not include target!
  #
  if (grep { $_ eq $target } @to_merge) {
    $db->dienice("VALIDATION FAILED: merge list (@to_merge) includes target variation ($target)");
  }

  #
  # Check that we are dealing exclusively with WBVar ids
  #
  foreach my $var ($target, @to_merge) {
    if ($var !~ /^WBVar/) {
      $db->dienice("VALIDATION FAILED: You have a non-WBVar id in your list ($var)");
    }
    if (scalar @{$db->find_variations($var)} == 0) {$db->dienice("VALIDATION FAILED: $var not found in Nameserver")};
  }

  push @merges, [$target, @to_merge];

}

#
# Validation complete. Can now do actual merges with (semi) confidence
#
foreach my $merge_arr (@merges) {
  my ($target, @to_merge) = @$merge_arr;

  # we have decided that merging variations is done so rarely that we
  # are not even going to make a REST entry-point for it and so we
  # simply deleted the variations that are to be merged into the
  # target
  my $ids = $db->kill_variations(\@to_merge);

}

$log->mail;

