#!/usr/bin/env perl
#
# script to merge variations in the form:
#   batch_variation_merge.pl -user mt3 -file FILE.TXT [-server test_wbgene_id;utlt-db;3307'] [-debug mt3] [-test]
#
#   the File is of format:
#   WBVar1	WBVar2	WBVar3
#   variation_to_stay	variations_to_merge


use Getopt::Long;
use IO::File;

use lib $ENV{CVS_DIR};
use lib "$ENV{CVS_DIR}/NAMEDB/lib";
use NameDB_handler;
use Log_files;

use strict;

my ($user,$pass,$test,$file,$server,$debug);
GetOptions(
           'user:s'     => \$user,
           'password:s' => \$pass,
           'test'       => \$test,
           'file:s'     => \$file,
           'server'     => \$server,
           'debug:s'    => \$debug,
)||die(@!);

my $log;
if ($user) {$log = Log_files->make_log("NAMEDB:$file", $user);}
elsif (defined $debug) {$log = Log_files->make_log("NAMEDB:$file", $debug);}
else {$log = Log_files->make_log("NAMEDB:$file");}

my ($DB,$db);

if ($server){
    $DB = $server;
}elsif ($test) {
    $DB = 'test_wbgene_id;utlt-db;3307';
}else {
    $DB = 'nameserver_live;web-wwwdb-core-02;3449';
}

$log->write_to("Contacting NameServer $DB.....\n");
$db = VariationDB_handler->new($DB,$user,$pass);

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
    $db->validate_id($var);
  }

  push @merges, [$target, @to_merge];

}

#
# Validation complete. Can now do actual merges with (semi) confidence
#
foreach my $merge_arr (@merges) {
  my ($target, @to_merge) = @$merge_arr;
  
  foreach my $v(@to_merge){
    if (my $ids = $db->merge_variations($target, $v)){
      $log->write_to("Merge complete, $v is DEAD and has been merged into variation $target\n");
    } else {
      $log->write_to("ERROR: Sorry, the variation merge of $v into $target failed\n");
    }
  }
}

$log->mail;

# to not trample all over the NameDB_handler
package VariationDB_handler;
use parent 'NameDB_handler';

# trimmed down constructor for variations
sub new{
    my ($class,$dsn,$name,$password) = @_;

    my $db = NameDB->connect($dsn,$name,$password);

    bless ($db, $class);

    # prevent it being used with anything but variations
    $db->setDomain('Variation');

    return $db;
}


sub merge_variations {
  my ($db, $variation, $merge_variation) = @_;
  
  $db->remove_all_names($merge_variation);

  if ($db->idMerge($merge_variation,$variation)) {
    return ([$variation,$merge_variation]);
  } else {
    $db->dienice("FAILED : merge failed");
    return undef;
  }
}
