#!/usr/bin/env perl
#
# script to merge variations in the form:
#   batch_variation_merge.pl -user mt3 -file FILE.TXT [-server test_wbgene_id;mcs12;3307'] [-debug mt3] [-test]
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
    $DB = 'test_wbgene_id;mcs12;3307';
}else {
    $DB = 'wbgene_id;shap;3303';
}

$log->write_to("Contacting NameServer $DB.....\n");
$db = VariationDB_handler->new($DB,$user,$pass);

my $inf = IO::File->new($file,'r')||die(@!);
while (<$inf>){
  chomp;
  my @variations = split;

  # minimalistic sanity check
  unless(scalar(@variations) == scalar(grep{/WBVar\d+/}@variations)){
      $log->write_to("ERROR: $_ contains a non WBVarID ... skipping the line\n");
      next;
  }

  my $variation = shift @variations;
  foreach my $v(@variations){
   if (my $ids = $db->merge_variations($variation, $v)){
       $log->write_to("Merge complete, $v is DEAD and has been merged into variation $variation\n");
   }else {
       $log->write_to("ERROR: Sorry, the variation merge of $v into $variation failed\n");
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

# the merge method
# my $var_id = $db->merge_variations($id_to_keep, $id_to_kill);
sub merge_variations {
  my ($db, $variation, $merge_variation) = @_;
  
  my $variation_id  = ($db->idGetByAnyName($variation)->[0]  || $variation); #allow use of any name or variatiom_id
  my $merge_id = ($db->idGetByAnyName($merge_variation)->[0] || $merge_variation); #allow use of any name or variation_id
  $db->validate_id($variation_id) or return undef;
  $db->validate_id($merge_id) or return undef;
  
  if ( $variation_id eq $merge_id) {
    $db->dienice("FAILED : $variation and $merge_variation are the same!");
    return undef;
  }
  
  #always remove names from merged id
  $db->remove_all_names($merge_id);

  #do the merge
  if ($db->idMerge($merge_id,$variation_id)) {
    return ([$variation_id,$merge_id]);
  } else {
    $db->dienice("FAILED : merge failed");
    return undef;
  }
}
