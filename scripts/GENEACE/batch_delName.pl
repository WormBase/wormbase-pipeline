#!/usr/bin/env perl
use lib $ENV{'CVS_DIR'};
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";

use Wormbase;
use NameDB_handler;
use Getopt::Long;
use Log_files;

use warnings;
use strict;

=pod

=head batch_delName.pl

=item Options:

  -user      username
  -password  password
  -file      file containing list of GeneIDs and CGC/Sequence names eg WBGene00008040 ttr-5
  -test      use the test nameserver

e.g. perl batch_delName.pl -u fred -p secret -file genenames.txt

=cut

my ($USER,$PASS, $test, $file, $type,);
GetOptions(
	   'user:s'     => \$USER,
	   'password:s' => \$PASS,
	   'test'       => \$test,
	   'file:s'     => \$file,
	   'type:s'     => \$type,
	  ) or die;

my $log = Log_files->make_log("NAMEDB:$0:$file", $USER);
my $DB;
if ($test) {
    $DB = 'test_wbgene_id;utlt-db;3307';
  } else {
    $DB = 'nameserver_live;web-wwwdb-core-02;3449';
}

$log->write_to("loading $file to $DB\n\n");
$log->write_to("TEST mode is ON!\n\n") if $test;

unless (($type eq  'CGC') || ($type eq 'Sequence')) {
  $log->log_and_die ("$type is not CGC or Sequence\n");
}

my $db = NameDB_handler->new($DB,$USER,$PASS);

$db->setDomain('Gene');

# open file and read
# file needs to be in such a format:
#
# WBGene007 flunk-1
# WBGene123 nope-2

open (FILE,"<$file") or die "can't open $file : $!\n";
my $count=0;
while(<FILE>) {
    my($id,$name,$cds) = split;
    eval{
	my $success = $db->delName($id,$type => $name);
	my $msg = defined $success ? 'deleted' : 'FAILED to delete';
	$log->write_to("$msg $type name $name for $id\n");
        $success = $db->_update_public_name($id);
        $msg = defined $success ? 'updated' : 'FAILED to update';
	$log->write_to("$msg public name to $success for $id\n");
    };
    $count++;
}

$log->write_to("=======================\nprocessed $count genes\n");
$log->mail;
