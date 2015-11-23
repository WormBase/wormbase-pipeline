#!/usr/local/bin/perl -w
use strict;
use lib $ENV{'CVS_DIR'};
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
use NameDB_handler;
use Getopt::Long;
use Log_files;

=pod

=head batch_addname.pl

=item Options:

  -user      username
  -password  password
  -file 	 file containing list of GeneIDs and CGC name eg WBGene00008040 ttr-5
  -species   what species these are for - default = elegans
  -test      use the test nameserver
  -force     bypass CGC name validation check eg to add Cbr-cyp-33E1; use with care!
  -type      allows you to specify Sequence names to add (defaults to CGC if not speciefied)

e.g. perl batch_addname.pl -u fred -p secret -file genenames.txt -species briggsae

=cut

my ($USER,$PASS, $test, $file, $species, $force, $type);
GetOptions(
	   'user:s'     => \$USER,
	   'password:s' => \$PASS,
	   'test'       => \$test,
	   'file:s'     => \$file,
	   'species:s'  => \$species,
	   'force'      => \$force,
	   'type:s'     => \$type,
	  ) or die;

$species = 'elegans' unless $species;
$type = 'CGC' unless $type;

if (($type ne "CGC") && ($type ne "Sequence")) {
  die "can't open $file : $!\n$type is not a valid type (CGC/Sequence)";
}

my $log = Log_files->make_log("NAMEDB:$file", $USER);
my $DB;
if ($test) {
    $DB = 'test_wbgene_id;utlt-db;3307';
  } else {
    $DB = 'wbgene_id;shap;3303';
}

$log->write_to("loading $file to $DB\n\n");
$log->write_to("FORCE mode is ON!\n\n") if $force;
$log->write_to("TEST mode is ON!\n\n") if $test;

my $db = NameDB_handler->new($DB,$USER,$PASS,'/nfs/WWWdev/SANGER_docs/cgi-bin');

$db->setDomain('Gene');

#open file and read

open (FILE,"<$file") or die "can't open $file : $!\n";
my $method = defined($force) ? 'force_name' : 'add_name';
my $count=0;
while(<FILE>) {
    my($id,$name) = split;
    my $success = $db->$method($id,$name,$type,$species);
    if (defined $success) {
      $log->write_to("$id\t$name\tok\n");
      print "$id\t$name\tok\n";
    }
    else {
      $log->write_to("$id\t$name\tFAILED NAMEDB assignment\n");
      print "$id\t$name\tFAILED NAMEDB assignment\n";
    }
  $count++;
}

$log->write_to("=======================\nprocessed $count genes\n");
$log->mail;
















