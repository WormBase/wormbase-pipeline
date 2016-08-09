#!/usr/local/bin/perl -w
use lib $ENV{'CVS_DIR'};
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
use lib '/software/worm/lib/perl/';

use strict;
use NameDB;
use Getopt::Long;
use Ace;
use Wormbase;
use Storable;


my $USAGE = <<END;
Usage: $0 <options> domain id nametype1=name1 nametype2=name2...
 Delete the indicated names.

Options:

  --database  DBI-style database dsn, e.g. mysql:test:host=localhost
  --user      username
  --password  password
  --seq_name Sequence_name that was removed from geneace

  id_resurrect.pl -user ar2 -password <ar2 pswd> -id WBGene00000001 -domain Gene -species brugia -load
  optional  (-test -database).  Will default to live nameserver, -test will change to test database.
  
  Valid DOMAINS are domains currently in the nameserver (ie 'Gene' and 'Variation' at time of writing)  
  
END

my ($DB,$USER,$PASS);
my ($domain, $id, $test, $output, $species, $load, $force,);
GetOptions('database:s'      => \$DB,
	   'user:s'          => \$USER,
	   'password:s'      => \$PASS,
	   'id:s'            => \$id,
	   'test'            => \$test,
	   'domain:s'        => \$domain,
	   'output:s'        => \$output,
	   'species:s'       => \$species,
	   'load'            => \$load,
	   'seq_name:s'      => \$force,
	  ) or die $USAGE;

unless ($DB) {
  $DB = 'nameserver_live;web-wwwdb-core-02;3449';
  $DB = 'test_wbgene_id;utlt-db;3307' if $test;
}
$USER ||= $ENV{NAMEDB_USER};
$PASS ||= $ENV{NAMEDB_PASS};

$DB or die "Must provide a --database option";
$USER ||= $ENV{USER};

die $USAGE unless ($id and $domain and $USER and $PASS and $DB and $species);

my $wormbase;
$wormbase = Wormbase->new( -test => $test,
			   -organism => $species
			 );

my $rundate     = `date +%y%m%d`; chomp $rundate;
my $tace = $wormbase->tace;

my $WBUSERS = {
	       # these are the users WBPerson id
	       'klh'           => 3111,
	       'pad'           => 1983,
	       'mt3'           => 2970,
	       'gw3'           => 4025,
	       'mh6'           => 4055,
	      };

my $acedb = "/nfs/wormpub/DATABASES/geneace";
my $db = NameDB->connect($DB,$USER,$PASS);

my $ace = Ace->connect('-path', $acedb) or die("cant open $acedb: $!\n");
my $idObj = $ace->fetch('Gene', $id);

my $ver = $idObj->Version->name;
$ver++;

my $clone;
my $seq_name;
if (defined $force) {
  if ($force =~  /(\S+)\.\d+/) {
    $clone = $1;
    $seq_name = $force;
  }
}

  unless  (defined $force)  {
    $seq_name = $idObj->Sequence_name->name;
    if ($seq_name =~ /(\S+)\.\d+/) {
      $clone = $1;
    }
    else {$clone = $seq_name;}
  }


unless ($output) {$output = "/tmp/idressurect_${id}.ace";}
open (OUT,">$output") or die "Can't open $output\n";
my @doms = $db->getDomains;
#die ("$domain is invalid\n Use one of".join(@doms,', ')."\n") unless (grep $domain @doms);

$db->setDomain($domain);

# Doesn't exist.  Try to look up
unless ($db->idExists($id)) {
  die "$id: not found\n";
}

#die( "$id is already alive!\n") if $db->idLive($id);

if ($db->idResurrect($id)) {
	print "$id resurrected\n";
}
print OUT "Gene : \"$id\"\n";
print OUT "Version $ver\n";
print OUT "Live\n";

print OUT "Sequence_name $seq_name\n";
if (defined $$WBUSERS{$USER}) {
print OUT "History Version_change $ver now WBPerson". $$WBUSERS{$USER}." Event Resurrected\n";
}
else {
print OUT "History Version_change $ver now WBPerson1983 Event Resurrected\n";
}

if ($species eq "elegans") {
  print OUT "Positive_clone $clone Inferred_automatically \"From sequence, transcript, pseudogene data\"\n";
}
print OUT "Remark \"[$rundate ${USER}] Gene Resurrected via id_ressurect.pl Additional comments:\"\n";
print OUT "Method Gene";

if ($load) {
  print "Output file created and loaded into $acedb : $output\n";
  &loadace("$output", "id_ressurect", $acedb) or die "Failed to load $output\n";
  #&loadace("$output", "id_ressurect", "/nfs/wormpub/DATABASES/TEST_DATABASES/NewModels") or die "Failed to load $output\n";
}
else {
print "Output file created : $output\nYou need to manuall load this into geneace\n\n";
}

exit(0);



sub loadace {
  my $filepath = shift;
  my $tsuser   = shift;
  my $ldb = shift;  

  my $command = "pparse $filepath\nsave\nquit\n";

  open (TACE, "| $tace $ldb -tsuser $tsuser") || die "Couldn't open $ldb\n";
  print TACE $command;
  close TACE;
  print ("SUCCESS! Loaded $filepath into ${ldb}\n\n");
  return 1;
}
