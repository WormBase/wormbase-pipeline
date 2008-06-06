#!/usr/local/bin/perl -w

use lib '/software/worm/lib/perl/';
use lib '/nfs/WWWdev/SANGER_docs/lib/Projects/C_elegans/';
use strict;
use lib '../blib/lib';
use lib '../lib';
use NameDB;
use Getopt::Long;

my $USAGE = <<END;
Usage: $0 <options> domain id nametype1=name1 nametype2=name2...
 Delete the indicated names.

Options:

  --database  DBI-style database dsn, e.g. mysql:test:host=localhost
  --user      username
  --password  password

  id_resurrect.pl -user ar2 -password <ar2 pswd> -id WBGene00000001 -domain Gene 
  optional  (-test -database).  Will default to live nameserver, -test will change to test database.
  
  Valid DOMAINS are domains currently in the nameserver (ie 'Gene' and 'Variation' at time of writing)  
  
END

my ($DB,$USER,$PASS);
my ($domain, $id, $test);
GetOptions('database:s' => \$DB,
	   'user:s'     => \$USER,
	   'password:s' => \$PASS,
	   'id:s'       => \$id,
	   'test'       => \$test,
	   'domain:s'   => \$domain
	  ) or die $USAGE;

unless ($DB) {
	$DB = 'wbgene_id;mcs2a;3305';
	$DB = 'test_wbgene_id;mcs2a;3305' if $test;
}
$USER ||= $ENV{NAMEDB_USER};
$PASS ||= $ENV{NAMEDB_PASS};

$DB or die "Must provide a --database option";
$USER ||= $ENV{USER};

die $USAGE unless ($id and $domain and $USER and $PASS and $DB);
my $db = NameDB->connect($DB,$USER,$PASS);

my @doms = $db->getDomains;
#die ("$domain is invalid\n Use one of".join(@doms,', ')."\n") unless (grep $domain @doms);

$db->setDomain($domain);

# Doesn't exist.  Try to look up
unless ($db->idExists($id)) {
  die "$id: not found\n";
}

die( "$id is already alive!\n") if $db->idLive($id);

if ($db->idResurrect($id)) {
	print "$id resurrected\n";
}

exit(0);

