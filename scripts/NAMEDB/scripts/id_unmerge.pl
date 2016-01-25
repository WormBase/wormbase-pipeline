#!/software/bin/perl -w

use strict;
use lib $ENV{'CVS_DIR'}."/NAMEDB/lib";
use NameDB;
use Getopt::Long;

my $USAGE = <<END;

Split the identifier given by source_id into target_id, where
target_id had previously been merged into source_id. This
has the effect of undoing a previous erroneous merge, resurrecting
the split target

Options:

  --database  DBI-style database dsn, e.g. mysql:test:host=localhost
  --user      username
  --password  password

  id_unmerge.pl -user ar2 -password <ar2 pswd> -sourceid WBGene00000001 -targetid WBGene00000002 -domain Gene 
  optional  (-test -database).  Will default to live nameserver, -test will change to test database.
  
  Valid DOMAINS are domains currently in the nameserver (ie 'Gene' and 'Variation' at time of writing)  
  
END

my ($DB,$USER,$PASS);
my ($domain, $source_id, $target_id, $test);
GetOptions('database:s'  => \$DB,
	   'user:s'      => \$USER,
	   'password:s'  => \$PASS,
	   'source_id:s' => \$source_id,
	   'target_id:s' => \$target_id,
	   'test'        => \$test,
	   'domain:s'    => \$domain
	  ) or die $USAGE;

unless ($DB) {
	$DB = 'wbgene_id;shap;3303';
	$DB = 'test_wbgene_id;utlt-db;3307' if $test;
}
$USER ||= $ENV{NAMEDB_USER};
$PASS ||= $ENV{NAMEDB_PASS};

$DB or die "Must provide a --database option";
$USER ||= $ENV{USER};

die $USAGE unless ($source_id and $target_id and $domain and $USER and $PASS and $DB);
my $db = NameDB->connect($DB,$USER,$PASS);

my @doms = $db->getDomains;
#die ("$domain is invalid\n Use one of".join(@doms,', ')."\n") unless (grep $domain @doms);

$db->setDomain($domain);

# Doesn't exist.  Try to look up
unless ($db->idExists($source_id)) {
  die "$source_id: not found\n";
}

die( "$source_id is dead!\n") if !$db->idLive($source_id);
die( "$target_id is already alive!\n") if $db->idLive($target_id);

if ($db->idUnmerge($source_id, $target_id, $domain)) {
	print "$target_id is unmerged\n";
}

exit(0);

