#!/usr/local/bin/perl -w


use strict;
use lib '../blib/lib';
use NameDB;
use Getopt::Long;

my $USAGE = <<END;
Usage: $0 <options> domain id nametype1=name1 nametype2=name2...
 Delete the indicated names.

Options:

  --database  DBI-style database dsn, e.g. mysql:test:host=localhost
  --user      username
  --password  password

Options can be abbreviated, as in:

  id_delname.pl -d mysql:test -u fred -p secret Gene P0000001 cgc=unc-3 wtp=3j45
END

my ($DB,$USER,$PASS);
GetOptions('database:s' => \$DB,
	   'user:s'     => \$USER,
	   'password:s' => \$PASS
	  ) or die $USAGE;

$DB   ||= $ENV{NAMEDB_DB};
$USER ||= $ENV{NAMEDB_USER};
$PASS ||= $ENV{NAMEDB_PASS};

$DB or die "Must provide a --database option";
$DB = "dbi:$DB" unless $DB =~ /^dbi:/;
$USER ||= $ENV{USER};
my $db = NameDB->connect($DB,$USER,$PASS);

my $domain = shift or die $USAGE;
my $id     = shift or die $USAGE;
$db->setDomain($domain);

# Doesn't exist.  Try to look up
unless ($db->idExists($id)) {
  my @names = $db->idGetByAnyName($id);
  die "$id: not a unique name\n" if @names > 1;
  die "$id: not found\n" if @names == 0;
  $id = $names[0];
}

while (my $name = shift) {
  $name=~ /(.+)=(.+)/ or die "$name: invalid type=name pair\n";
  $db->delName($id,$1=>$2) or die "$name: couldn't add\n";
}

