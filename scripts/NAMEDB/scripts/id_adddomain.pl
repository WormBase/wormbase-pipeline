#!/usr/local/bin/perl -w


use strict;
use NameDB;
use Getopt::Long;

my $USAGE = <<END;
Usage: $0 <options> domain [template]
  Add a domain and template.

Options:

  --database  DBI-style database dsn, e.g. mysql:test:host=localhost
  --user      username
  --password  password

Options can be abbreviated, as in:

  namedb_add_domain.pl -d mysql:test -u fred -p secret
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

my ($domain,$template) = @ARGV or die $USAGE;

$db->addDomain($domain=>$template);

