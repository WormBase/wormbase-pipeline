#!/usr/local/bin/perl -w


use strict;
use lib '../blib/lib';
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
use NameDB;
use Getopt::Long;

my $USAGE = <<END;
Usage: $0 <options> domain type isUnique isPrimary
  Add a type to the indicated domain.
  isUnique and isPrimary can be "yes" or "true".  Anything
  else (including empty) is "no" by default

Options:

  --database  DBI-style database dsn, e.g. mysql:test:host=localhost
  --user      username
  --password  password

Options can be abbreviated, as in:

  namedb_add_type.pl -d mysql:test -u fred -p secret Protein EMBL yes yes
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

@ARGV >= 2 or die $USAGE;

my ($domain,$type,$isUnique,$isPrimary) = @ARGV;
my $unique  = defined($isUnique) && $isUnique   =~ /yes|true/i;
my $primary = defined($isPrimary) && $isPrimary =~ /yes|true/i;

$db->setDomain($domain);
$db->addNameType($type,$unique,$primary);

