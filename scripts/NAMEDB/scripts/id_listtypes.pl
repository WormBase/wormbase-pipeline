#!/usr/local/bin/perl -w


use strict;
use lib '../blib/lib';
use NameDB;
use Getopt::Long;

my $USAGE = <<END;
Usage: $0 <options> domain
  List the name types in a naming domain.

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


my $domain = shift or die $USAGE;
$db->setDomain($domain);
my @types = $db->getNameTypes;
printf("%-10s %-7s %-7s\n",'Type','Unique','Primary');
printf("%-10s %-7s %-7s\n",'----','------','-------');
for my $t (@types) {
  my ($unique,$primary) = $db->getNameTypeAttributes($t);
  printf("%-10s %-7s %-7s\n",$t,$unique ? 'yes' : 'no',$primary ? 'yes' : 'no');
}

