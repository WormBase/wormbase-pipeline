#!/usr/local/bin/perl -w


use strict;
use lib "$ENV{'CVS_DIR'}/NAMEDB/lib";
use NameDB;
use Getopt::Long;

my $USAGE = <<END;
Usage: $0 <options> domain nametype1=name1 nametype2=name2...
  Create a new ID with the indicated names attached.

Options:

  --database  DBI-style database dsn, e.g. mysql:test:host=localhost
  --user      username
  --password  password

Options can be abbreviated, as in:

  namedb_newid.pl -d mysql:test -u fred -p secret Gene cgc=unc-3 wtp=3j45
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
my $id = $db->idCreate;
while (my $name = shift) {
  unless ($name=~ /(.+)=(.+)/) {
    warn "$name: invalid pattern, use type=name\n";
    next;
  }
  my($type,$name) = ($1,$2);
  unless (eval {$db->addName($id,$type,$name)}) {
    warn $@ =~ /Errcode:\s+6/ ? "bad type $type: not entered" : $@;
  };
}
print $id,"\n";

