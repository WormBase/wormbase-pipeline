#!/usr/local/bin/perl -w


use strict;
use lib '../blib/lib';
use NameDB;
use Getopt::Long;

my $USAGE = <<END;
Usage: $0 <options> domain name1 name2 nametype=name3...
  Look up a list of names

Options:

  --database  DBI-style database dsn, e.g. mysql:test:host=localhost
  --user      username
  --password  password

Options can be abbreviated, as in:

  namedb_lookup.pl -d mysql:test -u fred -p secret Gene unc-5 cgc=unc-3 wtp=3j45
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

while (my $name = shift) {
  my $id;
  if ($name=~ /(.+)=(.+)/) {
    $id = $db->idGetByTypedName($1,$2);
  } else {
    $id = $db->idGetByAnyName($name);
  }
  unless (@$id) {
    warn "$name: not found\n";
    next;
  }
  print_names($db,$_) foreach @$id;
}

sub print_names {
  my ($db,$id) = @_;
  my %typed_names = $db->idAllNames($id);

  print "$id:\n";
  for my $type (keys %typed_names) {
    my @values = ref $typed_names{$type} ? @{$typed_names{$type}} : $typed_names{$type};
    print "\t",$type,": ","@values","\n";
  }
}

