#!/usr/local/bin/perl -w


use strict;
use lib '../blib/lib';
use NameDB;
use Getopt::Long;

my $USAGE = <<END;
Usage: $0 <options> domain id1 id2 id3...
  Lookup a list of primary ids and print a status report.
  Will find identifier even if it is no longer live.
  use namedb_lookup.pl for name lookups.

Options:

  --database  DBI-style database dsn, e.g. mysql:test:host=localhost
  --user      username
  --password  password

Options can be abbreviated, as in:

  namedb_lookup.pl -d mysql:test -u fred -p secret Gene G0000012 G0000029
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

while (my $id = shift) {
  unless ($db->idExists($id)) {
    my @names = $db->idGetByAnyName($id);
    $id = $names[0] if @names == 1;
  }
  unless ($db->idExists($id)) {
    warn "$id: not found\n";
    next;
  }
  print_names($db,$id);
  print "\n";
  print "\t","Status: ",$db->idLive($id) ? "live" : "dead", "\n";
  print_history($db,$id);
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

sub print_history {
  my ($db,$obj) = @_;
  local $^W = 0;
  for my $event (@{$db->idGetHistory($obj)}) {
    print "\t",join "\t",'version',@$event,"\n";
  }
print "\n";
}

