#!/usr/local/bin/perl -w


use strict;
use lib '../lib';
use NameDB;
use Getopt::Long;

my $USAGE = <<END;
Usage: $0 <options> domain [nametype=]name_pattern ...
  Wildcard search for names.

Options:

  --database  DBI-style database dsn, e.g. mysql:test:host=localhost
  --user      username
  --password  password
  --type      name type
  --brief

Options can be abbreviated, as in:

  namedb_ls.pl -d mysql:test -u fred -p secret Gene CGC=unc*
END

my ($DB,$USER,$PASS,$BRIEF,$TYPE);
GetOptions('database:s' => \$DB,
	   'user:s'     => \$USER,
	   'password:s' => \$PASS,
	   'brief'      => \$BRIEF,
	   'type:s'       => \$TYPE,
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

unshift @ARGV,'*' unless @ARGV;

while (my $name = shift) {
  my @results;
  if ($name=~ /(.+)=(.+)/) {
    my ($type,$pattern) =  ($1,$2);
    @results = $db->idSearch($type=>$pattern);
  } else {
    @results = $db->idSearch($name);
  }
  unless (@results) {
    warn "$name: not found\n";
    next;
  }
  print_names($db,$_) foreach @results;
}

sub print_names {
  my ($db,$id) = @_;
  if ($TYPE) {
    my @names = $db->idTypedNames($TYPE => $id);
    push @names,$id unless @names;
    print "@names\n";
  } else {
    my %typed_names = $db->idAllNames($id);
    my @names = map {ref($_) ? @$_ : $_} values %typed_names;
    print $BRIEF ? "$id\n" : "$id @names\n";
  }
}

