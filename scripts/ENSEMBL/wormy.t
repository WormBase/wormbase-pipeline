#!/usr/bin/env perl
use strict;
use warnings;
use Test::More 'no_plan';
use Wormy;
use YAML;

use constant INFILE  => 'test.gff';

ok(1,'test Test::More');
ok(Wormy->new,'new Wormy');
my $testworm=Wormy->new;
ok($testworm->get_lines(INFILE),'trying to parse the file');
ok($testworm->parse_lines,'parsing into lvl1 features');

print YAML::DumpFile('test.yml',$testworm);
