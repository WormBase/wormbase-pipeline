#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';
use Path::Class;

my $fh = file($ARGV[0])->openr;
while (my $line = $fh->getline()) {
    chomp $line;
    my ($id) = $line =~ /:(PPA\d+)/;
    print $id . "\n";
}
