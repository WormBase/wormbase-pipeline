#!/usr/local/bin/perl -w

use strict;

while (<>) {
    my @f = split /\t/;
    if ($f[4] !~ /N/) { 
        chomp;
        print $_, "\t", $f[5]."_C", "\t", $f[1], "\t", $f[2], "\n";
    }
    elsif ($f[4] =~ /N/) {
        s/[+-]//;
        print $_;
    }
    else {
        print "Problem!\n";
    }
}
