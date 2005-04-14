#!/usr/local/bin/perl

# ms2@sanger.ac.uk

# extracts the lines from a worm acedb gff file being either cosmids or yacs

use strict;

while (<>) {
    chomp;
    if ((/^\S+\tGenomic_canonical\tregion\s+/) && (/\"(\S+)\"$/)) {
        print "$_\n";
    }
    elsif ((/^\S+\tLink\tregion\s+/) && (/\"(\S+)\"$/)) {
        print "$_\n";
    }
}

