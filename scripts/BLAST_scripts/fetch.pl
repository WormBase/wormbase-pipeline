#!/usr/local/bin/perl

# ms2@sanger.ac.uk

# extract the appropriate (sub)sequences from the database(es), and write them to STDOUT
# uses Sean Eddy's GSI indexing (requires GSI-indexed databases, called .gsi)
# -o option uses srs to get the species info for every sequence

# hsp file format:
#             (score, e-value, percent identity, query-start, query-end, query, query length,
#              target-start, target-end, target, target length, first iteration this target appeared)

BEGIN {
    unshift (@INC , "/nfs/acari/wormpipe/Pipeline");
}

use strict;
use GSI;
use Getopt::Std;
use vars qw($opt_g $opt_i $opt_p);

my $usage = "fetch.pl.pl\n";
$usage .= "-g [gsi index of database]\n";
$usage .= "-i [sequence id]\n";
$usage .= "-p return only the sequence string  DEFAULT: Fasta\n";

getopts ("g:i:p");

unless ($opt_g && $opt_i) {
    die "$usage";
}

########################
# open the GSI file
GSI::openGSI ("$opt_g");

########################
my $key = $opt_i;
my ($file , $fmt , $offset) = GSI::getOffset ($key);
open (DB , "$file") || die "cannot open db $file\n";
seek (DB , $offset , 0);
my $header = <DB>;
chomp $header;
my $seq = "";
while ((my $line = <DB>) =~ /^[^\>]/) {
    $seq .= $line;
}
close DB;
$seq =~ s/\n//g;
if ($opt_p) {
    print "$seq\n";
}
else {
    print "$header";
    my $rev = reverse $seq;
    my $count = 0;
    while (my $base = chop $rev) {
        if ($count++ % 50 == 0) {
            print "\n";
	}
        print "$base";
    }
    print "\n";
}

########################
# close the GSI file 
GSI::closeGSI
