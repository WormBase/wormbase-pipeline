#!/usr/local/ensembl/bin/perl

use strict;
my $ace_dir = "/wormsrv2/current_DB";
my $wquery_dir = "/wormsrv2/autoace/wquery";
my $tace = "/nfs/disk100/acedb/RELEASE.SUPPORTED/bin.ALPHA_4/tace";

$ENV{'ACEDB'} = $ace_dir;
my $command=<<EOF;
Table-maker -p "$wquery_dir/accession2clone.def"
quit
EOF
    
open (TACE, "echo '$command' | $tace |");
while (<TACE>) {
    chomp;
    if (/\"(\S+)\"\s+\"(\S+)\"\s*\d*/) {
        print "$2\t$1\n";
    }
}
close TACE;
