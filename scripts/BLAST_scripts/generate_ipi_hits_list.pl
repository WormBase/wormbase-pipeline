#!/usr/local/ensembl/bin/perl


BEGIN {
    unshift (@INC , "/nfs/acari/wormpipe/scripts/BLAST_scripts");
}

use strict;

my $wormpipe = glob("~wormpipe");
my $wormpipe_dump = "/acari/work2a/wormpipe/dumps";

#Pep_homol "TR:Q9BZ18" wublastp_human 36.337242 765 939 29 195 Align 886 142

open (DATA,"cat $wormpipe_dump/blastp_ensembl.ace $wormpipe_dump/blastx_ensembl.ace |");
my %uniquify;
while (<DATA>) {
  if( /human/ ) {
    next if /WP:CE/; # reciprocal match
    s/\"//g;
    my @data = split(/\s+/,$_);
    $uniquify{ (split(/:/,$data[1]))[1] } = 1;
  }
}
close DATA;  


open (LIST,">$wormpipe_dump/ipi_hits_list") or die "cant open output file $wormpipe_dump/ipi_hits_list";

foreach( keys %uniquify ) {
  print LIST "$_\n";
}
close LIST;

exit(0);
