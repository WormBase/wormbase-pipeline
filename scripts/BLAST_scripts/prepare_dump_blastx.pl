#!/usr/local/bin/perl5.6.0 -w

# Marc Sohrmann (ms2@sanger.ac.uk)
# Amended by Keith Bradnam (krb@sanger.ac.uk)
# runs the necessary TableMaker query to extract Sequence->Accession
# information from the current_DB database on wormsrv2

use strict;

####################################################################
# set some parameters
####################################################################

# AceDB database
my $ace_dir = "/wormsrv2/current_DB";
my $wquery_dir = "/wormsrv2/autoace/wquery";
my $tace = "/nfs/disk100/wormpub/ACEDB/bin_ALPHA/tace";

####################################################################
# connect to AceDB using TableMaker,
# populating %accession2name (maps embl accession to contig name)
####################################################################

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

exit 0;
