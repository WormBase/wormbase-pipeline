#!/usr/local/bin/perl5.6.1 -w 
# Last updated by: $Author: wormpipe $     
# Last updated on: $Date: 2004-01-12 10:03:20 $      

use strict;
use Getopt::Long;
use Carp;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $database, $ver);
my $maintainers = "All";
our $log;

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "version"     => \$ver,
	    "database=s"  => \$database);


# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

##########################
# MAIN BODY OF SCRIPT
##########################

# Do SlimSwissProt
`cat ~pubseq/data/swall/splitted_files/sprot.dat | ~/scripts/BLAST_scripts/swiss_trembl2dbm.pl -s`
`cat ~pubseq/blastdb/swall-1 ~pubseq/blastdb/swall-2 | ~/scripts/BLAST_scripts/swiss_trembl2slim.pl -s $ver`
`fasta2gsi.pl -f /acari/work2a/wormpipe/swall_data/slimswissprot`
`cp /acari/work2a/wormpipe/swall_data/slimswissprot ~/BlastDB/slimswissprot$ver.pep`


# SlimTrEMBL
`cat ~pubseq/data/swall/splitted_files/sprot.dat | ~/scripts/BLAST_scripts/swiss_trembl2dbm.pl -s`
`cat ~pubseq/blastdb/swall-1 ~pubseq/blastdb/swall-2 | ~/scripts/BLAST_scripts/swiss_trembl2slim.pl -s $ver`

# Close log files and exit

&mail_maintainer("script template",$maintainers,$log);
close(LOG);
exit(0);

