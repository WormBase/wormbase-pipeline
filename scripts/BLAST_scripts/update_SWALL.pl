#!/usr/local/ensembl/bin/perl -w 
# Last updated by: $Author: mh6 $     
# Last updated on: $Date: 2006-11-28 16:42:46 $      

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
`fasta2gsi.pl -f /lustre/work1/ensembl/wormpipe/swall_data/slimswissprot`
`cp /lustre/work1/ensembl/wormpipe/swall_data/slimswissprot ~/BlastDB/slimswissprot$ver.pep`


# SlimTrEMBL
`cat ~pubseq/data/swall/splitted_files/sprot.dat | ~/scripts/BLAST_scripts/swiss_trembl2dbm.pl -s`
`cat ~pubseq/blastdb/swall-1 ~pubseq/blastdb/swall-2 | ~/scripts/BLAST_scripts/swiss_trembl2slim.pl -s $ver`

# Close log files and exit

&mail_maintainer("script template",$maintainers,$log);
close(LOG);
exit(0);

