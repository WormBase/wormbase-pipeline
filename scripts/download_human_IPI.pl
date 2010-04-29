#!/usr/local/bin/perl5.6.1 -w
# 
# download_human_IPI.pl
#
# by Chao-Kung Chen
#
# Script to run consistency checks on the geneace database
#
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2010-04-29 14:29:59 $


use strict;

my $rundate = `date +%m_%d`; chomp $rundate;

chdir "/lustre/scratch101/ensembl/wormpipe/BlastDB";

`wget ftp://ftp.ebi.ac.uk/pub/databases/IPI/current/ipi.HUMAN.fasta.gz`;

system("gunzip ipi.HUMAN.fasta.gz; mv ipi.HUMAN.fasta ipi_human_$rundate");

exit(0);


