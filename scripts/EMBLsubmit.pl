#!/usr/local/bin/perl5.6.1 -w
#
# AcePublish.pl
#
# by Steve Jones, 1995
#
# Last updated on: $Date: 2002-12-09 14:33:30 $
# Last updated by: $Author: krb $
#
# This script will mail the contents of ~wormpub/analysis/TO_SUBMIT
# to celegans@ebi.ac.uk.

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;


my $total = 0;
my $dir = "/nfs/disk100/wormpub/analysis/TO_SUBMIT";

opendir(DIR,'$dir') || die "Couldn't open directory\n";
my @filenames=grep(/\.embl$/,readdir(DIR));
close DIR;

foreach my $file (@filenames) {
  print "$file\n";
  $total++;
}

print "\n$total cosmids found in ~wormpub/analysis/TO_SUBMIT\n";
print "\nAll of these will be submitted to celegans\@ebi.ac.uk\n\n";
print "Do you wish to proceed???(y or n)\n\n";

my $answer=<STDIN>;
if ($answer ne "y\n") {die "\nSubmission aborted\n";}


open(FILES, "ls  $dir  |");

while (<FILES>) {
  chop;
  if (/embl$/) { 
    print "Mailing /nfs/disk100/wormpub/analysis/TO_SUBMIT/$_\n";
    system ("/usr/ucb/Mail celegans\@ebi.ac.uk < $dir/$_");}
}
close(FILES);

exit(0);

