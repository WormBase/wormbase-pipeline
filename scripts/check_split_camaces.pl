#!/usr/local/bin/perl5.6.1 -w
#
# check_split_camaces.pl
#
# Cronjob integrity check controls for split camace databases.
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2002-12-05 10:55:10 $

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;

my $wormpub_dir = glob("~wormpub");

open (LOG,">/wormsrv2/scripts/check_split_camaces") or &mail_maintainer("check_split_camaces.pl didnt work","ar2\@sanger.ac.uk");

opendir (DIR,"$wormpub_dir") or die "cant read $wormpub_dir\n";
my @ll =  readdir (DIR);
close DIR;

my @databases;

foreach my $database (@ll) {
  if( $database =~ /camace_(\w+)/ ) {
    if(-e "$wormpub_dir/$database/database/ACEDB.wrm" ) {
      system("/wormsrv2/scripts/camcheck.pl -s $wormpub_dir/$database -l -e $1") unless ("$1" eq "orig");
    }
  }
}
