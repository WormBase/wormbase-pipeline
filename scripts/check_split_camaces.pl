#!/usr/local/bin/perl5.8.0 -w
#
# check_split_camaces.pl
#
# Cronjob integrity check controls for split camace databases.
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2003-11-28 11:31:38 $

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
my $rundate = &rundate;
my $runtime = &runtime;

my $wormpub_dir = glob("~wormpub");
my $log = "/wormsrv2/logs/check_split_camaces.$rundate.$$";
open (LOG,">$log") or &mail_maintainer("LOG failed in check_split_camaces.pl","wormbase\@sanger.ac.uk");

my @users = ("ar2", "dl1", "pad");

foreach my $name (@users) {
  if(-e "$wormpub_dir/camace_${name}/database/ACEDB.wrm" ) {
      system("/wormsrv2/scripts/camcheck.pl -s $wormpub_dir/camace_${name} -l -e $name");
  }

}

close LOG;
