#!/usr/local/bin/perl5.8.0 -w
#
# check_split_camaces.pl
#
# Cronjob integrity check controls for split camace databases.
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2004-07-06 09:54:12 $

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;

my $rundate = &rundate;
my $log = "/wormsrv2/logs/check_split_camaces.$rundate.$$";
my $path = "/nfs/disk100/wormpub";
my $age = 1;

# is today monday?  If so set age to be 3 to ignore weekend
$age = 3 if (`date +%u` == 1);

open (LOG,">$log") or &mail_maintainer("LOG failed in check_split_camaces.pl","wormbase\@sanger.ac.uk");


my @users = ("ar2", "dl1", "pad", "krb");

foreach my $user (@users) {

  my $dbpath = $path."/camace_".$user."/database";
  
  if(-e "${dbpath}/ACEDB.wrm" ) {
    my (@date) =`ls -ltr $dbpath/block*.wrm`;

    # get the last block file in @date and split to find date
    my @line = split(/\s+/,$date[$#date]);
    
    # Was file modified in last day?
    if (-M $line[8] >$age){
      print LOG "camace_${user} not modified in last $age day(s), no need to re-run camcheck\n";
    }
    else{
      print LOG "\n\n";
      print LOG &runtime, ": processing camace_${user}\n";
      system("/wormsrv2/scripts/camcheck.pl -s $path/camace_${user} -l -e $user");
      print LOG &runtime, ": Finished\n\n";
    }
  }
  else{
    print LOG "ERROR: ${dbpath}/ACEDB.wrm file not present\n";
  }
}

close(LOG);
exit(0);

  

