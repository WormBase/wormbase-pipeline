#!/usr/local/bin/perl5.8.0 -w
#
# check_split_camaces.pl
#
# Cronjob integrity check controls for split camace databases.
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2005-12-16 11:18:55 $

use strict;
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Wormbase;

my $rundate = &rundate;
my $log = "/nfs/disk100/wormpub/logs/check_split_camaces.$rundate.$$";
my $path = "/nfs/disk100/wormpub";
my $age = 1;
my $maintainers = "All";

# is today monday?  If so set age to be 3 to ignore weekend
$age = 3 if (`date +%u` == 1);

open (LOG,">$log") or &mail_maintainer("LOG failed in check_split_camaces.pl","wormbase\@sanger.ac.uk");

print LOG &runtime, ": script started\n\n";

my @users = ("gw3", "pad");

foreach my $user (@users) {

  print LOG "Processing camace_${user}:\n";

  my $dbpath = $path."/camace_".$user."/database";
  
  if(-e "${dbpath}/ACEDB.wrm" ) {
    my (@date) =`ls -ltr $dbpath/block*.wrm`;

    # get the last block file in @date and split to find date
    my @line = split(/\s+/,$date[$#date]);
    
    my $file_age = sprintf("%.1f", -M $line[8]);
    print LOG "last modified $file_age days ago: ";
    # Was file modified in last day?
    if (-M $line[8] >$age){
      print LOG "no need to re-run camcheck.pl\n\n";
    }
    else{
      print LOG "running camcheck.pl\n\n";
      system("$ENV{'CVS_DIR'}/camcheck.pl -s $path/camace_${user} -l -e $user");
    }
  }
  else{
    print LOG "ERROR: ${dbpath}/ACEDB.wrm file not present\n";
  }
}

print LOG &runtime, ": script finished\n";

close(LOG);

&mail_maintainer("check_split_camaces.pl Report:",$maintainers,$log);

exit(0);

  

