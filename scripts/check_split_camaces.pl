#!/usr/local/bin/perl5.8.0 -w
#
# check_split_camaces.pl
#
# Cronjob integrity check controls for split camace databases.
#
# Last updated by: $Author: pad $
# Last updated on: $Date: 2006-03-14 10:23:54 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Log_files;
use Getopt::Long;
use Storable;

my ($debug,$test,);

GetOptions(
	   'debug:s'        => \$debug,
	   'test'           => \$test,
	  );

my $wormbase = Wormbase->new(
    -test    => $test,
    -debug   => $debug,
);


if ($debug) {
my $maintainers = "$debug";
}
else {
my $maintainers = "All";
}


my $rundate = $wormbase->rundate;
my $runtime = $wormbase->runtime;
my $path = "/nfs/disk100/wormpub";
my $age = 1;
my $maintainers = "All";

# is today monday?  If so set age to be 3 to ignore weekend
$age = 3 if (`date +%u` == 1);

my $log = Log_files->make_build_log();
$log->write_to("$runtime : script Started\n");
my @users = ("gw3", "pad");

foreach my $user (@users) {

  $log->write_to("Processing camace_${user}:\n");

  my $dbpath = $path."/camace_".$user."/database";
  
  if(-e "${dbpath}/ACEDB.wrm" ) {
    my (@date) =`ls -ltr $dbpath/block*.wrm`;

    # get the last block file in @date and split to find date
    my @line = split(/\s+/,$date[$#date]);
    
    my $file_age = sprintf("%.1f", -M $line[8]);
    $log->write_to("last modified $file_age days ago: ");
    # Was file modified in last day?
    if (-M $line[8] >$age){
      $log->write_to("no need to re-run camcheck.pl\n\n");
    }
    else{
      $log->write_to("running camcheck.pl\n\n");
      $wormbase->run_script("camcheck.pl -db $path/camace_${user} -l -e $user" , $log);
    }
  }
  else{
    $log->write_to("ERROR: ${dbpath}/ACEDB.wrm file not present\n");
  }
}

$log->write_to("$runtime : script finished\n");

$log->mail;
exit(0);

  

