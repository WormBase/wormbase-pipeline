#!/software/bin/perl -w
#
# check_split_camaces.pl
#
# Cronjob integrity check controls for split camace databases.
#
# Last updated by: $Author: pad $
# Last updated on: $Date: 2007-06-13 11:38:14 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Log_files;
use Getopt::Long;
use Storable;

my ($debug,$test, $store);

GetOptions(
	   'debug:s'   => \$debug,
	   'test'      => \$test,
	   'store:s'   => \$store
	  );

my $wormbase;

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

my $log = Log_files->make_build_log($wormbase);

my $rundate = $wormbase->rundate;
my $runtime = $wormbase->runtime;
my $path = glob("~wormpub");
my $age = 1;

# is today monday?  If so set age to be 3 to ignore weekend
$age = 3 if (`date +%u` == 1);

my @users = ("pad", "gw3", "ar2");

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
      $wormbase->run_script("camcheck.pl -database /.automount/evs-users2/root/wormpub/camace_${user} -low -email $user" , $log);
    }
  }
  else{
    $log->write_to("ERROR: ${dbpath}/ACEDB.wrm file not present\n");
  }
}

$log->write_to("$runtime : script finished\n");

$log->mail;
exit(0);

  

