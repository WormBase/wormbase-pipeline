#!/usr/local/bin/perl5.8.0 -w 
#
# initiate_build.pl
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2004-11-24 15:06:07 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Coords_converter;
use File::Copy;
use File::Spec;

my ($test,$debug,$database);

GetOptions (
	    'test'       => \$test,
	    'debug:s'    => \$debug,
	    'database:s' => \$database
	   );

$test = 1 if ( defined $ENV{'TEST_BUILD'} and $ENV{'TEST_BUILD'} == 1);

# is this script being run as user wormpub???
&test_user_wormpub unless ( $test or $debug );

# establish log file.
my $log = Log_files->make_build_log($debug);

#################################################################################
# initiate autoace build                                                        #
#################################################################################

my $cvs_file = "$database/wspec/database.wrm";
    
# get old build version number, exit if no WS version is returned
# add 1 for new build number, but just use '666' if in test mode

my $WS_new_name;
my $WS_version;

if ($test) {
  $WS_version = "666";
  $WS_new_name = "666";
} else {
  $WS_version = &get_wormbase_version;
  $WS_new_name = $WS_version +1;
}
$log->log_and_die("No_WormBase_release_number") unless defined($WS_version);


# make new build_in_process flag
    
# make sure that the database.wrm file is younger
sleep 10;

# update database.wrm using cvs
&run_command("sed 's/WS${WS_version}/WS${WS_new_name}/' < $cvs_file > ${cvs_file}.new");
my $status = move("$database/wspec/database.wrm.new", "$cvs_file") or $log->write_to("ERROR: renaming file: $!\n");
$log->write_to("ERROR: Couldn't move file: $!\n") if ($status == 0);

# update CVS
unless ($test) {
  system("cd $database/wspec; cvs update");
  system("cd ${database}_config ;cvs update");
}
else {
  print "NOT updating cvs as in test mode\n";
}


# check that top doesn't reveal strange processes running
my $top = `top`;

# add lines to the logfile
my $msg = "Updated WormBase version number to WS$WS_new_name\n";
$msg .= "You are ready to build another WormBase release\n";
$msg .= "Please tell camace and geneace curators to update their database to use the new models!!!\n\n";
$msg .= "Please also check following 'top' output to see if there are stray processes that should\n";
$msg .= "be removed:\n$top\n\n";

$log->write_to($msg);

# this will force a refresh of the coordinate files.
#my $coords = Coords_converter->invoke($database,1);

$log->mail;
exit(0);

#__ end initiate_build __#
__END__
__BUILD_INFO__

<p class=indent>
<a href="http://intweb.sanger.ac.uk/Projects/C_elegans/PERL/autoace_minder.txt.shtml">
<font color=red>autoace_minder.pl</a></font><font color=red> -initial</font> (-test)<br>
+ Starts a new build.<br>
+ Updates WS version number in <font color=green>/wormsrv2/autoace/wspec/database.wrm</font>.<br>
+ Writes <font color=green>/wormsrv2/autoace/logs/A1:Build_in_progress</font> lock file<br>
+ Writes <font color=green>/wormsrv2/autoace/logs/A2:Updated_WS_version</font> lock file on completion<br>
<b>~1 minute</b>
</p>

<p class=check>
<font color=blue><B>CHECK:</B></font><br>
The <font color=green>A1:Build_in_progress</font> lock file produced must predate the version file 
<font color=green>/wormsrv2/autoace/wspec/database.wrm</font> for all operations in this build. 
Fatal failures (which lead to the modification of the flag file) will prevent any further runs of 
autoace_minder.pl. The build procedure must be started from the beginning by manual roll-back of the 
version file and running <font color=red>autoace_minder.pl -initial</font> again.
</p>
