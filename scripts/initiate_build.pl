#!/usr/local/bin/perl5.8.0 -w 
#
# initiate_build.pl
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2006-01-05 17:27:34 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Coords_converter;
use File::Copy;
use File::Spec;
use Storable;

my ($test,$debug,$database, $version);
my ($store, $wormbase);
GetOptions (
	    'test'       => \$test,
	    'debug:s'    => \$debug,
	    'database:s' => \$database,
	    'version:s'  => \$version,
	    'store:s'    => \$store
	   );

if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			     -version => $version
			   );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$log->log_and_die( "version to build not specified\n") unless $wormbase->version;

my $old_ver = $wormbase->version - 1;

#################################################################################
# initiate autoace build                                                        #
#################################################################################

## update CVS wspec, wquery and autoace_config from CVS
$wormbase->run_command("cd ".$wormbase->autoace.";cvs -d :ext:cvs.sanger.ac.uk:/nfs/ensembl/cvsroot/ checkout -d wspec wormbase/wspec", $log);
$wormbase->run_command("cd ".$wormbase->autoace.";cvs -d :ext:cvs.sanger.ac.uk:/nfs/ensembl/cvsroot/ checkout -d wquery wormbase/wquery", $log);
$wormbase->run_command("cd ".$wormbase->basedir.";cvs -d :ext:cvs.sanger.ac.uk:/nfs/ensembl/cvsroot/ checkout -d autoace_config wormbase/autoace_config",$log);

## make new build_in_process flag ( not done yet in rebuild )

## update database.wrm using cvs
#$wormbase->run_command("sed 's/WS${old_version}/WS${version}/' < $cvs_file > ${cvs_file}.new");
#my $status = move("$database/wspec/database.wrm.new", "$cvs_file") or $log->write_to("ERROR: renaming file: $!\n");
#$log->write_to("ERROR: Couldn't move file: $!\n") if ($status == 0);

# update CVS
unless ($wormbase->test) {
  system("cd $database/wspec; cvs update");
  system("cd ${database}_config ;cvs update");
}
else {
  print "NOT updating cvs as in test mode\n";
}

# add lines to the logfile
my $msg = "Updated WormBase version number to WS".$wormbase->version."\n";
$msg .= "You are ready to build another WormBase release\n";
$msg .= "Please tell camace and geneace curators to update their database to use the new models!!!\n\n";

$log->write_to($msg);
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
