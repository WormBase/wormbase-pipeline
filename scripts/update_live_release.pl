#!/usr/local/bin/perl5.8.0 -w
#
# update_live_release.pl
#
# by Anthony Rogers
#
# Updates the local webpages in synch with the main website
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2006-05-30 14:00:15 $


use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Carp;
use Getopt::Long;
use Storable;
use Log_files;

my($test, $debug, $store);
my $release;
my $errors = 0;

GetOptions ("release=i" => \$release,
	    "debug:s"   => \$debug,
	    "store:s"   => \$store,
	    "test"      => \$test
	   );

die "you must give a release version ( just numbers eg -release 125 )\n" unless $release;

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}
# establish log file.
my $log = Log_files->make_build_log($wormbase);


my $www = "/nfs/WWWdev/SANGER_docs/htdocs/Projects/C_elegans";

############################################
# update wormpep files
#############################################

my $wormpub_dir      = $wormbase->wormpub."/WORMPEP";
my $wormpep_ftp_root = glob("~ftp/pub/databases/wormpep");
my $wp_ftp_dir       = "$wormpep_ftp_root/wormpep${release}";
my $wormbase_ftp_dir = glob("~ftp/pub/wormbase"); 

# grab from Wormbase.pm subroutine
my @wormpep_files = $wormbase->wormpep_files;

# update new live wormpep release from ftp_site to /disk100/wormpub
foreach my $file ( @wormpep_files ) {
  unlink("$wormpub_dir/${file}_current") || $log->write_to("ERROR: Cannot delete file $wormpub_dir/${file}_current :\t$!\n");
  $wormbase->run_command("scp $wp_ftp_dir/${file}${release} $wormpub_dir/${file}_current", $log);
}
$wormbase->run_command("xdformat -p $wormpub_dir/wormpep_current", $log);
$wormbase->run_command("setdb $wormpub_dir/wormpep_current", $log);
# create a symbolic link from top level wormpep into the release on the wormpep ftp site
foreach my $file ( @wormpep_files ) {
  $wormbase->run_command("cd $wormpep_ftp_root; ln -fs wormpep${release}/${file}${release} $file", $log);
}

# delete the old symbolic link and make the new one to wormpep.prev
my $prev_release = $release -1;
$wormbase->run_command("cd $wormpep_ftp_root; ln -fs wormpep${prev_release}/wormpep${prev_release} wormpep.prev",$log);


##################################################################
# Update Sanger WORMBASE ftp site to change live_release symlink
##################################################################
$wormbase->run_command("cd $wormbase_ftp_dir; rm -f live_release; ln -fs WS${release} live_release",$log);




##############################################
# Update pages using webpublish
# Separate webpublish commands (for safety!) on the two top level directories that need updating
##############################################

my $webpublish = "/usr/local/bin/webpublish";

# Now update WORMBASE current link
$wormbase->run_command("cd $www/WORMBASE; rm -f current; ln -fs WS${release} current", $log) && $log->write_to("Couldn't update 'current' symlink\n", $log);
$wormbase->run_command("cd $www/WORMBASE; webpublish -q -r current", $log) && $log->write_to("Couldn't run webpublish on current symlink files\n", $log);

# Now need to update big dbcomp output in data directory
$wormbase->run_command("cd /nfs/WWWdev/SANGER_docs/data/Projects/C_elegans; $webpublish -q WS.dbcomp_output", $log) && $log->write_to("Couldn't webpublish data directory\n");


# The end
$log->mail;
exit(0);

