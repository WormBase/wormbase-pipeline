#!/usr/local/bin/perl5.8.0 -w
#
# update_live_release.pl
#
# by Anthony Rogers
#
# Updates the local webpages in synch with the main website
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2010-02-02 10:43:38 $


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
$wormbase->run_command("/software/farm/bin/xdformat -p $wormpub_dir/wormpep_current", $log);
$wormbase->run_command("/software/farm/bin/setdb $wormpub_dir/wormpep_current", $log);
# create a symbolic link from top level wormpep into the release on the wormpep ftp site
foreach my $file ( @wormpep_files ) {
  $wormbase->run_command("cd $wormpep_ftp_root && ln -fs wormpep${release}/${file}${release} $file", $log);
}

# delete the old symbolic link and make the new one to wormpep.prev
my $prev_release = $release -1;
$wormbase->run_command("cd $wormpep_ftp_root && ln -fs wormpep${prev_release}/wormpep${prev_release} wormpep.prev",$log);


##################################################################
# Update Sanger WORMBASE ftp site to change live_release symlink
##################################################################
#Should point at the current WSXXX release folder eg. WS167/, basically one less than the development_release link which will appear above it.
$wormbase->run_command("cd $wormbase_ftp_dir && rm -f live_release && ln -fs WS${release} live_release",$log);




##############################################
# Update pages using webpublish
# Separate webpublish commands (for safety!) on the two top level directories that need updating
##############################################

my $webpublish = "/software/bin/webpublish";

# symlinks are not allowed on the website so we have a "current" directory containing a file for each file we want to link to that is just "include virtual" 
# the place that these link to need to be updated.
my $currentdir = "$www/WORMBASE/current";
opendir(DIR,"$currentdir") or $log->log_and_die("cant opendir $currentdir\n :$!\n");
while(my $file = readdir(DIR) ){
    next unless ($file =~ /shtml/);
    open (OLD,"<$currentdir/$file") or $log->log_and_die("cant opendir $currentdir/$file\n :$!\n");
    open (NEW,">/tmp/$file") or $log->log_and_die("cant opendir /tmp/$file\n :$!\n");
    while(<OLD>){
	s/WS\d+/WS$release/;
	print NEW;
    }
    close OLD;
    close NEW;
    $wormbase->run_command("mv /tmp/$file $currentdir/$file", $log);
    $wormbase->run_command("cd $currentdir && webpublish $file", $log);
}
close DIR;

$wormbase->worm_webpublish("-file" => "$www/WORMBASE/current") or $log->error("Couldn't run webpublish on current symlink files\n");
# The end
$log->mail;
exit(0);

