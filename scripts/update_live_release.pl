#!/usr/local/bin/perl5.8.0 -w
#
# update_live_release.pl
#
# by Anthony Rogers
#
# Updates the local webpages in synch with the main website
# Last updated by: $Author: krb $
# Last updated on: $Date: 2004-06-18 15:22:07 $


use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Carp;
use Getopt::Long;

my $release;
my $errors;
Getoptions ("release:s" => \$release );

die "you must give a release version ( just numbers eg -release 125 )\n" unless $release;

my $www = "/nfs/WWWdev/htdocs/Projects/C_elegans";

############################################
# update wormpep files
#############################################

my $wormpub_dir = "/nfs/disk100/wormpub/WORMPEP";
my $wormpep_ftp_root = glob("~ftp/pub/databases/wormpep");
my $wp_ftp_dir = "$wormpep_ftp_root/wormpep${release}";
my @wormpep_files = &wormpep_file;

# update new live wormpep release from ftp_site to /disk100/wormpub
foreach my $file ( @wormpep_files ) {
  unlink("$wormpub_dir/${file}_current") || print LOG "ERROR: Cannot delete file $wormpub_dir/${file}_current :\t$!\n";
  &run_command("scp $wp_ftp_dir/$file $wormpub_dir/${file}_current");
}
&run_command("/usr/local/pubseq/bin/setdb $wormpub_dir/wormpep_current");

# create a symbolic link from top level wormpep into the release on the live website
foreach my $file ( @wormpep_files ) {
  &run_command("cd $wormpep_ftp_root; ln -fs $wp_ftp_dir/$file $file");
}

# delete the old symbolic link and make the new one to wormpep.prev
my $prev_release = $release -1;
&run_command("cd $wormpep_ftp_root; ln -fs wormpep/wormpep${prev_release} wormpep.prev");




##############################################
# Update pages using webpublish
# Separate webpublish commands (for safety!) on the two top level directories that need updating
##############################################

my $webpublish = "/usr/local/bin/webpublish";

#First update wormpep subdirectory
chdir($www) || print LOG "Couldn't run chdir\n";
system("$webpublish -f -q -r wormpep") && print LOG "Couldn't run webpublish on wormpep files\n";

# Now update WORMBASE pages
chdir("$www/WORMBASE") || print LOG "Couldn't run chdir\n";
system("$webpublish -f -q -r $release") && print LOG "Couldn't run webpublish on release directory\n";
system("$webpublish -f -q -r current") && print LOG "Couldn't run webpublish on current symlink files\n";

# Now need to update big dbcomp output in data directory
chdir("/nfs/WWWdev/SANGER_docs/data/Projects/C_elegans") || print LOG "Couldn't chdir to data directory\n";
system("$webpublish -f -q WS.dbcomp_output") && print LOG "Couldn't webpublish data directory\n";


##################################################################################
#
# Simple routine which will run commands via system calls but also check the 
# return status of a system call and complain if non-zero, increments error check 
# count, and prints a log file error
#
##################################################################################

sub run_command{
  my $command = shift;
  print LOG &runtime, ": started running $command\n";
  my $status = system($command);
  if($status != 0){
    $errors++;
    print LOG "ERROR: $command failed\n";
  }

  # for optional further testing by calling subroutine
  return($status);
}

