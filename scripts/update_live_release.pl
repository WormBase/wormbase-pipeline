#!/usr/local/bin/perl5.8.0 -w
#
# update_live_release.pl
#
# by Anthony Rogers
#
# Updates the local webpages in synch with the main website
# Last updated by: $Author: krb $
# Last updated on: $Date: 2004-07-20 09:21:50 $


use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Carp;
use Getopt::Long;

my $release;
my $errors;
GetOptions ("release:s" => \$release );

die "you must give a release version ( just numbers eg -release 125 )\n" unless $release;

my $www = "/nfs/WWWdev/htdocs/Projects/C_elegans";

############################################
# update wormpep files
#############################################

my $wormpub_dir = "/nfs/disk100/wormpub/WORMPEP";
my $wormpep_ftp_root = glob("~ftp/pub/databases/wormpep");
my $wp_ftp_dir = "$wormpep_ftp_root/wormpep${release}";
my @wormpep_files = &wormpep_files;
my $log  = "/wormsrv2/logs/update_live_release.${release}.$$";

open (LOG, ">$log") || die "Couldn't open log file\n";;
print LOG &runtime, " : starting script\n";


# update new live wormpep release from ftp_site to /disk100/wormpub
foreach my $file ( @wormpep_files ) {
  unlink("$wormpub_dir/${file}_current") || print LOG "ERROR: Cannot delete file $wormpub_dir/${file}_current :\t$!\n";
  &run_command("scp $wp_ftp_dir/${file}${release} $wormpub_dir/${file}_current");
}
&run_command("/usr/local/pubseq/bin/setdb $wormpub_dir/wormpep_current");

# create a symbolic link from top level wormpep into the release on the wormpep ftp site
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
&run_command("$webpublish -f -q -r wormpep") && print LOG "Couldn't run webpublish on wormpep files\n";

# Now update WORMBASE current link
chdir("$www/WORMBASE") || print LOG "Couldn't run chdir\n";
&run_command("$webpublish -f -q -r current") && print LOG "Couldn't run webpublish on current symlink files\n";

# Now need to update big dbcomp output in data directory
chdir("/nfs/WWWdev/SANGER_docs/data/Projects/C_elegans") || print LOG "Couldn't chdir to data directory\n";
&run_command("$webpublish -f -q WS.dbcomp_output") && print LOG "Couldn't webpublish data directory\n";


# The end

print LOG &runtime, " : finished\n";
close(LOG);


my $maintainers = "All";
# warn about errors in subject line if there were any
if($errors == 0){
  &mail_maintainer("BUILD REPORT: update_live_release.pl",$maintainers,$log);
}
elsif ($errors ==1){
  &mail_maintainer("BUILD REPORT: update_live_release.pl : $errors ERROR!",$maintainers,$log);
}
else{
  &mail_maintainer("BUILD REPORT: update_live_release.pl : $errors ERRORS!!!",$maintainers,$log);
}
exit(0);


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

