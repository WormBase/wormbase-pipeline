#!/usr/local/bin/perl5.6.0 -w
#
# archive_dbs.pl
# 
# by Keith Bradnam aged 12 and a half
#
# Usage : archive_dbs.pl [-options]
#
# A simple script to be run at the start of the weekly rebuild.
# Checks to see if there are three existing (and unpacked) WS releases in /wormsrv2
# If there are, then it archives the oldest release away into /wormsrv2/wormbase_archive
# Does a similar thing with Wormpep releases in /wormsrv2/WORMPEP

$| = 1;
use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use IO::Handle;

#################################################################################
# variables                                                                     #
#################################################################################

my $maintainers = "dl1\@sanger.ac.uk krb\@sanger.ac.uk kj2\@sanger.ac.uk";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
our $log        = "/wormsrv2/logs/archive_dbs.$rundate";
my $cvs_version = &get_cvs_version("$0");


my $WS_name     = &get_wormbase_version_name;
my $WS_current  = &get_wormbase_version;
my $WS_oldest   = $WS_current - 3; # the version that *should* be the oldest in /wormsrv2
my $WS_old_name = "WS".$WS_oldest;
my $db_path     = "/wormsrv2";
my $WS_old      = "$db_path"."/$WS_old_name";


&create_log_file;
&delete_autoace_files;

# turn the old release into a tarball, move into wormbase_archive and remove old directory
print LOG "\nCreating $WS_old.tar.gz\n";
system ("tar -cvf $WS_old.tar $WS_old/") && die "Couldn't create tar file\n";
system ("gzip $WS_old.tar") && die "Couldn't create gzip file\n";
print LOG "Moving archive to wormbase_archive and removing $WS_old\n";
system ("mv $WS_old.tar.gz $db_path/wormbase_archive/") && die "Couldn't move to wormbase_archive\n";
system ("rm -rf $WS_old") && die "Couldn't remove old directory\n";


# now do the same thing to archive the old wormpep version
my $wormpep_ver = &get_wormpep_version;
my $old_wormpep = "$db_path/WORMPEP/wormpep".($wormpep_ver-3);

if(-d $old_wormpep){
  print LOG "\nCreating $old_wormpep archive\n";
  system ("tar -cvf $old_wormpep.tar $old_wormpep") && die "Couldn't create wormpep tar file\n";
  system ("gzip $old_wormpep.tar") && die "Couldn't create gzip wormpep tar file\n";
  system ("rm -rf $old_wormpep") && die "Couldn't remove old directory\n";
}

print LOG "C'est finis.\n";
close(LOG);

&mail_maintainer("WormBase Report: archive_dbs.pl",$maintainers,$log);


exit(0);







#################################################################################
# set up log file                                                               #
#################################################################################

sub create_log_file{

  open (LOG,">$log") || die "Cannot open logfile $!\n";
  LOG->autoflush();
  
  print LOG "# update_website.pl\n\n";     
  print LOG "# run details    : $rundate $runtime\n";
  print LOG "# WormBase version : $WS_name\n";
  print LOG "# cvs version      : $cvs_version\n";
  print LOG "\n\n";

}


#################################################################################
# remove non-essential files from old database directory                        #
#################################################################################

sub delete_autoace_files{
  print LOG "\n\n";

  if (-d "$WS_old"){
    
    if (-d "$WS_old/database"){
      print LOG "Removing $WS_old/database\n";
      system("rm -rf $WS_old/database") && die "Couldn't remove $WS_old/database\n";
    }
    
    # remove wspec, wgf etc.
    my @files = glob("$WS_old/w*");
    foreach my $file (@files){
      print LOG "Removing $WS_old/$file\n";
      system("rm -rf $file") && die "Couldn't remove $file\n";
    }
    
    if (-d "$WS_old/pictures"){
      print LOG "Removing $WS_old/pictures\n";
      system("rm -rf $WS_old/pictures") && die "Couldn't remove $WS_old/pictures\n";
    }
    
  }
}
