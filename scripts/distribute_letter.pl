#!/usr/local/bin/perl5.6.1 -w
#
# distribute_letter.pl
#
# by Anthony Rogers
#
# copies release letter to ~ftp/pub/wormbase/WSxx
#                          /wormsrv2/autoace/release/
#                          /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/current/release_notes.txt/
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2004-12-03 17:00:28 $


use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;


##############
# variables  #
##############

# Most checking scripts should produce a log file that is a) emailed to us all 
# and b) copied to /wormsrv2/logs

my $maintainers = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
my $release   = &get_wormbase_version_name(); # e.g. WS89
my $release_number = &get_wormbase_version; # e.g. 89
my $log        = "/wormsrv2/logs/distribute_letter.${release}.$$";
my $www = "/nfs/WWWdev/htdocs/Projects/C_elegans";
my $errors = 0; 

open (LOG, ">$log");
print LOG "about to spread the word . . . \n";

# copy the letter around
print LOG "copying to ftp site . . . . ";
my $ftp_dir = glob("~ftp/pub/wormbase");
`cp /wormsrv2/autoace/REPORTS/letter.${release} $ftp_dir/${release}/letter.${release}` and die "couldnt copy to $ftp_dir\n";
print LOG "DONE.\n";
print LOG "copying to autoace/release . . . . ";
`cp /wormsrv2/autoace/REPORTS/letter.${release} /wormsrv2/autoace/release/letter.${release}` and die "couldnt copy to autoace/release";
print LOG "DONE.\n";
print LOG "copying to intranet . . . . ";
`cp /wormsrv2/autoace/REPORTS/letter.${release} ${www}/WORMBASE/${release}/release_notes.txt`
  and die "couldnt copy to ${www}/WORMBASE/${release}/\n";
print LOG "DONE.\n";


# Send email
print "\n\nMailing to wormbase-dev . . ";

my $to = "wormbase-dev\@wormbase.org";
my $name = "Wormbase ${release} release";
my $release_letter = "/wormsrv2/autoace/REPORTS/letter.${release}";
if( &mail_maintainer($name,$to,$release_letter) == 1 ) {
  print LOG "DONE.\n\n\n  Go and have a cuppa !\n";}
else {
  print LOG "! mailing failed !\n\n\nIs this divine intervention?  A last chance to fix something?\n
whatever - something is wrong.\n"; }

	 
################################### 	 
# Make data on FTP site available 	 
################################### 	 
  	 
# FTP site data is there but sym link needs to be updated so people can easily point to it 	 
  	 
print LOG "Updating symlink on FTP site\n"; 	 
  	 
my $targetdir = "/nfs/disk69/ftp/pub/wormbase";  # default directory, can be overidden 	 
  	 
# delete the old symbolic link and make the new one 	 
&run_command("rm -f $targetdir/development_release"); 	 
&run_command("cd $targetdir; ln -s $release development_release"); 	 
 

# update wormpep_dev symbolic link in wormpep ftp site
my $wormpep_dir = glob("~ftp/pub/databases/wormpep");
&run_command("rm -f $wormpep_dir/wormpep_dev");
&run_command("ln -fs $wormpep_dir/wormpep${release_number}/wormpep${release_number} $wormpep_dir/wormpep_dev");


#######################################
# Webpublish to live site
#######################################

print LOG "Updating some WormBase webpages to live site\n"; 	 

# update development_release symbolic link 
chdir("$www/WORMBASE");
&run_command("rm -f development_release");
&run_command("ln -fs $release development_release");

# Now update WORMBASE pages
# these won't be seen until current symlink is also updated
my $webpublish = "/usr/local/bin/webpublish";
&run_command("$webpublish -f -q -r $release") && print LOG "Couldn't run webpublish on release directory\n";
&run_command("$webpublish -f -q -r development_release") && print LOG "Couldn't run webpublish on dev sym link\n";



# say goodnight Brian

print LOG "$0 finished at ",`date`,"\n\n";
close(LOG);


# warn about errors in subject line if there were any
if($errors == 0){
  &mail_maintainer("BUILD REPORT: distribute_letter.pl",$maintainers,$log);
}
elsif ($errors ==1){
  &mail_maintainer("BUILD REPORT: distribute_letter.pl : $errors ERROR!",$maintainers,$log);
}
else{
  &mail_maintainer("BUILD REPORT: distribute_letter.pl : $errors ERRORS!!!",$maintainers,$log);
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
    print LOG "ERROR: $command failed\n";
    $errors++;
  }

  # for optional further testing by calling subroutine
  return($status);
}



########################################################



__END__

=pod

=head2 NAME - distribute_letter.pl

=head1 USAGE

=over 4

=item distribute_letter.pl

=back

This script:

copies the release letter to the ftp site, website and autoace/release

mails release letter to wormbase-dev

Then, when release letter is emailed it updates the symlink on the FTP
site to make current_release point to the latest release directory.

Finally, it runs webpublish in one place to update the dev site with the live site.
This is just so that WashU will be able to see the latest database checks
statistics even when they are not yet live.

script_template.pl MANDATORY arguments:

=over 4

=item none

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
