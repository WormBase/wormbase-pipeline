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
# Last updated by: $Author: dl1 $                       
# Last updated on: $Date: 2004-04-26 08:36:52 $         


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
our $log        = "/wormsrv2/logs/distribute_letter.$rundate";
my $release   = &get_wormbase_version_name(); # e.g. WS89
my $www = "/nfs/WWWdev/htdocs/Projects/C_elegans";

open (LOG, ">$log");
print LOG "about to spread the word . . . \n";

# copy the letter around
print LOG "copying to ftp site . . . . ";
my $ftp_dir = glob("~ftp/pub/wormbase");
`cp /wormsrv2/autoace/RELEASE_LETTERS/letter.${release} $ftp_dir/${release}/letter.${release}` and die "couldnt copy to $ftp_dir\n";
print LOG "DONE.\n";
print LOG "copying to autoace/release . . . . ";
`cp /wormsrv2/autoace/RELEASE_LETTERS/letter.${release} /wormsrv2/autoace/release/letter.${release}` and die "couldnt copy to autoace/release";
print LOG "DONE.\n";
print LOG "copying to intranet . . . . ";
`cp /wormsrv2/autoace/RELEASE_LETTERS/letter.${release} ${www}/WORMBASE/${release}/release_notes.txt`
  and die "couldnt copy to ${www}/WORMBASE/${release}/\n";
print LOG "DONE.\n";

print "\n\nMailing to wormbase-dev . . ";

my $to = "wormbase-dev\@wormbase.org";
my $name = "Wormbase ${release} release";
my $release_letter = "/wormsrv2/autoace/RELEASE_LETTERS/letter.${release}";
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
system "rm -f $targetdir/development_release";
system "cd $targetdir; ln -s $release development_release";

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

print LOG "$0 finished at ",`date`,"\n\n";

exit(0);




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

Finally, it runs webpublish in a few separate places to update the dev
site with the live site.

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
