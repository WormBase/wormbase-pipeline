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
# Last updated by: $Author: krb $                       
# Last updated on: $Date: 2003-04-07 09:12:34 $         


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
my $ver         = &get_wormbase_version;


open (LOG, ">$log");
print LOG "about to spread the word . . . \n";

# copy the letter around
print LOG "copying to ftp site . . . . ";
my $ftp_dir = glob("~ftp/pub/wormbase");
`cp /wormsrv2/autoace/RELEASE_LETTERS/letter.WS$ver $ftp_dir/WS$ver/letter.WS$ver` and die "couldnt copy to $ftp_dir\n";
print LOG "DONE.\n";
print LOG "copying to autoace/release . . . . ";
`cp /wormsrv2/autoace/RELEASE_LETTERS/letter.WS$ver /wormsrv2/autoace/release/letter.WS$ver` and die "couldnt copy to autoace/release";
print LOG "DONE.\n";
print LOG "copying to intranet . . . . ";
`cp /wormsrv2/autoace/RELEASE_LETTERS/letter.WS$ver /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/current/release_notes.txt`
  and die "couldnt copy to /nfs/WWW/htdocs/Projects/C_elegans/WORMBASE/current/\n";
print LOG "DONE.\n";

print "\n\nMailing to wormbase-dev . . ";

my $to = "wormbase-dev\@wormbase.org";
my $name = "Wormbase WS$ver release";
my $release_letter = "/wormsrv2/autoace/RELEASE_LETTERS/letter.WS$ver";
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
my $release   = &get_wormbase_version_name(); # e.g. WS89

# delete the old symbolic link and make the new one
system "rm -f $targetdir/current_release";
system "cd $targetdir; ln -s $release current_release";


# now update pages using webpublish
# Separate webpublish commands (for safety!) on the two top level directories that need updating
my $www = "/nfs/WWWdev/htdocs/Projects/C_elegans";
chdir($www) || print LOG "Couldn't run chdir\n";

system("/usr/local/bin/webpublish -f -q -r wormpep") && print LOG "Couldn't run webpublish on wormpep files\n";

chdir("$www/WORMBASE") || print LOG "Couldn't run chdir\n";
system("/usr/local/bin/webpublish -f -q -r $release") && print LOG "Couldn't run webpublish on release directory\n";
system("/usr/local/bin/webpublish -f -q -r current") && print LOG "Couldn't run webpublish on current symlink files\n";

# Now need to unpack dbcomp file on live site
chdir("$www/WORMBASE/${release}") || print LOG "Couldn't chdir to $www/WORMBASE/${release}\n";
system("/bin/gunzip WS.dbcomp_output.gz") && print LOG "Couldn't gunzip WS.dbcomp_output.gz\n";
system("rm -f WS.dbcomp_output.gz") && print LOG "Couldn't remove gzip fil\n"; 

print LOG "$0 finished at ",`date`,"\n\n";


# Always good to cleanly exit from your script
exit(0);


# Add perl documentation in POD format
# This should expand on your brief description above and add details of any options
# that can be used with the program.  Such documentation can be viewed using the perldoc
# command.


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

Finally, when release letter is emailed it updates the symlink on the FTP
site to make current_release point to the latest release directory.

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
