#!/usr/local/bin/perl5.8.0 -w
#
# update_live_release.pl
#
# by Anthony Rogers
#
# Updates the local webpages in synch with the main website
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2004-05-05 10:28:40 $


use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Carp;



my $release   = &get_wormbase_version_name(); # e.g. WS89
my $www = "/nfs/WWWdev/htdocs/Projects/C_elegans";
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



###########################################################################################################
# Create new symbolic link
###########################################################################################################
my $wwwlive          = "/nfs/WWW/htdocs/Projects/C_elegans/WORMBASE";
print LOG "\nChanging 'current symbolic link to point to new release\n";
system("rm -f $www/current") && croak "Couldn't remove 'current' symlink\n";
system("ln -s $wwwlive/$release/ $www/current") && croak "Couldn't create new symlink\n";

