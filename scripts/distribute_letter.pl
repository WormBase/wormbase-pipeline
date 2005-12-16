#!/usr/local/bin/perl5.6.1 -w
#
# distribute_letter.pl
#
# by Anthony Rogers
#
# copies release letter to ~ftp/pub/wormbase/WSxx
#                          /wormsrv2/autoace/release/
#                          /nfs/WWW/SANGER_docs/htdocs/Projects/C_elegans/WORMBASE/current/release_notes.txt/
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2005-12-16 11:18:55 $

use strict;
use warnings;
use File::Copy;
use Getopt::Long;
use lib "/wormsrv2/scripts/";
use Wormbase;

#####################################
# variables and command-line options #
######################################
my $maintainers = 'All';

my $test;       # In test build mode
my $help;       # Help/Usage page
my $debug;      # Debug mode
my $store;      # configuration file

GetOptions(
    'help'    => \$help,
    'test'    => \$test,
    'debug=s' => \$debug,
    'store=s' => \$store
);

# Help pod if needed
&usage(0) if ($help);

############################
# recreate configuration   #
############################
my $wb;
if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug, -test => $test, ) }

###########################################
# Variables Part II (depending on $wb)    #
###########################################
$test  = $wb->test  if $wb->test;     # Test mode
$debug = $wb->debug if $wb->debug;    # Debug mode, output only goes to one user
my $repdir = $wb->reports;            # Reports path
my $acedir = $wb->autoace;             # base directory of the build

# Use debug mode?
if ($debug) {
    print "DEBUG = \"$debug\"\n\n";
    ( $maintainers = $debug . '\@sanger.ac.uk' );
}

# create log
my $log = Log_files->make_build_log($debug);

##############
# variables  #
##############

# Most checking scripts should produce a log file that is a) emailed to us all
# and b) copied to /wormsrv2/logs

<<<<<<< distribute_letter.pl
my $maintainers = "All";
my $rundate     = `date +%y%m%d`; chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;
my $release   = &get_wormbase_version_name(); # e.g. WS89
my $release_number = &get_wormbase_version; # e.g. 89
my $log        = "/wormsrv2/logs/distribute_letter.${release}.$$";
my $www = "/nfs/WWWdev/SANGER_docs/htdocs/Projects/C_elegans";
my $errors = 0; 
=======
my $rundate        = $wb->rundate;
my $runtime        = $wb->runtime;
my $release        = $wb->get_wormbase_version_name;
my $release_number = $wb->get_wormbase_version;
my $www            = "/nfs/WWWdev/SANGER_docs/htdocs/Projects/C_elegans";
my $errors         = 0;
>>>>>>> 1.19.4.4

$log->write_to("about to spread the word . . . \n");

# copy the letter around
$log->write_to("copying to ftp site . . . . ");
my ($ftp_dir) = glob("~ftp/pub/wormbase");       
      # ftp-site
      &_copy( "$repdir/letter.${release}", "$ftp_dir/${release}/letter.${release}" ) || die "couldnt copy to $ftp_dir\n";
      $log->write_to("DONE.\n");

      # local
      $log->write_to("copying to autoace/release . . . . ");
      &_copy( "$repdir/REPORTS/letter.${release}", "$acedir/release/letter.${release}" )
      || die "couldnt copy to autoace/release\n";
      $log->write_to("DONE.\n");

      # web fluff
      $log->write_to("copying to intranet . . . . ");
      &_copy( "$repdir/letter.${release}", "${www}/WORMBASE/${release}/release_notes.txt" )
      || die "couldnt copy to ${www}/WORMBASE/${release}/\n";
      $log->write_to("DONE.\n");

  # Send email
  print "\n\nMailing to wormbase-dev . .\n";

my $to             = $debug?$maintainers:"wormbase-dev\@wormbase.org";
my $name           = "Wormbase ${release} release";
my $release_letter = "$repdir/letter.${release}";
$wb->mail_maintainer( $name, $to, $release_letter);

###################################
# Make data on FTP site available
###################################

# FTP site data is there but sym link needs to be updated so people can easily point to it
$log->write_to("Updating symlink on FTP site\n");

my $targetdir = "/nfs/disk69/ftp/pub/wormbase";    # default directory, can be overidden

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

$log->write_to("Updating some WormBase webpages to live site\n");

<<<<<<< distribute_letter.pl
# update development_release symbolic link 
chdir("$www/WORMBASE");
&run_command("rm -f development_release");
&run_command("ln -fs $release development_release");
=======
# update development_release symbolic link
chdir("$www/WORMBASE");
&run_command("rm -f development_release");
&run_command("ln -fs $release development_release");
>>>>>>> 1.19.4.4

# Now update WORMBASE pages
# these won't be seen until current symlink is also updated
my $webpublish = "/usr/local/bin/webpublish";
<<<<<<< distribute_letter.pl
&run_command("$webpublish -f -q -r $release") && print LOG "Couldn't run webpublish on release directory\n";
&run_command("$webpublish -f -q -r development_release") && print LOG "Couldn't run webpublish on dev sym link\n";


=======
&run_command("$webpublish -f -q -r $release")            && $log->write_to("Couldn't run webpublish on release directory\n");
&run_command("$webpublish -f -q -r development_release") && $log->write_to("Couldn't run webpublish on dev sym link\n");
>>>>>>> 1.19.4.4

# say goodnight Brian
$log->write_to("$0 finished at ", `date`, "\n\n");

# warn about errors 
$log->mail( "$maintainers", "BUILD REPORT: $0" );

exit(0);

##################################################################################
#
# Simple routine which will run commands via system calls but also check the
# return status of a system call and complain if non-zero, increments error check
# count, and prints a log file error
#
##################################################################################

sub _copy {
	my ($from,$to)=@_;
	if ($debug||$test) {
		print $from,' => ',$to,"\n";
		return 1;
	}
	else {
		return copy($from,$to);
	}
}

sub run_command {
    my $command = shift;
    $log->write_to($wb->runtime.": started running $command\n");
    my $status = 0;
    if ( $test||$debug ) { print $command,"\n\n" }
    else {
        $status = system($command)
    }
    if ( $status != 0 ) {
        $log->write_to("ERROR: $command failed\n");
        $errors++;
    }

    # for optional further testing by calling subroutine
    return ($status);
}

sub usage {
    my $error = shift;

    if ( $error eq "Help" ) {
        # Normal help menu
        system( 'perldoc', $0 );
        exit(0);
    }
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

=item -d username, debug username ,will only print stuff and send mails to the specified user

=item -t, test , basically the same as debug, but will use the default email addresses

=item -s filename, specify a configuration file

=back

=back

=head1 AUTHOR

=over 4

=item Anthony Rogers (ar2@sanger.ac.uk)

=back

=cut
