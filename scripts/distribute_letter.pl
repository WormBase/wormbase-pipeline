#!/usr/local/bin/perl5.6.1 -w
#
# distribute_letter.pl
#
# by Anthony Rogers
#
# copies release letter to ~ftp/pub/wormbase/WSxx
#                          ~wormpub/BUILD/autoace/release/
#                          /nfs/WWW/SANGER_docs/htdocs/Projects/C_elegans/WORMBASE/current/release_notes.txt/
#
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2006-05-17 08:57:14 $


use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use warnings;
use File::Copy;



######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my $maintainers = 'All';



GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"      => \$store,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);
  $maintainers = `whoami` . '\@sanger.ac.uk';
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

# Use debug mode?
if ($debug) {
  print "DEBUG = \"$debug\"\n\n";
  $maintainers = $debug . '\@sanger.ac.uk';
}


################
# Variables    #
################

my $repdir = $wormbase->reports;            # Reports path
my $acedir = $wormbase->autoace;             # base directory of the build

##############
# variables  #
##############

my $release        = $wormbase->get_wormbase_version_name;
my $release_number = $wormbase->get_wormbase_version;
my $www            = "/nfs/WWWdev/SANGER_docs/htdocs/Projects/C_elegans";
my $errors         = 0;

$log->write_to("about to spread the word . . . \n");

# copy the letter around
$log->write_to("copying to ftp site . . . . ");
my ($ftp_dir) = glob("~ftp/pub/wormbase");       
      # ftp-site
      &_copy( "$repdir/letter.${release}", "$ftp_dir/${release}/letter.${release}" ) || die "couldnt copy to $ftp_dir\n";
      $log->write_to("DONE.\n");

      # local
      $log->write_to("copying to autoace/release . . . . ");
      &_copy( "$repdir/letter.${release}", "$acedir/release/letter.${release}" )
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
$wormbase->mail_maintainer( $name, $to, $release_letter);

###################################
# Make data on FTP site available
###################################

# FTP site data is there but sym link needs to be updated so people can easily point to it
$log->write_to("Updating symlink on FTP site\n");

my $targetdir = "/nfs/disk69/ftp/pub/wormbase";    # default directory, can be overidden

# delete the old symbolic link and make the new one
$wormbase->run_command("rm -f $targetdir/development_release", $log);
$wormbase->run_command("cd $targetdir; ln -s $release development_release", $log);

# update wormpep_dev symbolic link in wormpep ftp site
my $wormpep_dir = glob("~ftp/pub/databases/wormpep"); 
$wormbase->run_command("rm -f $wormpep_dir/wormpep_dev", $log);
$wormbase->run_command("ln -fs $wormpep_dir/wormpep${release_number}/wormpep${release_number} $wormpep_dir/wormpep_dev", $log);

#######################################
# Webpublish to live site
#######################################

$log->write_to("Updating some WormBase webpages to live site\n");

# update development_release symbolic link
chdir("$www/WORMBASE");
$wormbase->run_command("rm -f development_release", $log);
$wormbase->run_command("ln -fs $release development_release", $log);

# Now update WORMBASE pages
# these won't be seen until current symlink is also updated
my $webpublish = "/usr/local/bin/webpublish";
$wormbase->run_command("$webpublish  -q -r $release", $log)            && $log->write_to("Couldn't run webpublish on release directory\n");
$wormbase->run_command("$webpublish  -q -r development_release", $log) && $log->write_to("Couldn't run webpublish on dev sym link\n");


$log->mail();
print "Finished.\n" if ($verbose);
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

##########################################


sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################



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
