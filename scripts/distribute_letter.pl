#!/usr/local/bin/perl5.6.1 -w
#
# distribute_letter.pl
#
# by Anthony Rogers
#
# copies release letter to ~ftp/pub2/wormbase/WSxx
#                          ~wormpub/BUILD/autoace/release/
#                          /nfs/WWW/SANGER_docs/htdocs/Projects/C_elegans/WORMBASE/current/release_notes.txt/
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-05-15 20:35:15 $


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

my ($help, $debug, $test, $verbose, $store, $wormbase, $mail_dev);
my $maintainers = 'All';



GetOptions ('help'       => \$help,
            'debug=s'    => \$debug,
            'test'       => \$test,
            'verbose'    => \$verbose,
            'store:s'    => \$store,
            'maildev'    => \$mail_dev,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage('Help') if ($help);

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
my $ftpdir = $wormbase->ftp_site;

##############
# variables  #
##############

my $release        = $wormbase->get_wormbase_version_name;
my $release_number = $wormbase->get_wormbase_version;
my $errors         = 0;

$log->write_to("about to spread the word . . . \n");

# copy the letter around
$log->write_to("copying to ftp site . . . . ");

# ftp-site
&_copy( "$repdir/letter.${release}", "$ftpdir/releases/${release}/letter.${release}" ) || $log->log_and_die("couldnt copy $repdir/letter.${release} to $ftpdir/releases/${release}/letter.${release}\n");
$log->write_to("DONE.\n");

# local
$log->write_to("copying to autoace/release . . . . ");
&_copy( "$repdir/letter.${release}", "$acedir/release/letter.${release}" )
    || $log->log_and_die("couldnt copy to autoace/release\n");
$log->write_to("DONE.\n");


##################
# Check the files
##################
$wormbase->check_file("$ftpdir/releases/${release}/letter.${release}", $log,
		      samesize => "$repdir/letter.${release}");
$wormbase->check_file("$acedir/release/letter.${release}", $log,
		      samesize => "$repdir/letter.${release}");


if ($mail_dev) {
# Send email
  print "\n\nMailing to wormbase-dev . .\n";
  
  my $to             = $debug?$maintainers:'dev@wormbase.org';
  my $name           = "Wormbase ${release} release";
  my $release_letter = "$repdir/letter.${release}";
  $wormbase->mail_maintainer( $name, $to, $release_letter);
}

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
