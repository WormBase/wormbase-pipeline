#!/usr/local/bin/perl5.8.0 -w
#
# finish_build.pl
# 
# by Keith Bradnam
#
# Usage : finish_build.pl [-options]
#
# 1) checks to see if there are three existing (and unpacked) WS releases in /wormsrv2
#    If there are, then it archives the oldest release away into /nfs/wormarchive
# 2) Does a similar thing with Wormpep releases in /wormsrv2/WORMPEP
# 3) Archives old GFF_SPLITS directory
# 4) Makes current_DB (copy of latest release) in ~wormpub/DATABASES
#
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2006-12-01 10:45:26 $


use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Cwd;
use File::Glob ':glob';
use File::Path;
use Carp;
use Log_files;
use Storable;

##############################
# command-line options       #
##############################

my ($help, $debug, $test, $verbose, $store, $wormbase);


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
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


my $lockfile = $wormbase->autoace."/FTP_LOCK";
$log->log_and_die("FTP_LOCK present indicating that make_FTP_sites.pl failed/nfix this before continuing/n/n$lockfile\n") if -e ($lockfile);

#################################################################################
# variables                                                                     #
#################################################################################

# the following point to various BUILD directories (or TEST_BUILD if $test)
my $basedir     = $wormbase->basedir;
my $ace_dir     = $wormbase->autoace; # AUTOACE DATABASE DIR
my $chromosomes = $wormbase->chromosomes; # AUTOACE CHROMSOMES
my $logs        = $wormbase->logs; # AUTOACE LOGS

# need to change this if in test mode
my $WS_name;
my $WS_current;

$WS_current = $wormbase->get_wormbase_version;
$WS_name    = $wormbase->get_wormbase_version_name;


my $WS_new      = $WS_current + 1;
my $WS_new_name = "WS".$WS_new;
my $WS_oldest   = $WS_current - 3; # the version that *should* be the oldest
my $WS_old_name = "WS".$WS_oldest;
my $WS_old_path = $wormbase->database("$WS_old_name");
my $old_wormpep = "$basedir/WORMPEP/wormpep".($WS_current-3);

#####################################################################################

&archive_old_releases ($WS_old_name);

my $new_dir = $wormbase->wormpub."/DATABASES/".$wormbase->get_wormbase_version_name;

# move autoace to DATABASES/WSXXX
$wormbase->run_command("mv ".$wormbase->autoace ." ". $new_dir, $log);

# update symlink current_DB to new build
$wormbase->run_command("unlink ".$wormbase->database('current'),$log);
$wormbase->run_command("ln -sf ".$new_dir." ".$wormbase->database('current'),$log);


# Transfer autoace to ~wormpub/DATABASES/current_DB - first remove existing files
$log->write_to("Removing $new_dir/acefiles/\n");
$wormbase->delete_files_from("$new_dir/acefiles","*","+") unless ($test);

$log->write_to("Unzipping any gzipped chromosome files\n");
if (!$test) {
  $wormbase->run_command("/bin/gunzip $new_dir/CHROMOSOMES/*.gz", $log);
}

# Remove redundant files and directories in $ace_dir
$log->write_to("Removing old files in $ace_dir/release/\n");
$wormbase->delete_files_from("$new_dir/release","*","-") if (-e "$new_dir/release");

$log->write_to("Removing files in $new_dir/database/new/\n");
$wormbase->delete_files_from("$new_dir/database/new","*","+") if (-e "$new_dir/database/new");

$log->write_to("Removing files in $new_dir/database/touched/\n");
$wormbase->delete_files_from("$new_dir/database/touched","*","+") if (-e "$new_dir/database/touched");

#remove large backups created by loading routine
$log->write_to("Removing files in $new_dir/database/backup*\n");
$wormbase->delete_files_from("$new_dir/database/","backup","-");

##################
# End
##################

$log->write_to("Finished.\n");
$log->mail();
print "Finished.\n" if ($verbose);
exit(0);




#################################################################################
# remove non-essential files from old database directory                        #
#################################################################################



sub archive_old_releases{
  my $WS_old_name = shift;
  my $WS_old_path = $wormbase->database("$WS_old_name");
  unless ($WS_old_path) {
    $log->write_to("cant find database for $WS_old_path - not archiving\n");
    $log->error;
    return 1;
  }
  $log->write_to("backing up $WS_old_path to /nfs/wormarchive/$WS_old_name.tar.gz\n");

  # turn the old release into a tarball, move into /nfs/wormarchive and remove old directory
  $log->write_to("\nCreating $WS_old_name.tar.gz\n");

  my $tar = $wormbase->database("$WS_old_name").".tar";
  $wormbase->run_command("tar -cvf $tar ".$wormbase->database("$WS_old_name"), $log) && 
      $log->log_and_die("ERROR in tar -cvf $tar\n");
  $wormbase->run_command("gzip $tar", $log) && 
      $log->log_and_die("ERROR in gzip $tar\n");
  $wormbase->run_command("mv $tar.gz /nfs/wormarchive/", $log) && 
      $log->log_and_die("ERROR in mv $tar.gz /nfs/wormarchive/\n");
  $wormbase->run_command("rm -rf ".$wormbase->database("$WS_old_name"), $log) && 
      $log->log_and_die("ERROR in rm -rf " . $wormbase->database("$WS_old_name") . "\n");
}

#################################################################################
# Prints help and disappears                                                    #
#################################################################################

sub usage {
    exec ('perldoc',$0);
}


__END__

=pod

=head2 NAME - finish_build.pl

=head1 USAGE

=over 4

=item finish_build.pl  [-options]

=back

This script:

 1) checks to see if there are three existing (and unpacked) WS releases 
 in ~wormpub/BUILD. If there are, then it archives the oldest release away into 
 ~wormpub/BUILD/wormbase_archive
 2) Does a similar thing with Wormpep releases in ~wormpub/BUILD/WORMPEP
 but 
 3) Runs GFFsplitter.pl -a to archive away the last GFF_SPLITS directory
 4) Copies autoace into a separate WSxx directory
 5) updates the ~wormpub/BUILD/current_DB symlink to point to the directory created
    in step 4.

finish_build.pl MANDATORY arguments:

=over 4

=item none

=back

finish_build.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Keith Bradnam (krb@sanger.ac.uk)

=back

=cut
