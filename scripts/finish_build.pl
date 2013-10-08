#!/usr/bin/env perl
#
# finish_build.pl
#
# Usage : finish_build.pl [-options]
#
# 1) Archives release CURRENT - 3
# 2) Moves selected sub-folders from autoace to elegans folder
# 3) Moves remainder of autoace to DATABASES
# 4) Cleans out old(er) dbs in DATABASES, to keep the bare minimum
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-10-08 14:45:07 $


use strict;
use lib $ENV{CVS_DIR};
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Cwd;
use File::Glob ':glob';
use File::Path;
use Carp;
use Log_files;
use Storable;

my $ARCHIVE_DIR = '/warehouse/wormbase01/wormarchive';

my @SPECIES_SPECIFIC_FOLDERS = qw(
 blat
 gff_splits
 spell
 misc_output
 acefiles
 transcripts
    );

my @DO_NOT_KEEP_LONG_TERM = qw(
  TMP
  CHECKS
  logs
);

##############################
# command-line options       #
##############################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my $noarchive;

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "noarchive:s"=> \$noarchive
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
$log->log_and_die("FTP_LOCK present indicating that make_FTP_sites.pl failed\nfix this before continuing\n\n$lockfile\n") if -e ($lockfile);

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

#############################################################################
# Step 0: archive and delete the oldest release
#############################################################################
my $WS_oldest   = $WS_current - 3; # the version that *should* be the oldest
my $WS_old_name = "WS".$WS_oldest;
my $WS_old_path = $wormbase->database("$WS_old_name");
&archive_old_release ($WS_old_name) unless (defined $noarchive);

#############################################################################
# Step 1: move elegans-specific working folders into the "elegans" build dir
#############################################################################
my $elegans_dir = $wormbase->basedir . "/elegans";

foreach my $meth (@SPECIES_SPECIFIC_FOLDERS) {
  my $source = $wormbase->$meth;
  my $target = "$elegans_dir/";
  $wormbase->run_command("mv -f $source $target", $log);
}

#############################################################################
# Step 2: move the remaining build dir to the DATABASES directory
#############################################################################
my $new_dir = $wormbase->wormpub."/DATABASES/".$wormbase->get_wormbase_version_name;
warn "mv ".$wormbase->autoace." $new_dir";
$wormbase->run_command("mv ".$wormbase->autoace." $new_dir", $log);

#############################################################################
# Step 3: remove redundant files that we do not need any more
#############################################################################
$log->write_to("Removing files in $new_dir/database/new/\n");
$wormbase->delete_files_from("$new_dir/database/new","*","+") if (-e "$new_dir/database/new");

$log->write_to("Removing files in $new_dir/database/touched/\n");
$wormbase->delete_files_from("$new_dir/database/touched","*","+") if (-e "$new_dir/database/touched");

#remove large backups created by loading routine
$log->write_to("Removing files in $new_dir/database/backup*\n");
$wormbase->delete_files_from("$new_dir/database/","backup","-");

#############################################################################
# Step 4: update the current_db symlink
#############################################################################
$wormbase->run_command("rm ".$wormbase->database('current'), $log);
$wormbase->run_command("ln -sf ".$new_dir." ".$wormbase->database('current'),$log);

#############################################################################
# Step 5: clean up stuff from older directories
#############################################################################
my $old_version = $WS_current -1;
my $old_dir = $wormbase->wormpub."/DATABASES/WS".$old_version;

foreach my $folder (@DO_NOT_KEEP_LONG_TERM) {
  $log->write_to("Removing $old_dir/$folder\n");
  $wormbase->delete_files_from("$old_dir/$folder","*","+") unless ($test);
}

$log->write_to("Zipping files in CHROMOSOMES and SEQUENCES folders\n");
foreach my $file (glob("$old_dir/SEQUENCES/*.*"), glob("$old_dir/CHROMOSOMES/*.*")) {
  next if $file =~ /\.gz$/;
  $wormbase->run_command("gzip $file", $log);
}


##################
# End
##################

$log->write_to("Finished.\n");
$log->mail();
print "Finished.\n" if ($verbose);
exit(0);

########################################################
sub archive_old_release {
  my $WS_old_name = shift;
  my $WS_old_path = $wormbase->database("$WS_old_name");


  unless ($WS_old_path) {
    $log->write_to("cant find database for $WS_old_path - not archiving\n");
    $log->error;
    return 1;
  }
  my $tar = "${ARCHIVE_DIR}/${WS_old_name}.tar.gz";
  $log->write_to("backing up $WS_old_path to $tar\n");

  # turn the old release into a tarball, move into $archive_dir and remove old directory
  $log->write_to("\nCreating $tar\n");

  my $tar_cmd = "tar -C $WS_old_path/.. -cvzf $tar $WS_old_path";
  $wormbase->run_command($tar_cmd, $log) and 
      $log->log_and_die("ERROR in tar cmd: $tar_cmd\n");
  my $minsize = 7500000000; # 7 Gb
  if ($wormbase->check_file($tar, $log,
			    minsize => $minsize,
			   )) {
    $log->log_and_die("ERROR copying $WS_old_path to $tar - the result is less than $minsize\n");
  } else {
    $wormbase->run_command("rm -rf $WS_old_path", $log) && 
        $log->log_and_die("ERROR in rm -rf $WS_old_path\n");
  }
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
