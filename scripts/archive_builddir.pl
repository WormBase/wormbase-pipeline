#!/usr/bin/env perl
#
# archive_build.pl
#
# Archives the given build dir to archive storage
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-10-11 14:42:29 $


use strict;
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;
use File::Basename;

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($sourcedir, $targetdir, $target_fname, $delete, @exclude);

$targetdir = "/warehouse/wormbase01/wormarchive";

GetOptions ("help"         => \$help,
            "debug=s"      => \$debug,
	    "test"         => \$test,
	    "verbose"      => \$verbose,
	    "store:s"      => \$store,
            "sourcedir=s"     => \$sourcedir,
            "targetdir=s"     => \$targetdir,
            "targetfname=s"   => \$target_fname,
            "delete"       => \$delete,
            "exclude=s@"   => \@exclude,
	    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

my $log = Log_files->make_build_log($wormbase);

if (not $sourcedir) {
  $log->log_and_die("You must supply a source build/release dir with -source\n");
}

$target_fname = basename($sourcedir) . ".tar.gz" if not defined $target_fname;
$log->write_to("Will archive $sourcedir to $targetdir/$target_fname\n");

@exclude = map { "--exclude=$sourcedir/$_" } @exclude;
my $tar_cmd = "tar -cvzf $targetdir/$target_fname -C $sourcedir/.. @exclude $sourcedir";
$wormbase->run_command($tar_cmd, $log) 
    and $log->log_and_die("ERROR: failed to produce tar archive for $sourcedir\n");

if ($delete) {
  $wormbase->run_command("rm -fr $sourcedir", $log)
}

$log->write_to("Finished\n");
$log->mail();
exit(0);

