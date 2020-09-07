#!/usr/bin/env perl
#
# finish_build.pl
#
# Usage : finish_build.pl [-options]
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-10-14 09:33:14 $


use strict;
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

#
# The following folders are treated as "scratch" folders
# that are not really part of the build output. We move them into
# the "elegans" build dir where they will hang around until 
# the merge point in the next build
#
my @BUILD_SCRATCH_FOLDERS = qw(
 acefiles
 BLAT
 CHECKS
 GFF_SPLITS
 TMP
 TRANSCRIPTS
    );

##############################
# command-line options       #
##############################

my ($help, $debug, $test, $verbose, $store, $wormbase);

GetOptions ("help"         => \$help,
            "debug=s"      => \$debug,
	    "test"         => \$test,
	    "verbose"      => \$verbose,
	    "store:s"      => \$store,
	    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new()
}

# the following point to various BUILD directories (or TEST_BUILD if $test)
my $basedir     = $wormbase->basedir;
my $autoace     = "$basedir/autoace";
my $WS_current = $wormbase->get_wormbase_version;
my $WS_name    = $wormbase->get_wormbase_version_name;

#############################################################################
# Step 1: move elegans-specific working folders into the "elegans" build dir
#############################################################################
my $elegans_dir = $wormbase->basedir . "/elegans";

print STDERR "Moving Build scratch folders into elegans dir\n";
foreach my $folder (@BUILD_SCRATCH_FOLDERS) {
  my $source = "$autoace/$folder";
  my $target = "$elegans_dir/";
  my $cmd = "mv -f $source $target";
  if ($test) {
    print STDERR "TEST: would have run: $cmd\n";
  } else {    
    $wormbase->run_command($cmd, 'no_log');
  }
}

#############################################################################
# Step 2: move the remaining build dir to the DATABASES directory
#############################################################################
print STDERR "Moving remainder of autoace to DATABASE directory\n";
my $new_dir = $wormbase->wormpub."/DATABASES/".$wormbase->get_wormbase_version_name;
my $cmd = "mv $autoace $new_dir";
if ($test) {
  print STDERR "TEST: would have run: $cmd\n";
  $new_dir = $autoace;
} else {
  $wormbase->run_command($cmd, 'no_log');
}

#############################################################################
# Step 3: remove redundant files that we do not need any more
#############################################################################
print STDERR "Removing files in $new_dir/database/new/\n";
if ($test) {
  print STDERR "TEST: would have deleted $new_dir/database/new/*\n";
} else {
  $wormbase->delete_files_from("$new_dir/database/new","*","+") if (-e "$new_dir/database/new");
}

print STDERR "Removing files in $new_dir/database/touched/\n";
if ($test) {
  print STDERR "TEST: would have deleted $new_dir/databases/touched/*\n";
} else {
  $wormbase->delete_files_from("$new_dir/database/touched","*","+") if (-e "$new_dir/database/touched");
}

print STDERR "Removing files in $new_dir/database/backup*\n";
if ($test) {
  print STDERR "TEST: would have deleted $new_dir/database/backup*\n";
} else {
  $wormbase->delete_files_from("$new_dir/database/","backup","-");
}

print STDERR "Removing files in $new_dir/release/*\n";
if ($test) {
  print STDERR "TEST: would have deleted $new_dir/release/*\n";
} else {
  $wormbase->delete_files_from("$new_dir/release/","*","-");
}


#############################################################################
# Step 4: update the current_db symlink
#############################################################################
my $rm_cmd = "rm " . $wormbase->database('current');
my $ln_cmd = "ln -sf $new_dir ". $wormbase->database('current');
if ($test) {
  print STDERR "TEST: would have run : $rm_cmd\n";
  print STDERR "TEST: would have run : $ln_cmd\n";
} else {
  $wormbase->run_command("rm ".$wormbase->database('current'), 'no_log');
  $wormbase->run_command("ln -sf ".$new_dir." ".$wormbase->database('current'), 'no_log');
}

#############################################################################
# Step 5: remove very old data from DATABASE directory
#############################################################################
#print STDERR "Removing old releases\n";
#for(my $ws = $WS_current - 5; $ws >=1 ; $ws--) {
#  my $old_name = "WS{$ws}";
#  my $old_path = $wormbase->wormpub . "/DATABASES/$old_name";
#  if (-d $old_path) {
#    if ($test) {
#      print STDERR "TEST: would have deleted $old_path\n";
#    } else {
#      $wormbase->run_command("rm -fr $old_path", 'no_log');
#    }
#  }
#}

#############################################################################
# Step 6: zip up files from remaining old directories
#############################################################################
print STDERR "Zipping up .fa and .gff files in older releases\n";
for(my $x = $WS_current - 1; $x >=1 ; $x--) {
  my $old_dir = $wormbase->wormpub."/DATABASES/WS${x}";
  if (-d $old_dir) {
    foreach my $file (glob("$old_dir/*/*.gff"), 
                      glob("$old_dir/*/*.gff2"),
                      glob("$old_dir/*/*.gff3"),
                      glob("$old_dir/*/*.fa")) {
      if ($test) {
        print STDERR "TEST: would have gzipped $file\n";
      } else {
        $wormbase->run_command("gzip $file", 'no_log');
      }
    }
  }
}

##################
# End
##################

$wormbase->mail_maintainer("BUILD: FINISHED", 'All') if not $debug;
print STDERR "Finished\n";

exit(0);



__END__

