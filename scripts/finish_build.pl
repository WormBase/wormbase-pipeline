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
# Last updated on: $Date: 2005-12-16 13:41:11 $


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
	    "store"      => \$store,
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
# $debug = `whoami` unless ($debug);


}

# establish log file.
my $log = Log_files->make_build_log($wormbase);




#################################################################################
# variables                                                                     #
#################################################################################

my $currentdb   = glob("~wormpub/DATABASES/current_DB"); 

# the following point to various BUILD directories (or TEST_BUILD if $test)
my $basedir     = $wormbase->basedir;
my $ace_dir     = $wormbase->autoace; # AUTOACE DATABASE DIR
my $chromosomes = $wormbase->chromosomes; # AUTOACE CHROMSOMES
my $logs        = $wormbase->logs; # AUTOACE LOGS

# need to change this if in test mode
my $WS_name;
my $WS_current;

if ($test) {
    $WS_current = "666";
    $WS_name    = "WS666";
} else {
    $WS_current = $wormbase->get_wormbase_version;
    $WS_name    = $wormbase->get_wormbase_version_name;
}

my $WS_new      = $WS_current + 1;
my $WS_new_name = "WS".$WS_new;
my $WS_oldest   = $WS_current - 3; # the version that *should* be the oldest
my $WS_old_name = "WS".$WS_oldest;
my $WS_old_path = "$basedir"."/$WS_old_name";
my $old_wormpep = "$basedir/WORMPEP/wormpep".($WS_current-3);

#####################################################################################


&archive_old_releases;

# update all Common_data files - see Commom_data.pm
$wormbase->run_script("update_Common_data.pl -build -all", $log) 
    && die "Couldn't run update_Common_data.pl -update -in_build -all\n";


# Transfer autoace to WSxx
$log->write_to("Transferring autoace into $basedir/$WS_name\n");
$wormbase->run_script("TransferDB.pl -start $ace_dir -end $basedir/$WS_name -database -release -wspec -chromosomes -acefiles -name $WS_name", $log) 
    && die "couldn't run TransferDB for autoace\n";

# Transfer autoace to ~wormpub/DATABASES/current_DB - first remove existing files
$log->write_to("Removing $currentdb/database/\n");
$wormbase->delete_files_from("$currentdb/database","*","+") unless ($test);

$log->write_to("Removing $currentdb/database/CHROMOSOMS/\n");
$wormbase->delete_files_from("$currentdb/CHROMOSOMES","*","+") unless ($test);


$log->write_to("Removing $currentdb/database/COMMON_DATA/\n");
$wormbase->delete_files_from("$currentdb/COMMON_DATA","*","+") unless ($test);

$log->write_to("Running TransferDB.pl to copy autoace to ~wormpub/DATABASES/current_DB\n");
if (!$test) {
  $wormbase->run_script("TransferDB.pl -start $ace_dir -end $currentdb -database -chromosomes -common -wspec -name $WS_name", $log)  
      && die "couldn't run TransferDB for wormpub\n";
}
$log->write_to("Unzipping any gzipped chromosome files\n");
if (!$test) {
  system("/bin/gunzip $currentdb/CHROMOSOMES/*.gz") 
      && die "Couldn't gunzip CHROMOSOMES/*.gz\n";
}

# Remove redundant files and directories in $ace_dir
$log->write_to("Removing old files in $ace_dir/release/\n");
$wormbase->delete_files_from("$ace_dir/release","*","-");

$log->write_to("Removing old files in $ace_dir/acefiles/\n");
$wormbase->delete_files_from("$ace_dir/acefiles","*","-");

$log->write_to("Removing old CHROMOSOME files in $chromosomes/\n");
$wormbase->delete_files_from("$chromosomes","*","-");

$log->write_to("Removing *.wrm files in $ace_dir/database/\n");
$wormbase->delete_files_from("$ace_dir/database",".\.wrm","-") or $log->write_to("ERROR: Problems removing files from autoace/database\n");

$log->write_to("Removing files in $ace_dir/database/new/\n");
$wormbase->delete_files_from("$ace_dir/database/new","*","+");
$log->write_to("Removing files in $ace_dir/database/touched/\n");
$wormbase->delete_files_from("$ace_dir/database/touched","*","+");

$log->write_to("Removing old files in $logs\n");
$wormbase->delete_files_from("$logs",":","-");
# Exception needed because we need to keep one file (Primary_databases) and this log file uses different name, so
# *:* won't remove it.
unlink("$logs/UTR_gff_dump");

#remove coordinate files for Coords_converter ready for next build.
unlink("$ace_dir/clone_coords") if -e "$ace_dir/clone_coords";
unlink("$ace_dir/SL_coords") if -e "$ace_dir/SL_coords";

#remove ace files created at start of build by make_acefiles.
$log->write_to("\nRemoving ace files created by make_acefiles from . .\n");
my @source_dbs = qw( briggsae caltech camace csh geneace stlace );
foreach my $db (@source_dbs) {
  $log->write_to("\t$db\n");
  $wormbase->delete_files_from("$basedir/wormbase/$db/","ace\$","-");
}

# archive old GFF splits directory'
$log->write_to("Archiving GFFsplits directory using GFFsplitter.pl -a\n\n");
$wormbase->run_script("GFFsplitter.pl -archive", $log) 
    && die "Couldn't run GFFsplitter.pl -archive\n";



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
  $log->write_to("\n\n");

  my @list;
  my $file;

  # remove some unnecessary files before archiving
  if (-d "$WS_old_path") {

      if (-d "$WS_old_path/database") {
	  $log->write_to("Removing $WS_old_path/database\n");
	  
	  # remove contents of the database folder
	  $wormbase->delete_files_from("$WS_old_path/database/new","*","+");
	  $wormbase->delete_files_from("$WS_old_path/database/oldlogs","*","+");
	  $wormbase->delete_files_from("$WS_old_path/database/readlocks","*","+");
	  $wormbase->delete_files_from("$WS_old_path/database/touched","*","+");
	  $wormbase->delete_files_from("$WS_old_path/database","*","+");
      }
      
      # remove wspec, wgf etc.
      $log->write_to("Removing $WS_old_path/wspec\n");
      $wormbase->delete_files_from("$WS_old_path/wspec","*","+");
            
      if (-d "$WS_old_path/pictures"){
	  $log->write_to("Removing $WS_old_path/pictures\n");
	  $wormbase->delete_files_from("$WS_old_path/pictures","*","+");
      }    
  }
  

  # turn the old release into a tarball, move into /nfs/wormarchive and remove old directory
  $log->write_to("\nCreating $WS_old_path.tar.gz\n");
  system ("tar -P /$basedir/ -cvf $WS_old_path.tar $WS_old_path/") 
      && die "Couldn't create tar file\n";
  system ("gzip $WS_old_path.tar") 
      && die "Couldn't create gzip file\n";
  $log->write_to("Moving archive to /nfs/wormarchive and removing $WS_old_path\n");
  system ("mv $WS_old_path.tar.gz /nfs/wormarchive/") 
      && die "Couldn't move to /nfs/wormarchive\n";
  $wormbase->delete_files_from("$WS_old_path","*","+");
 
  # archive old wormpep version
  if (-d $old_wormpep) {
      $log->write_to("\nCreating $old_wormpep archive\n");
      system ("tar -P /$basedir/ -cvf $old_wormpep.tar $old_wormpep") 
	  && die "Couldn't create wormpep tar file\n";
      system ("gzip $old_wormpep.tar") 
	  && die "Couldn't create gzip wormpep tar file\n";
      $wormbase->delete_files_from("$old_wormpep","*","+");
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
