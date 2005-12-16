#!/usr/local/bin/perl5.8.0 -w
#
# TransferDB.pl
#
# by ag3 [991221]
#
# Last updated on: $Date: 2005-12-16 11:18:54 $
# Last updated by: $Author: ar2 $

# transferdb moves acedb database files across filesystems.
# Creates a temporary database.BCK 
# Checks after every copy on return value and file size
# Updates display.wrm

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;

use Carp;
use IO::Handle;
use File::Find;
use File::Path;
use Getopt::Long;
use Cwd;
use Log;
use Storable;


##############################
# command-line options       #
##############################

my $help;              # Help/Usage page
my $verbose;           # turn on extra command-line output
my $debug;             # For sending log file output to just one person
my $mail;              # -mail: turn flag on to email log file
my $dbname;            # -name: name of the database (will be overwritten in displays.wrm)
my $srcdir;            # -start: Location of source database, mandatory option
my $enddir;            # -end: Location of target database, mandatory option
my $backup;            # -bck: will specify that script keep .BCK directory
my $S_acefiles;        # -acefiles: will copy acefiles dir
my $S_database;        # -database: will copy database dir
my $S_wspec;           # -wspec
my $S_wgf;             # -wgf
my $S_wscripts;        # -wscripts
my $S_wquery;          # -wquery
my $S_wtools;          # -wtools
my $S_whelp;           # -whelp
my $S_wdata;           # -wdata ??? what is this
my $S_chromosomes;     # -chromosomes: copies CHROMOSOMES dir
my $S_release;         # -release 
my $S_common;          # -common : copies COMMON_DATA
my $S_all;             # -all: all of the above
my $file;              # ???
my $retry = 5;         # for making repeat attempts to copy a file
my $test;              # test mode, logs go to ~wormpub/TEST_BUILD/logs
my $split;             # -split: changes cp function for models.wrm if splitting camace
my $store;
my $wormbase;

GetOptions (
	    "start=s"     => \$srcdir,
	    "end=s"       => \$enddir,
	    "name=s"      => \$dbname,
	    "bck"         => \$backup,
	    "database"    => \$S_database,
	    "all"         => \$S_all,
	    "acefiles"    => \$S_acefiles,
	    "wspec"       => \$S_wspec,
	    "wgf"         => \$S_wgf,
	    "wscripts"    => \$S_wscripts,
	    "wquery"      => \$S_wquery,
	    "wtools"      => \$S_wtools,
	    "whelp"       => \$S_whelp,
	    "wdata"       => \$S_wdata,
	    "chromosomes" => \$S_chromosomes,
	    "release"     => \$S_release,
	    "common"      => \$S_common,
	    "help"        => \$help,
	    "verbose"     => \$verbose,
	    "debug=s"     => \$debug,
	    "mail"        => \$mail,
	    "retry=i"     => \$retry,
	    "test"        => \$test,
            "split"       => \$split,
	    "store:s"     => \$store
);

if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}

##################################
#set up log files and debug mode
##################################

my $log = Log_files->make_build_log( $wormbase );


#########################################################
# Help pod if needed
# or warn if mandatory command line arguments not set
##########################################################
&usage("Help") if ($help);
&usage("1")    if (!$srcdir || !$enddir);
&usage("2")    if (($dbname && !$S_wspec) && ($dbname && !$S_all));

if ($srcdir =~ /^\~/) {
  my $tmp = glob ($srcdir);  $srcdir = $tmp;
} 
&usage("1")    if (!-d $srcdir);



###########################
# catch signal interrupts
###########################
my %SIG;

$SIG{INT}  = sub {$log->log_and_die(".. INTERRUPT ..\n");};
$SIG{TERM} = sub {$log->log_and_die(".. INTERRUPT ..\n");};
$SIG{QUIT} = sub {$log->log_and_die(".. INTERRUPT ..\n");};


#######################################################################
# set target directory - $enddir
# if leading ~ specified, need to use glob
# create directory if it doesn't exist
#######################################################################

if ($enddir =~ /^\~/) {
  my $tmp = glob ($enddir);
  $enddir = $tmp;
} 
else {
  if ($enddir !~ /^\//) {
    my $cwd;
    $enddir="$cwd"."/$enddir";
  }
}
# does it exist?
if (!-d $enddir){ 
  print "$enddir doesn't exist, will try to create it\n";
  $log->write_to( "$enddir doesn't exist, will try to create it\n");
  mkdir($enddir,07777) or &usage("3");
} 


$log->write_to( "SOURCE DIR: $srcdir\n");
$log->write_to( "TARGET DIR: $enddir\n\n");

my $new_subdir  = "$enddir"."/database";
my $bck_subdir  = "$enddir"."/database.BCK";
my $database    = "$srcdir"."/database";
my $wspec       = "$srcdir"."/wspec";
my $wgf         = "$srcdir"."/wgf";
my $wscripts    = "$srcdir"."/wscripts";
my $wquery      = "$srcdir"."/wquery";
my $wtools      = "$srcdir"."/wtools";
my $whelp       = "$srcdir"."/whelp";
my $wdata       = "$srcdir"."/data";
my $chromosomes = "$srcdir"."/CHROMOSOMES";
my $release     = "$srcdir"."/release";
my $acefiles    = "$srcdir"."/acefiles";
my $common_data = "$srcdir"."/COMMON_DATA";

# set what is to be copiedin @TOBEMOVED
my @TOBEMOVED;

push (@TOBEMOVED,"$database")    if ($S_database    || $S_all);
push (@TOBEMOVED,"$wspec")       if ($S_wspec       || $S_all);
push (@TOBEMOVED,"$wgf")         if ($S_wgf         || $S_all);
push (@TOBEMOVED,"$wscripts")    if ($S_wscripts    || $S_all);
push (@TOBEMOVED,"$wquery")      if ($S_wquery      || $S_all);
push (@TOBEMOVED,"$wtools")      if ($S_wtools      || $S_all);
push (@TOBEMOVED,"$whelp")       if ($S_whelp       || $S_all);
push (@TOBEMOVED,"$wdata")       if ($S_wdata       || $S_all);
push (@TOBEMOVED,"$chromosomes") if ($S_chromosomes || $S_all);
push (@TOBEMOVED,"$release")     if ($S_release     || $S_all);
push (@TOBEMOVED,"$acefiles")    if ($S_acefiles    || $S_all);
push (@TOBEMOVED,"$common_data") if ($S_common      || $S_all);

$log->write_to( "Directories to be copied: @TOBEMOVED \n");



#############################################################################
# Make a backup database subdirectory
# Only do this if -database/-all specified and target database subdir exists
#############################################################################

if (($S_database || $S_all) && -d $new_subdir) {
  $log->write_to( "Making backup copy of old database ...\n");
  find (\&backup_db,$database); 
}


#################################################################
# Move the actual acedb tree structure, unless it doesn't exist!
#################################################################
my @failed_files;
foreach my $dir (@TOBEMOVED){
  if (!-d $dir){
    $log->write_to( "$dir doesn't exist, skipping to next category\n");
  }
  else{
    find (\&process_file,$dir);
  }
}


######################################################################
# Remove the backup database directory, unless bck switch specified
######################################################################

if ((!$backup) && (-d $bck_subdir)) {
  $log->write_to( "REMOVING BACKUP TREE $bck_subdir\n");
  rmtree ($bck_subdir);
} 
elsif ($backup &&(-d $bck_subdir)) {
  $log->write_to( "Backup database directory is in $bck_subdir\n");
}

########################
# Finish script cleanly
########################

$log->write_to( "\n=============================================\n");
if( @failed_files ) {
  $log->write_to( "TransferDB process $$ ended in FAILURE  at ",$wormbase->runtime,"\n\n");
  $log->write_to( "\n\nThe following ",scalar @failed_files," files were NOT copied correctly \n");
  foreach ( @failed_files ) {
    $log->write_to( "\t$_\n");
  }
}
else {
  $log->write_to( "TransferDB process $$ ended SUCCESSFULLY at ",$wormbase->runtime,"\n");
}
$log->mail;
exit 0;


###############################################################################
#
#
#     T H E    S U B R O U T I N E S
#
#
###############################################################################


##############################################################################
# backup_db
#
# Makes a backup copy of the old database directory. Checks on return value and 
# file size
###############################################################################

sub backup_db {
  my $filename=$_;
  chomp $filename;
  my $file="$File::Find::name";
  my ($newfile,$bk_chk,$bk_val,$O_SIZE,$N_SIZE); 
  if (-d $file) {
    $file=~s/$database/$bck_subdir/;
    if (!-d $file){
      mkdir($file,07777) or $log->write_to( "Could not mkdir backup directory $file: $!");
      $log->write_to( "CREATED backup directory $file\n");
    }
  } 
  else {
    $newfile=$file;
    $newfile=~s/$database/$bck_subdir/;
    if ($file !~ /^\./) {
      $bk_chk="0";
      $bk_val=system("\/usr/apps/bin/scp $file $newfile");
      $bk_chk=$bk_val >> 8;
      $O_SIZE = (stat($file))[7];
      $N_SIZE = (stat($newfile))[7];
      if (($bk_chk != 0)||($O_SIZE != $N_SIZE)) {
	$log->write_to( "Copy of backup $file FAILED\n");
      } else {
	$log->write_to( "CREATED BACKUP copy - $newfile ..\n");
      }      
    } else {
      $log->write_to( "SKIPPING BACKUP copy - $file ..\n");
    }
  }
}

############################################################################################
# 
# process_file
#
# General copy-over subroutine. Checks return value and 
# file size comparing with source. Also updates display.wrm
#
############################################################################################

sub process_file {
  my $filename=$_;

  my $s_subdir="$File::Find::dir";

  if (!-d $s_subdir) {
    $log->write_to( "ERROR: Could not read $s_subdir\n");
  }
  $s_subdir =~ s/$srcdir//;

  my $s_file="$File::Find::name";
  my $e_subdir="$enddir"."$s_subdir";
  $e_subdir =~ s/\/\//\//;

  if (!-d $e_subdir){
    unless(mkdir($e_subdir,07777)){
      $log->write_to( "ERROR: Could not mkdir subdir $e_subdir: $!\n");
      close(LOG);
      croak "ERROR: Could not mkdir subdir $e_subdir: $!\n";
    }
    $log->write_to( "CREATED SUBDIR $e_subdir\n");
  }
  my $e_file="$e_subdir"."/"."$filename";


  if ((-d $s_file) || ($filename =~ /^\./) || ($filename =~ /lock.wrm/)){ 
    $log->write_to( " .. SKIPPING file $filename \n");
  } 
  else {
    # if you are copying displays.wrm and -name was specified, you can update
    # the contents of the file itself to ad6~d the new name
    if (($filename =~ m/displays.wrm$/) && $dbname){
      $log->write_to( "Updating displays.wrm ...\n");
      open (INFILE,"cat $s_file|");
      open (OUTFILE,">$e_file");
      while (<INFILE>) {
	if ($_ =~ /DDtMain/) {
	  my $rundate = `date +%y%m%d`; chomp $rundate;
	  my $string="_DDtMain -g TEXT_FIT -t \"$dbname $rundate\"  -w .46 -height .32 -help acedb";
	  print OUTFILE $string;
	} 
	else {
	  print OUTFILE $_;
	}      
      }
      close INFILE;
      close OUTFILE;
    }
    elsif ($filename !~ /^\./) {
      #loop to retry copying on failure
      my $success = 0;
    ATTEMPT:
      for(my $i = 1; $i <= $retry ; $i++) {
	my $cp_chk = "0";
	my $cp_val;

	if( ($filename =~ m/models\.wrm$/) && (!$split) ){
	  $cp_val = system("\/usr/bin/cp -R $s_file $e_file") 
	}
	else{
	  $cp_val = system("\/usr/apps/bin/scp $s_file $e_file");
	}
	$cp_chk = $cp_val >> 8;
	
	my $S_SIZE = (stat($s_file))[7];
	my $E_SIZE = (stat($e_file))[7];
	if (($cp_chk != 0)||($S_SIZE != $E_SIZE)){
	  $log->write_to( "ERROR: COPY of $s_file FAILED - attempt $i\n");
	}
	else {
	  $log->write_to( "COPIED - $e_file .. \n");
	  $success = 1;
	  last ATTEMPT;
	}     
      } 
      # end retry loop
      push(@failed_files, $s_file) unless $success == 1;
    }
    else {
      $log->write_to( "SKIPPING COPY of $s_file .. \n");
    };
  }
}


#######################################################################
# Help and error trap outputs                                         #
#######################################################################

sub usage {
  my $error = shift;
  if ($error eq "Help") {
    # Normal help menu
    exec ('perldoc',$0);
  }
  elsif($error == 1){
    print "You must specify:\n";
    print "-start <database_path> and\n";
    print "-end <database_path>\n";
    print "Both database paths must be valid ones\n\n";
    $log->write_to( "-start and/or -end not specified, or -start describes an invalid database path\n");
    $log->write_to( "TransferDB prematurely quitting at".$wormbase->runtime."\n");
    exit(1);
  }
  elsif($error == 2){
    print "If -name is specified you must also specify -wspec or -all\n";
    $log->write_to( "-name specified but -wspec or -all not specified\n");
    $log->write_to( "TransferDB prematurely quitting at".$wormbase->runtime."\n");
    exit(1);
  }
  elsif($error == 3){
    print "ERROR: Could not run mkdir for directory specified by -end\n";
    $log->write_to( "ERROR: Could not run mkdir for directory specified by -end\n");
    $log->write_to( "TransferDB prematurely quitting at".$wormbase->runtime."\n");
    exit(1);
		 }



}

#############################################################################################





__END__

=pod

=head2   NAME - TransferDB

=head1 USAGE

TransferDB copies over AceDB database directory structure and contents.
It creates a temporary backup (database.BCK) of the pre-existing target 
database subdirectory. This directory will be removed when the copy process 
is completed, but it can be kept using the optional -bck switch.
The database backup directory will not be created unless the -all or
-database switches are included.

TransferDB mandatory arguments:

=over 4

=item -start source_database (argument required)

=item -end target_database_dir (argument required)

=item -all, if you want to copy the whole database structure

=item -name new display name, mandatory if you ask to copy all the database

=back

log_filename is assumed to be a filename under current directory

TransferDB OPTIONAL arguments:

=over 4

=item -bck, if you want to keep a backup copy of the database directory

=item -debug <username>, if you want the log mailed to username

=item -verbose, not really used yet, but will give more output to command line

=item -mail, if set will email log file to everyone (default) or just to user specified by -debug
=back

=item -test, if building in -test mode this script will write logs to ~wormpub/TEST_BUILD/logs

Choose one or more of the following switches if you 
want to copy only a subset of the acedb directories, 
in alternative to the -all switch:

=over 4

=item -database  =>  database subdir

=item -wspec     =>  wspec subdir

=item -wgf       =>  wgf subdir

=item -wscripts  =>  wscripts subdir

=item -wquery    =>  wquery subdir

=item -wtools    =>  wtools subdir

=item -whelp     =>  whelp subdir

=item -wdata     =>  wdata subdir 

=item -chromosomes =>  CHROMOSOMES subdir

=item -release =>  release subdir

=item -acefiles => acefiles subdir

=back

=cut
