#!/usr/local/bin/perl5.8.0 -w
#
# TransferDB.pl
#
# by ag3 [991221]
#
# Last updated on: $Date: 2003-09-22 22:43:34 $
# Last updated by: $Author: krb $


# transferdb moves acedb database files across filesystems.
# Creates a temporary database.BCK 
# Checks after every copy on return value and file size
# Updates display.wrm

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;

use Carp;
use IO::Handle;
use File::Find;
use File::Path;
use Getopt::Long;
use Cwd;

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
my $S_all;             # -all: all of the above
my $file;              # ???

GetOptions (
	    "start=s"     => \$srcdir,
	    "end=s"       => \$enddir,
	    "name=s"      => \$dbname,
	    "bck"         => \$backup,
	    "database"    => \$S_database,
	    "all"         => \$S_all,
	    "wspec"       => \$S_wspec,
	    "wgf"         => \$S_wgf,
	    "wscripts"    => \$S_wscripts,
	    "wquery"      => \$S_wquery,
	    "wtools"      => \$S_wtools,
	    "whelp"       => \$S_whelp,
	    "wdata"       => \$S_wdata,
	    "chromosomes" => \$S_chromosomes,
	    "release"     => \$S_release,
	    "help"        => \$help,
	    "verbose"     => \$verbose,
	    "debug=s"     => \$debug,
	    "mail"        => \$mail
	    );

##################################
#set up log files and debug mode
##################################
my $log; 
my $maintainers = "All";

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}
&create_log_files;



#########################################################
# Help pod if needed
# or warn if mandatory command line arguments not set
##########################################################
&usage("Help") if ($help);
&usage("1")    if (!$srcdir || !$enddir);
&usage("2")    if ($S_all && !$dbname);

if ($srcdir =~ /^\~/) {
  my $tmp = glob ($srcdir);  $srcdir = $tmp;
} 
&usage("1")    if (!-d $srcdir);



###########################
# catch signal interrupts
###########################
my %SIG;

$SIG{INT}  = sub {print LOG ".. INTERRUPT ..\n"; close LOG; die()};
$SIG{TERM} = sub {print LOG ".. INTERRUPT ..\n"; close LOG; die()}; 
$SIG{QUIT} = sub {print LOG ".. INTERRUPT ..\n"; close LOG; die()};


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
  print LOG "$enddir doesn't exist, will try to create it\n";
  mkdir($enddir,07777) or &usage("3");
} 


print LOG "SOURCE DIR: $srcdir\n";
print LOG "TARGET DIR: $enddir\n\n";

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

print LOG "Directories to be copied: @TOBEMOVED \n";



#############################################################################
# Make a backup database subdirectory
# Only do this if -database/-all specified and target database subdir exists
#############################################################################

if (($S_database || $S_all) && -d $new_subdir) {
  print LOG "Making backup copy of old database ...\n";
  find (\&backup_db,$database); 
}


#################################################################
# Move the actual acedb tree structure, unless it doesn't exist!
#################################################################
foreach my $dir (@TOBEMOVED){
  if (!-d $dir){
    print LOG "$dir doesn't exist, skipping to next category\n";
  }
  else{
    find (\&process_file,$dir);
  }
}


######################################################################
# Remove the backup database directory, unless bck switch specified
######################################################################

if ((!$backup) && (-d $bck_subdir)) {
  print LOG "REMOVING BACKUP TREE $bck_subdir\n";
  rmtree ($bck_subdir);
} 
elsif ($backup &&(-d $bck_subdir)) {
  print LOG "Backup database directory is in $bck_subdir\n";
}

########################
# Finish script cleanly
########################
print LOG "\n=============================================\n";
print LOG "TransferDB process $$ ended SUCCESSFULLY at ",&runtime,"\n";
close(LOG);
&mail_maintainer("TransferDB.pl",$maintainers,$log) if ($mail);
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
      mkdir($file,07777) or print LOG "Could not mkdir backup directory $file: $!";
      print LOG "CREATED backup directory $file\n";
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
	print LOG "Copy of backup $file FAILED\n";
      } else {
	print LOG "CREATED BACKUP copy - $newfile ..\n";
      }      
    } else {
      print LOG "SKIPPING BACKUP copy - $file ..\n";
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
    print LOG "ERROR: Could not read $s_subdir\n";
  }
  $s_subdir =~ s/$srcdir//;

  my $s_file="$File::Find::name";
  my $e_subdir="$enddir"."$s_subdir";
  $e_subdir =~ s/\/\//\//;

  if (!-d $e_subdir){
    unless(mkdir($e_subdir,07777)){
      print LOG "ERROR: Could not mkdir subdir $e_subdir: $!\n";
      close(LOG);
      croak "ERROR: Could not mkdir subdir $e_subdir: $!\n";
    }
    print LOG "CREATED SUBDIR $e_subdir\n";
  }
  my $e_file="$e_subdir"."/"."$filename";


  if ((-d $s_file) || ($filename =~ /^\./) || ($filename =~ /lock.wrm/)){ 
    print LOG " .. SKIPPING file $filename \n";
  } 
  else {
    # if you are copying displays.wrm and -name was specified, you can update
    # the contents of the file itself to ad6~d the new name
    if (($filename =~ m/displays.wrm$/) && $dbname){
      print LOG "Updating displays.wrm ...\n";
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
      my $cp_chk = "0";
      my $cp_val;
      if($filename =~ m/models\.wrm$/){
	$cp_val = system("\/usr/bin/cp -R $s_file $e_file");
      }
      else{
	$cp_val = system("\/usr/apps/bin/scp $s_file $e_file");
      }
      $cp_chk = $cp_val >> 8;

      my $S_SIZE = (stat($s_file))[7];
      my $E_SIZE = (stat($e_file))[7];
      if (($cp_chk != 0)||($S_SIZE != $E_SIZE)){
	print LOG "ERROR: COPY of $s_file FAILED\n";
	croak "ERROR: COPY of $s_file FAILED\n";
      }
      else {
	print LOG "COPIED - $e_file .. \n";
      }     
    } 
    else {
      print LOG "SKIPPING COPY of $s_file .. \n";
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
    print LOG "-start and/or -end not specified, or -start describes an invalid database path\n";
    print LOG "TransferDB prematurely quitting at",&runtime,"\n";
    close(LOG);
    exit(1);
  }
  elsif($error == 2){
    print "If -all is specified you must also specify -name\n";
    print LOG "-all specified but -dbname not specified\n";
    print LOG "TransferDB prematurely quitting at",&runtime,"\n";
    close(LOG);
    exit(1);
  }
  elsif($error == 3){
    print "ERROR: Could not run mkdir for directory specified by -end\n";
    print LOG "ERROR: Could not run mkdir for directory specified by -end\n";
    print LOG "TransferDB prematurely quitting at",&runtime,"\n";
    close(LOG);
    exit(1);
  }



}

#############################################################################################

sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate     = `date +%y%m%d`; chomp $rundate;
  $log        = "/wormsrv2/logs/$script_name.$rundate.$$";

  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name process $$ started at ",&runtime,"\n";
  print LOG "=============================================\n";
  print LOG "\n";

}




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


=back

=cut
