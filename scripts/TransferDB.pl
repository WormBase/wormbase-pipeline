#!/usr/local/bin/perl5.6.1 -w

# TransferDB.pl

# by ag3 [991221]
#
# Last updated on: $Date: 2003-05-02 16:18:03 $
# Last updated by: $Author: krb $


# transferdb moves acedb database files across filesystems.
# Creates a temporary database.BCK 
# Checks after every copy on return value and file size
# Updates display.wrm

# touch logfile for run details
$0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

use strict;
use IO::Handle;
use File::Find;
use File::Path;
use Mail::Mailer;
use Getopt::Long;
use Cwd;

my $start="";
my $end="";
my $name="";
my $log="";
my $bck="";
my $mail="";
my $database="";
my $all="";
my $wspec="";
my $wgf="";
my $wscripts="";
my $wquery="";
my $wtools="";
my $whelp="";
my $wdata="";
my $chromosomes="";
my $release="";

my ($srcdir,$enddir,$backup,$printlog,$S_database,
    $S_all,$S_wspec,$S_wgf,$S_wscripts,$S_wquery,$S_wtools,$S_whelp,
    $S_wdata,$S_chromosomes,$S_release,$file);   

our $dbname;

GetOptions (
	    "start=s"   => \$srcdir,
	    "end=s"     => \$enddir,
	    "name=s"    => \$dbname,
	    "bck"       => \$backup,
	    "mail=s"    => \$mail,
	    "printlog"  => \$printlog,
	    "database"  => \$S_database,
	    "all"       => \$S_all,
	    "wspec"     => \$S_wspec,
	    "wgf"       => \$S_wgf,
	    "wscripts"  => \$S_wscripts,
	    "wquery"    => \$S_wquery,
	    "wtools"    => \$S_wtools,
	    "whelp"     => \$S_whelp,
	    "wdata"     => \$S_wdata,
	    "chromosomes"     => \$S_chromosomes,
	    "release"     => \$S_release,
);

# warn if mandatory command line arguments not set
&PrintHelp if(!$srcdir);
&PrintHelp if(!$enddir);
&PrintHelp if ($S_all && !$dbname);



my $DEBUG = 0;
my $OK=0;
our $DB=0;
my $date = `date +%y%m%d`; chomp $date;
my $LOGFILE = "/wormsrv2/logs/TransferDB.${date}.$$";
my %SIG;

$SIG{INT}= sub {print LOGFILE ".. INTERRUPT ..\n"; close LOGFILE; die()};
$SIG{TERM}= sub {print LOGFILE ".. INTERRUPT ..\n"; close LOGFILE; die()}; 
$SIG{QUIT}= sub {print LOGFILE ".. INTERRUPT ..\n"; close LOGFILE; die()};

open (LOGFILE,">$LOGFILE");
LOGFILE->autoflush;

my $datestamp=&GetTime();
my ($s_dir, $e_dir);

print LOGFILE "TransferDB run $$ started at $datestamp\n\n"; 

# set source directory - $s_dir
if ($srcdir =~ /^\~/) {
  $s_dir = glob ($srcdir);
  print LOGFILE "SOURCE DIR : $s_dir\n";
} 
else {
  $s_dir=$srcdir;
  print LOGFILE "SOURCE DIR: $s_dir\n";
}
if (!-d $s_dir) {&SendMail ("ERROR: Could not found directory $s_dir\n");}


# set target directory - $e_dir
if ($enddir =~ /^\~/) {
  $e_dir = glob ($enddir);
  print LOGFILE "TARGET DIR: $e_dir\n\n";
} 
else {
  $e_dir=$enddir;
  if ($e_dir !~ /^\//) {
    my $cwd;
    $e_dir="$cwd"."/$e_dir";
  }
  print LOGFILE "TARGET DIR: $e_dir\n\n";
}

my $new_subdir="$e_dir"."/database";
my $bck_subdir="$e_dir"."/database.BCK";
$database = "$s_dir"."/database";
$wspec = "$s_dir"."/wspec";
$wgf = "$s_dir"."/wgf";
$wscripts = "$s_dir"."/wscripts";
$wquery = "$s_dir"."/wquery";
$wtools = "$s_dir"."/wtools";
$whelp = "$s_dir"."/whelp";
$wdata = "$s_dir"."/data";
$chromosomes = "$s_dir"."/CHROMOSOMES";
$release = "$s_dir"."/release";

my @TOBEMOVED;

$S_all && do {@TOBEMOVED=("$database","$wspec","$wgf","$wscripts","$wtools","$whelp","$wdata","$chromosomes","$release");$OK=1;$DB=1;};
$S_database && do {push (@TOBEMOVED,"$database");$OK=1;};
$S_wspec && do {push (@TOBEMOVED,"$wspec");$OK=1;$DB=1;};
$S_wgf && do {push (@TOBEMOVED,"$wgf");$OK=1;};
$S_wscripts && do {push (@TOBEMOVED,"$wscripts");$OK=1;};
$S_wquery && do {push (@TOBEMOVED,"$wquery");$OK=1;};
$S_wtools && do {push (@TOBEMOVED,"$wtools");$OK=1;};
$S_whelp && do {push (@TOBEMOVED,"$whelp");$OK=1;};
$S_wdata && do {push (@TOBEMOVED,"$wdata");$OK=1;};
$S_chromosomes && do {push (@TOBEMOVED,"$chromosomes");$OK=1;};
$S_release && do {push (@TOBEMOVED,"$release");$OK=1;};
if ($OK == "0") {&PrintHelp;};

print LOGFILE "Directories to be copied: @TOBEMOVED \n";

# Make the backup directory
# if we requested to copy over database directory also

my @OLDDATABASE;

if ($DB==1) {
  if (!-d $e_dir){ 
    mkdir($e_dir,07777) or &SendMail("ERROR 1: Could not mkdir $e_dir: $!\n");
    print LOGFILE "CREATED TARGET DIR: $e_dir\n";
  } 
  else {
    if (-d $new_subdir) {
      @OLDDATABASE = ("$database");
      print LOGFILE "Making backup copy of old database ...\n";
      find (\&backup_db,@OLDDATABASE);
    } 
  }
}


# Move the actual acedb tree structure
find (\&process_file,@TOBEMOVED);

# Remove the backup database directory, unless bck switch specified
if ((!$backup) && (-d $bck_subdir)) {
  print LOGFILE "REMOVING BACKUP TREE $bck_subdir\n";
  rmtree ($bck_subdir) if !$DEBUG;
  print LOGFILE "$0 ended SUCCESSFULLY\n";
  close LOGFILE;
} 
elsif ($backup &&(-d $bck_subdir)) {
  print LOGFILE "Backup database directory is in $bck_subdir\n";
}


my $body = "SUCCESS: Your transferdb process $$ has ended succesfully.\n";
&SendMail($body);

exit 0;

#----------------------------------------
# Makes a backup copy of the old database
# directory. Checks on return value and 
# file size
#
sub backup_db {
  my $filename=$_;
  chomp $filename;
  my $file="$File::Find::name";
  my ($newfile,$bk_chk,$bk_val,$O_SIZE,$N_SIZE); 
  if (-d $file) {
    $file=~s/$database/$bck_subdir/;
    if (!-d $file){
      mkdir($file,07777) or &SendMail ("Could not mkdir backup directory $file: $!");
      print LOGFILE "CREATED backup directory $file\n";
    }
  } 
  else {
    $newfile=$file;
    $newfile=~s/$database/$bck_subdir/;
    if ($file !~ /^\./) {
      $bk_chk="0";
      $bk_val=system("\/usr/apps/bin/scp $file $newfile") if !$DEBUG;
      $bk_chk=$bk_val >> 8;
      $O_SIZE = (stat($file))[7];
      $N_SIZE = (stat($newfile))[7];
      if (($bk_chk != 0)||($O_SIZE != $N_SIZE)) {
	\&SendMail ("Copy of backup $file FAILED") if !$DEBUG;
      } else {
	print LOGFILE "CREATED BACKUP copy - $newfile ..\n";
      }      
    } else {
      print LOGFILE "SKIPPING BACKUP copy - $file ..\n";
    }
  }
}

#---------------------------------
# General copy-over subroutine.
# Checks return value and 
# file size comparing with source.
# Also updates display.wrm
#
sub process_file {
  my $filename=$_;

#  print "File is $file\n";
  my $s_subdir="$File::Find::dir";
  if (!-d $s_subdir) {
    &SendMail ("ERROR: Could not read $s_subdir\n");
  }
  $s_subdir=~s/$s_dir//;

  my $s_file="$File::Find::name";
  my $e_subdir="$e_dir"."$s_subdir";
  if (!-d $e_subdir){
    mkdir($e_subdir,07777) or &SendMail ("ERROR: Could not mkdir subdir $e_subdir: $!\n");
    print LOGFILE "CREATED SUBDIR $e_subdir\n";
  }
  my $e_file="$e_subdir"."/"."$filename";

  if ((-d $s_file) || ($filename =~ /^\./) || ($filename =~ /lock.wrm/)){ 
    print LOGFILE " .. SKIPPING file $filename \n";
  } 
  else {
    if (($filename =~ m/displays.wrm$/) && ($DB==1)){
      print LOGFILE "Updating displays.wrm ...\n";
      open (INFILE,"cat $s_file|");
      open (OUTFILE,">$e_file");
      while (<INFILE>) {
	if ($_ =~ /DDtMain/) {
	  my $string="_DDtMain -g TEXT_FIT -t \"$dbname $datestamp\"  -w .46 -height .32 -help acedb";
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
	$cp_val = system("\/usr/apps/bin/scp -R $s_file $e_file") if !$DEBUG;
      }
      else{
	$cp_val = system("\/usr/apps/bin/scp $s_file $e_file") if !$DEBUG;
      }
      $cp_chk = $cp_val >> 8;
      my $S_SIZE = (stat($s_file))[7];
      my $E_SIZE = (stat($e_file))[7];
      if (($cp_chk != 0)||($S_SIZE != $E_SIZE)){
	\&SendMail ("ERROR: COPY of $s_file FAILED\n") if !$DEBUG;
      }
      else {
	print LOGFILE "COPIED - $e_file .. \n";
      }     
    } 
    else {
      print LOGFILE "SKIPPING COPY of $s_file .. \n";
    };
  }
}

#-------------------------------
# Update log with error message;
# Send  mail with log
#
sub SendMail {
  my $mailbody = shift @_;
  $mailbody.="\n\n";
  my $maintainer = "$mail\@sanger.ac.uk";
  my $from = "TransferDB";
  my $subj = "TransferDB Job $$";
  print LOGFILE $mailbody;
  close LOGFILE;
  if (length ($mail)!=0) {
    open (LOGFILE,"cat $LOGFILE |");
    while (<LOGFILE>) {
      $mailbody .= $_;
    }
    close LOGFILE;
    my $mailer = Mail::Mailer->new();
    $mailer->open({From => $from,
		   To   => $maintainer,
		   Subject => $subj,
		  }) or die "Can't open: $!\n";
    print $mailer $mailbody;
    $mailer->close();
  }
  if ($printlog) {
    open (LOGFILE,"cat $LOGFILE |");
    while (<LOGFILE>) {
      print $_;
    }
    close LOGFILE;
    unlink $LOGFILE;
  } 
exit 0;
}

#------------------------
# Get time coordinates
#
sub GetTime {
  my @time = localtime();
  my ($SECS,$MINS,$HOURS,$DAY,$MONTH,$YEAR)=(localtime)[0,1,2,3,4,5];
  if ($time[1]=~/^\d{1,1}$/) {
    $time[1]="0"."$time[1]";
  }
  my $REALMONTH=$MONTH+1;
  my $REALYEAR=$YEAR+1900;
  my $TODAY = "$DAY $REALMONTH $REALYEAR at $HOURS:$MINS";
  my $NOW="$time[2]"."$time[1]";
  return $TODAY;
}

#---------------------------
# Prints help and disappears
#
sub PrintHelp {
   exec ('perldoc',$0);
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

=item -mail username, if you want the log mailed to username

=item -printlog, if you want to dump the log to STDOUT (useful in scripts)

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
