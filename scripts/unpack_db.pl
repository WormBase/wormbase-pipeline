#!/usr/local/bin/perl5.8.0 -w
#
# unpack_db.pl
# 
# by Keith Bradnam aged 12 and half
#
# Usage : unpack_db.pl [-options]
#
# A PERL wrapper to automate the extraction and building of:
# the C. briggsae database (brigace)
# the St. Louis database (stlace)
# the Cold Spring Harbor Laboratory database (cshace)
# the Caltech database (citace)
#
# Last updated by: $Author: ck1 $
# Last updated on: $Date: 2004-05-10 14:26:41 $


#################################################################################
# variables                                                                     #
#################################################################################

use IO::Handle;
use Getopt::Long;
use Cwd;
use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use File::Copy;

##############################
# Script variables (run)     #
##############################

my $maintainers = "All";
my $rundate    = &rundate;
my $runtime    = &runtime;

##############################
# Paths for I/O files        #
##############################

my $tace   = &tace;
my $giface = &giface;


##############################
# command-line options       #
##############################

my $help;
my $brigace;
my $citace;
my $cshace;
my $stlace;
my $test;

GetOptions (
            "help"      => \$help,
            "brigace=s" => \$brigace,
            "citace=s"  => \$citace,
            "cshace=s"  => \$cshace,
            "stlace=s"  => \$stlace,
            "test"      => \$test
	    );

&usage if ($help);
&usage if (!defined($brigace) && !defined($citace) && !defined($stlace) && !defined($cshace));

# Set top level path for build, different if in test mode
my $basedir = "/wormsrv2";
$basedir    = glob("~wormpub")."/TEST_BUILD" if ($test); 




##############################################################
# loop through main routine for each database to be unpacked #
##############################################################

&unpack_stuff("brigace",$brigace) if ($brigace);
&unpack_stuff("citace",$citace)   if ($citace);
&unpack_stuff("cshace",$cshace)   if ($cshace);
&unpack_stuff("stlace",$stlace)   if ($stlace);


sub unpack_stuff{
  my $database = shift;
  my $opt = shift;
  my (@filenames,$filename,$ftp,$dbdir,$logfile,$dbname);


  # set up correct locations
  if ($database eq "brigace"){
    $ftp     = "/nfs/privateftp/ftp-wormbase/pub/incoming/stl";
    $dbdir   = "$basedir/brigace";
    $logfile = "$basedir/logs/unpack_briggsae.$rundate.$$";
    $dbname  = "brigdb";
  }

  if ($database eq "cshace"){
    $ftp     = "/nfs/privateftp/ftp-wormbase/pub/incoming/csh";
    $dbdir   = "$basedir/cshace";
    $logfile = "$basedir/logs/unpack_cshace.$rundate.$$";
    $dbname  = "cshl_dump";
  }

  if ($database eq "stlace"){
    $ftp     = "/nfs/privateftp/ftp-wormbase/pub/incoming/stl";
    $dbdir   = "$basedir/stlace";
    $logfile = "$basedir/logs/unpack_stlace.$rundate.$$";
    $dbname  = "stlace";
  }


  if ($database eq "citace"){
    $ftp     = "/nfs/privateftp/ftp-wormbase/pub/incoming/caltech";
    $dbdir   = "$basedir/citace";
    $logfile = "$basedir/logs/unpack_citace.$rundate.$$";
    $dbname  = "citace_dump";
  }

  ##############################
  # open logfile               #
  ##############################

  # year 2000 pragmatic programming, remember to alter script in 2100 !
  my $today = "20" . substr($opt,0,2) . "-" . substr($opt,2,2) . "-" . substr($opt,4,2);

  system ("/bin/touch $logfile") && die "Couldn't touch $logfile\n";
  open (LOGFILE,">>$logfile") || die "Couldn't append to $logfile";
  LOGFILE->autoflush;


  print LOGFILE "# unpack_db.pl - $database\n#\n";
  print LOGFILE "# run details         : $rundate $runtime\n";
  print LOGFILE "# Source directory    : $ftp\n";
  print LOGFILE "# Target directory    : $dbdir\n#\n";
  print LOGFILE "# Source file         : ".$dbname."_$today.tar.gz\n";
  print LOGFILE "\n";


  # check that command line argument is correct and form new date string
  &error(1,$database,$logfile,$opt) if ((length $opt) != 6);

 ##########################################
 # copy the tar.gz file from the ftp site #
 ##########################################


  # make temp unpack directory
  my $unpack_dir = "$dbdir"."/temp_unpack_dir";
  system("/bin/mkdir $unpack_dir") && print LOGFILE "Couldn't make temp unpack directory\n";

  chdir $unpack_dir;
  my $dir = cwd();
  print LOGFILE "Move to directory: '$dir'\n";


  # copy database.tar.gz file & check size
  my $status = copy("$ftp/".$dbname."_$today.tar.gz", ".");
  print "ERROR: Couldn't copy file: $!\n" if ($status == 0);


  my $match = &copy_check("${ftp}/${dbname}_${today}.tar.gz","${dbname}_${today}.tar.gz");
  print LOGFILE "Copy '".$dbname."_$today.tar.gz' to $unpack_dir successful\n" if ($match == 1); 
  print LOGFILE "Copy '".$dbname."_$today.tar.gz' to $unpack_dir failed\n"     if ($match == 0);

  system ("/bin/gzip -d ${dbname}_${today}.tar.gz") && die "Couldn't run gzip command\n";
  print LOGFILE "uncompress file\n";

  system ("/bin/tar -xvf ${dbname}_${today}.tar");
  print LOGFILE "untar file\n\n";


  # add list of ace files to be loaded into array
  print LOGFILE "Database files to be loaded:\n";
  open (LIST, "/bin/ls *.ace |") || die "Couldn't pipe";
  while (<LIST>) {
    chomp;
    push (@filenames,"$_");
    print LOGFILE "$_\n";
  }
  close LIST;
  print LOGFILE "\n\n";




  ###################################
  # modify displays.wrm for rebuild #
  ###################################

  $status = move("$dbdir/wspec/displays.wrm", "$dbdir/wspec/displays.old");
  print "ERROR: Couldn't move file: $!\n" if ($status == 0);


  open (FILE_OLD,"$dbdir/wspec/displays.old") || do { 
    print LOGFILE "failed to open $dbdir/wspec/displays.old\n"; 
    die(1);
  };
  open (FILE_NEW,">$dbdir/wspec/displays.wrm") || do { 
    print LOGFILE "failed to open $dbdir/wspec/displays.wrm\n"; 
    die(1);
  };
  while (<FILE_OLD>) { 
    if (/^_DDtMain/) {
      print FILE_NEW "_DDtMain -g TEXT_FIT -t \"$database $rundate\"  -w .43 -height .23 -help acedb\n";
    } else {
      print FILE_NEW $_;}
  }
  close FILE_OLD;
  close FILE_NEW;
  unlink "$dbdir/wspec/displays.old" ;

 ####################################
 # Re-initialise the ACEDB database #
 ####################################

  if (-e "$dbdir/database/lock.wrm") {
      print LOGFILE "*Reinitdb error - lock.wrm file present..\n";
      close LOGFILE;
      die();
  }

  unlink glob("$dbdir/database/new/*") or print LOGFILE "ERROR: Couldn't unlink file: $!\n";
  unlink glob("$dbdir/database/touched/*") or print LOGFILE "ERROR: Couldn't unlink file: $!\n";

  $status = move("$dbdir/database/log.wrm", "$dbdir/database/log.old");
  print "ERROR: Couldn't move file: $!\n" if ($status == 0);
  unlink glob("$dbdir/database/*.wrm") or print LOGFILE "ERROR: Couldn't run rm command: $!\n";

  my $command="y\n";
  print LOGFILE "* Reinitdb: reinitializing the database ..\n";
  &DbWrite($command,$tace,$dbdir,"ReInitDB",$logfile);


 ##############################
 # Upload the .ace dump files #
 ##############################

  foreach $filename (@filenames) { 
    my $command=<<END;
pparse $filename
save 
quit
END
    if (-e $filename) {
      print LOGFILE "* Reinitdb: reading in new database  $filename\n";
      &DbWrite($command,$tace,$dbdir,"ParseFile",$logfile);
    } else {
      print LOGFILE "* Reinitdb: $filename is not existent - skipping ..\n";
      next;
    }
  }
 
 
  
  ###############################
  # Tidy up old ace files       #
  ###############################
  unlink glob("$unpack_dir/*") or print LOGFILE "ERROR: Couldn't remove $unpack_dir/*\n";
  rmdir("$unpack_dir") or print LOGFILE "ERROR: Could't remove $unpack_dir\n";

  $runtime = &runtime;
  print LOGFILE "\n$database build complete at $rundate $runtime\n";
  close LOGFILE;


  ###############################
  # Mail log to curator         #
  ###############################
  my $subject_line = "BUILD REPORT: unpack_$database";
  $subject_line = "TEST BUILD REPORT: WormBase Report: unpack_$database" if ($test);  
  &mail_maintainer("$subject_line",$maintainers,$logfile);

}

# say goodnight Brian
exit(0);



###################################################
# Subroutine for writing to a given database      #   
###################################################

sub DbWrite {
  my ($command,$exec,$dir,$name,$logfile)=@_;
  open (WRITEDB,"| $exec $dir >> $logfile") or do {print LOGFILE "$name DbWrite failed\n";close LOGFILE; die();};
  print WRITEDB $command;
  close WRITEDB;
}

###################################################
# Errors and pod documentation

sub error {
  my $error = shift;
  my $database = shift;
  my $logfile = shift;
  my $opt = shift;
  # Error 1 - date directory ($opt) is of incorrect length 
  if ($error == 1) {
    print "The database.tar.gz file date for $database is incorrect\n";
    print "'$opt' is not a correct (six figure) date format\n\n";
    $runtime = `date +%H:%M:%S`; chomp $runtime;
    print LOGFILE "The database.tar.gz file date is incorrect.\n";
    print LOGFILE "'$opt' is not a correct (six figure) date format\n\n";
    print LOGFILE "Exiting early at $rundate $runtime\n";
    close LOGFILE;
    &mail_maintainer("WormBase Report: unpack_db.pl",$maintainers,$logfile);
  }
  exit(1);
}

sub usage {
    system ('perldoc',$0);
    exit;	
}

__END__

=pod

=head1 NAME - unpack_db.pl

=head2 USAGE

unpack_db.pl is an all-in-one replacement for what was previously
performed by three separate scripts (unpack_briggsae, unpack_cshace,
and unpack stlace). It also unpacks the new citace dump.

This script can be used to replace any or all of the aforementioned
scripts.  For each of the databases that the user specifies it will
move the database dump files from the appropriate incoming FTP site 
to the appropriate subdirectory of /wormsrv2 (or ~wormpub/TEST_BUILD/
if in test mode).  It then unpacks the tar.gz file and re-initialises 
the ACEDB database. After uploading the dump files the directory is 
cleaned of all uploaded .ace files and the .tar file.

unpack_db.pl arguments:

=over 4

=item *

-brigace <date> unpack St. Louis C. briggsae data and read into brigace

=back

=over 4

=item *

-cshace <date> unpack CSHL data and read into cshace

=back

=over 4

=item *

-stlace <date> unpack St. Louis C. elegans data and read into stlace

=back

=over 4

=item *

-citace <date> unpack Caltech data and read into citace

=back

=over 4

=item *

-help, help page (what you are reading now).

=back

=over 4

=item *

-test, run in test mode, outputs to ~wormpub/TEST_BUILD/

=back

-brigace, -citace, -stlace, and -citace options each require a 6-figure datestamp 
(e.g. 001106).  You can work with individual databases or work with all three in one 
go (note though that the script will always process each database in the following order: 
brigace, cshace, stlace, citace).

Example usage:

unpack_db.pl -stlace 010719 -brigace 010725

This will unpack the St. Louis C. elegans data and the St. Louis C. briggsae data and 
read them into their respective databases (stlace and brigace).  In this example it will 
try to unpack files which contain the corresponding dates as part of their filename.

=cut
