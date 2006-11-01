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
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2006-11-01 14:06:11 $


#################################################################################
# variables                                                                     #
#################################################################################

use IO::Handle;
use Getopt::Long;
use Cwd;
use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use File::Copy;
use File::Path;
use Storable;

##############################
# command-line options       #
##############################

my $help;
my $brigace;
my $citace;
my $cshace;
my $stlace;
my ($debug,$test,$database,$basedir);
my $store;

GetOptions (
            "help"      => \$help,
            "citace=s"  => \$citace,
            "cshace=s"  => \$cshace,
            "stlace=s"  => \$stlace,
            "brigace=s" => \$brigace,
            "test"      => \$test,
	    "debug:s"   => \$debug,
	    "database|basedir:s"  => \$basedir,
	    "store:s"   => \$store,
	   );

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}

my $rundate = $wormbase->rundate;

&usage if ($help);
&usage if (!defined($brigace) && !defined($citace) && !defined($stlace) && !defined($cshace));

# establish log file.
my $log = Log_files->make_build_log($wormbase);

##############################################################
# loop through main routine for each database to be unpacked #
##############################################################

&unpack_stuff("citace",$citace)   if ($citace);
&unpack_stuff("cshace",$cshace)   if ($cshace);
&unpack_stuff("stlace",$stlace)   if ($stlace);
&unpack_stuff("brigace",$brigace) if ($brigace);


sub unpack_stuff{
  my $database = shift;
  my $opt = shift;
  my (@filenames,$filename,$ftp,$dbdir,$logfile,$dbname);

  $log->write_to("Unpacking $database $opt\n");

  my $ftp_dir = $wormbase->ftp_upload;
  my $primaries = $wormbase->primaries;
  my $logs    = $wormbase->logs;

  # set up correct locations
  if ($database eq "cshace"){
    $ftp     = "$ftp_dir/csh";
    $dbdir   = "$primaries/cshace";
    $logfile = "$logs/unpack_cshace.$rundate.$$";
    $dbname  = "cshl_dump";
  }

  if ($database eq "stlace"){
    $ftp     = "$ftp_dir/stl";
    $dbdir   = "$primaries/stlace";
    $logfile = "$logs/unpack_stlace.$rundate.$$";
    $dbname  = "stlace";
  }


  if ($database eq "citace"){
    $ftp     = "$ftp_dir/caltech";
    $dbdir   = "$primaries/citace";
    $logfile = "$logs/unpack_citace.$rundate.$$";
    $dbname  = "citace_dump";
  }

  if ($database eq "brigace"){
    $ftp     = "$ftp_dir/stl";
    $dbdir   = "$primaries/brigace";
    $logfile = "$logs/unpack_brigace.$rundate.$$";
    $dbname  = "brigace";
  }
  ##############################
  # open logfile               #
  ##############################

  # year 2000 pragmatic programming, remember to alter script in 2100 !
  my $today = "20" . substr($opt,0,2) . "-" . substr($opt,2,2) . "-" . substr($opt,4,2);

  my $msg = "# unpack_db.pl - $database\n#\n".
    "# run details         : $rundate ".$wormbase->runtime."\n".
      "# Source directory    : $ftp\n".
	"# Target directory    : $dbdir\n#\n".
	  "# Source file         : ".$dbname."_$today.tar.gz\n".
	    "\n";

  $log->write_to("$msg");


  # check that command line argument is correct and form new date string
  &error(1,$database,$opt) if ((length $opt) != 6);

 ##########################################
 # copy the tar.gz file from the ftp site #
 ##########################################


  # make temp unpack directory
  my $unpack_dir = "$dbdir"."/temp_unpack_dir";
  eval { mkpath($unpack_dir) };
  if ($@) {
    print "Couldn't create $unpack_dir: $@";
    $log->write_to("Couldn't make temp unpack directory\n");
  } 

  chdir $unpack_dir;
  my $dir = cwd();
  $log->write_to("Move to directory: '$dir'\n");


  # copy database.tar.gz file & check size
  my $file_to_copy = "$ftp/".$dbname."_$today.tar.gz";
  my $status = copy("$file_to_copy", ".");
  $log->write_to("ERROR: Couldn't copy file $file_to_copy : $!\n") if ($status == 0);


  my $match = $wormbase->copy_check("${ftp}/${dbname}_${today}.tar.gz","${dbname}_${today}.tar.gz");
  $msg = "Copy '".$dbname."_$today.tar.gz' to $unpack_dir successful\n" if ($match == 1);
  $msg = "ERROR : Copy '".$dbname."_$today.tar.gz' to $unpack_dir FAILED\n"     if ($match == 0);
  $log->write_to("$msg");

  $wormbase->run_command("/usr/local/bin/gtar zxvf ${dbname}_${today}.tar.gz", $log);
  $log->write_to("unzip and untar file\n\n");


  # add list of ace files to be loaded into array
  $log->write_to("Database files to be loaded:\n");
  open (LIST, "/bin/ls *.ace |") || die "Couldn't pipe";
  while (<LIST>) {
    chomp;
    push (@filenames,"$_");
    $log->write_to("$_\n");
  }
  close LIST;
  $log->write_to("\n\n");




  ###################################
  # modify displays.wrm for rebuild #
  ###################################

  unless( -e "$dbdir/wspec/displays.wrm" )  {
    $log->write_to("display.wrm missing");
    my $db = $wormbase->autoace;
    eval system("cp -R $db/wspec $dbdir");
    if (@!) {
      $log->log_and_die("ERROR $db $dbdir copy ");
    }
  }
  system("mv $dbdir/wspec/displays.wrm $dbdir/wspec/displays.old");
  open (FILE_OLD,"$dbdir/wspec/displays.old")  || $log->log_and_die("failed to open $dbdir/wspec/displays.old\n");
  open (FILE_NEW,">$dbdir/wspec/displays.wrm") || $log->log_and_die("failed to open $dbdir/wspec/displays.wrm\n");

  while (<FILE_OLD>) { 
    if (/^_DDtMain/) {
	print FILE_NEW "_DDtMain -g TEXT_FIT -t \"$database $rundate\"  -w .43 -height .23 -help acedb\n";
      } else {
	print FILE_NEW $_;}
  }
  close FILE_OLD;
  close FILE_NEW;
  unlink "$dbdir/wspec/displays.old" ;

  #############################################################################
  # need to edit the passwrd file to allow re-init if not wormpub ( ie test ) #
  #############################################################################
  if ( $wormbase->test ) {
    open( PASS,"<$dbdir/wspec/passwd.wrm") or $log->log_and_die("cant open password file");
    open( NEWPASS,">$dbdir/wspec/passwd.new") or $log->log_and_die("cant open new password file");
    while ( <PASS> ) {
      if ( /^wormpub/ ) { 
	my $user = `whoami`;
	chomp $user;
	print NEWPASS "$user\n";
      }
      print NEWPASS;
    }
    move("$dbdir/wspec/passwd.new","$dbdir/wspec/passwd.wrm") 
  }

 ####################################
 # Re-initialise the ACEDB database #
 ####################################

  $log->log_and_die("*Reinitdb error - lock.wrm file present..\n") if (-e "$dbdir/database/lock.wrm");

  if( -e "$dbdir/database" ) {
    unlink glob("$dbdir/database/new/*") or $log->write_to("ERROR: Couldn't unlink file $dbdir/database/new/ : $!\n");
    unlink glob("$dbdir/database/touched/*") or $log->write_to( "ERROR: Couldn't unlink file $dbdir/database/touched/ : $!\n");

    if( -e "$dbdir/database/log.wrm") {
      $status = move("$dbdir/database/log.wrm", "$dbdir/database/log.old");
      print "ERROR: Couldn't move file: $!\n" if ($status == 0);
    }
    unlink glob("$dbdir/database/*.wrm") or $log->write_to("ERROR: Couldn't run rm command $dbdir/database/*.wrm: $!\n");
  }
  else {
    eval{
      mkdir("$dbdir/database");
    };
    $log->log_and_die(@!) if @!;
  }
  my $command="y\n";
  $log->write_to("* Reinitdb: reinitializing the database ..\n");
  &DbWrite($command,$wormbase->tace,$dbdir,"ReInitDB");


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
     $log->write_to( "* Reinitdb: reading in new database  $filename\n");
      &DbWrite($command,$wormbase->tace,$dbdir,"ParseFile");
    } else {
      $log->write_to("* Reinitdb: $filename is not existent - skipping ..\n");
      next;
    }
  }
  $log->write_to( "\n$database build complete at ".$wormbase->runtime."\n");
}


$log->mail;

# say goodnight Brian
exit(0);



###################################################
# Subroutine for writing to a given database      #   
###################################################

sub DbWrite {
  my ($command,$exec,$dir,$name,$logfile)=@_;
  open (WRITEDB,"| $exec $dir ") or $log->log_and_die("$name DbWrite failed\n");
  print WRITEDB $command;
  close WRITEDB;
}

###################################################
# Errors and pod documentation

sub error {
  my $error = shift;
  my $database = shift;
  my $opt = shift;
  # Error 1 - date directory ($opt) is of incorrect length 
  if ($error == 1) {
    my $runtime = $wormbase->runtime;
    my $rundate = $wormbase->rundate;

    my $msg = "The database.tar.gz file date for $database is incorrect\n" .
      "'$opt' is not a correct (six figure) date format\n\n" . 
	"Exiting early at $rundate $runtime\n";

    $log->log_and_die("$msg");
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
