#!/usr/local/bin/perl5.6.0 -w
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


#################################################################################
# variables                                                                     #
#################################################################################

$|=1;
use IO::Handle;
use Getopt::Std;
use Cwd;
use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;

##############################
# Script variables (run)     #
##############################

my $maintainers = "dl1\@sanger.ac.uk krb\@sanger.ac.uk kj2\@sanger.ac.uk";
my $rundate    = `date +%y%m%d`; chomp $rundate;
my $runtime    = `date +%H:%M:%S`; chomp $runtime;
my $version = &get_cvs_version("$0");

##############################
# Paths for I/O files        #
##############################

my $tace   = "/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/tace";
my $giface = "/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/giface";



##############################
# command-line options       #
##############################

our ($opt_h,$opt_b,$opt_c,$opt_s,$opt_i);

#options with arguments
getopt("bcsi");
#boolean options
getopts("h");

&usage if ($opt_h);
&usage if (!defined($opt_b) && !defined($opt_c) && !defined($opt_s) && !defined($opt_i));

##############################################################
# loop through main routine for each database to be unpacked #
##############################################################


&unpack_stuff("brigace",$opt_b) if ($opt_b);
&unpack_stuff("cshace",$opt_c) if ($opt_c);
&unpack_stuff("stlace",$opt_s) if ($opt_s);
&unpack_stuff("citace",$opt_i) if ($opt_i);

sub unpack_stuff{
  my $database = shift;
  my $opt = shift;
  my (@filenames,$filename,$ftp,$dbdir,$logfile,$dbname);


  # set up correct locations
  if ($database eq "brigace"){
    $ftp    = "/nfs/privateftp/ftp-wormbase/pub/incoming/stl";
    $dbdir  = "/wormsrv2/brigace";
    $logfile = "/wormsrv2/logs/unpack_briggsae.$rundate.$$";
    $dbname = "brigdb";
  }

  if ($database eq "cshace"){
    $ftp    = "/nfs/privateftp/ftp-wormbase/pub/incoming/csh";
    $dbdir  = "/wormsrv2/cshace";
    $logfile = "/wormsrv2/logs/unpack_cshace.$rundate.$$";
    $dbname = "cshl_dump";
  }

  if ($database eq "stlace"){
    $ftp    = "/nfs/privateftp/ftp-wormbase/pub/incoming/stl";
    $dbdir  = "/wormsrv2/stlace";
    $logfile = "/wormsrv2/logs/unpack_stlace.$rundate.$$";
    $dbname = "stlace";
  }


  if ($database eq "citace"){
    $ftp    = "/nfs/privateftp/ftp-wormbase/pub/incoming/caltech";
    $dbdir  = "/wormsrv2/citace";
    $logfile = "/wormsrv2/logs/unpack_stlace.$rundate.$$";
    $dbname = "citace_dump";
  }

  ##############################
  # open logfile               #
  ##############################

  my $today = "20" . substr($opt,0,2) . "-" . substr($opt,2,2) . "-" . substr($opt,4,2);

  system ("/bin/touch $logfile") && die "Couldn't touch $logfile\n";
  open (LOGFILE,">>$logfile") || die "Couldn't append to $logfile";
  LOGFILE->autoflush;


  print LOGFILE "# unpack_db.pl - $database\n#\n";
  print LOGFILE "# version             : $version\n";
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

  chdir $dbdir;
  my $dir = cwd();
  print LOGFILE "Move to directory: '$dir'\n";

  # copy database.tar.gz file & check size
  system ("cp -f $ftp/".$dbname."_$today.tar.gz .") && die "Couldn't run cp command\n";



  my $match = &copy_check("$ftp/".$dbname."_$today.tar.gz","$dbdir/".$dbname."_$today.tar.gz");
  print LOGFILE "Copy '".$dbname."_$today.tar.gz' to $dbdir successful\n" if ($match == 1); 
  print LOGFILE "Copy '".$dbname."_$today.tar.gz' to $dbdir failed\n" if ($match == 0);

  system ("/bin/gzip -d ".$dbname."_$today.tar.gz") && die "Couldn't run gzip command\n";
  print LOGFILE "uncompress file\n";

  system ("/bin/tar -xvf ".$dbname."_$today.tar") && die "Couldn't run tar command\n";
  print LOGFILE "untar file\n\n";

  print LOGFILE "Database files to be loaded:\n";


  # add list of ace files to be loaded into array
#  open (LIST, "/bin/ls dump*$today*ace |") || die "Couldn't pipe";
  open (LIST, "/bin/ls *ace |") || die "Couldn't pipe";
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

  system ("mv $dbdir/wspec/displays.wrm $dbdir/wspec/displays.old") && die "Couldn't run mv command\n";

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

  system ("\\rm $dbdir/database/new/*") && die "Couldn't run rm command\n";

  system ("\\rm $dbdir/database/touched/*") && die "Couldn't run rm command\n";

  if (-e "$dbdir/database/lock.wrm") {
    print LOGFILE "*Reinitdb error - lock.wrm file present..\n";
    close LOGFILE;
    die();
  }
  system ("\\mv $dbdir/database/log.wrm $dbdir/database/log.old") && die "Couldn't run rm command\n";

  system ("\\rm $dbdir/database/*.wrm") && die "Couldn't run rm command\n";

  unlink "$dbdir/database/ACEDB.wrm";

  my $command=<<EOF;
y
EOF
    
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

  system ("/bin/rm -f ".$dbname."_$today.tar") && die "Couldn't run rm command\n";

  print LOGFILE "\ndelete tar file from current directory\n\n";

  print LOGFILE "Database files to be removed:\n";
  foreach $filename (@filenames) {
    print LOGFILE "$filename\n";
    unlink $filename;
  }

  $runtime = `date +%H:%M:%S`; chomp $runtime;
  print LOGFILE "\n$database build complete at $rundate $runtime\n";
  close LOGFILE;


  ###############################
  # Mail log to curator         #
  ###############################

  &mail_maintainer("WormBase Report: unpack_$database",$maintainers,$logfile);

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
to the appropriate subdirectory of /wormsrv2.  It then unpacks the 
tar.gz file and re-initialises the ACEDB database. After uploading 
the dump files the directory is cleaned of all uploaded .ace files 
and the .tar file.

unpack_db.pl arguments:

=over 4

=item *

-b <date> unpack St. Louis C. briggsae data and read into brigace

=back

=over 4

=item *

-c <date> unpack CSHL data and read into cshace

=back

=over 4

=item *

-s <date> unpack St. Louis C. elegans data and read into stlace

=back

=over 4

=item *

-i <date> unpack Caltech data and read into citace

=back

=over 4

=item *

-h help page (what you are reading now).

=back

-b, -c, -s, and -i options each require a 6-figure datestamp (e.g. 001106).  You can 
work with individual databases or work with all three in one go (note though that
the script will always process each database in the following order: brigace, cshace, 
stlace, citace).

Example usage:

unpack_db.pl -s 010719 -b 010725

This will unpack the St. Louis C. elegans data and the St. Louis C. briggsae data and 
read them into their respective databases (stlace and brigace).  In this example it will 
try to unpack files which contain the corresponding dates as part of their filename.

=cut
