#!/usr/local/bin/perl5.8.0 -w
#
# make_acefiles.pl 
#
# dl1
#
# Generates the .acefiles from the primary databases as a prelim for building
# autoace.
#
# Last updated by: $Author: krb $                     
# Last updated on: $Date: 2004-03-02 13:34:24 $       

#################################################################################
# Variables                                                                     #
#################################################################################

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use IPC::Open2;
use IO::Handle;
use Getopt::Long;
use File::Copy;

##############################
# command-line options       #
##############################

our $help;       # Help perdoc
our $db;         # Database name for single db option
our $debug;      # Debug mode, verbose output to runner only
my $test;        # If set, script will use TEST_BUILD directory under ~wormpub

GetOptions ("debug=s"   => \$debug,
	    "help"      => \$help,
	    "db=s"      => \$db,
	    "test"      => \$test);


##############################
# Script variables (run)     #
##############################

my $maintainers = "All";
my $rundate     = &rundate;
my $runtime     = &runtime;
my $tace           = &tace;


# set paths to take account of whether -test is being used
my $basedir     = "/wormsrv2";
$basedir        = glob("~wormpub")."/TEST_BUILD" if ($test); 

my $wormbasedir = "$basedir/wormbase";
my $miscdir     = "$wormbasedir/misc";

my $autoace_config = "$basedir/autoace_config/autoace.config";


# set WS version depending on whether in test mode
my $WS_version;
if($test){
  $WS_version  = "WS666";
}
else{
 $WS_version = &get_wormbase_version_name;
}


our $log;


# help page
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

# single database mode
if ($db) {
    &usage("Bad_database_name") unless (($db eq "camace") || ($db eq "stlace") || ($db eq "briggsae") || ($db eq "cshace") || ($db eq "citace") || ($db eq "geneace"));
}


&create_log_files;

# Main Loop   
&mknewacefiles;


##############################
# tidy up                    #
##############################


print LOG &runtime, ": finished running make_acefiles.pl\n\n";
print LOG "================================================================\n\n";


close(STDOUT);
close(STDERR);
close(LOG);

&mail_maintainer("WormBase Report: make_acefiles.pl",$maintainers,$log);

exit (0);





#################################################################################
### Subroutines                                                               ###
#################################################################################


sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch $basedir/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate  = &rundate;
  $log         = "$basedir/logs/$script_name.$WS_version.$rundate.$$";

  open (LOG, ">$log") or die "cant open $log";

  # also want to capture STDOUT AND STERR to logfile
  LOG->autoflush();
  open(STDOUT,">>$log");
  STDOUT->autoflush();
  open(STDERR,">>$log");
  STDERR->autoflush();

  print LOG "# make_acefiles.pl started at: $rundate ",&runtime,"\n";
  print LOG "# TEST MODE!!! Using $basedir/autoace\n" if ($test);
  print LOG "# WormBase/Wormpep version: WS${WS_version}\n\n";
  print LOG "======================================================================\n\n";
 

  if ($debug) {
    print "# make_acefiles.pl\n\n";
    print "# run details    : $rundate $runtime\n";
    print "\n";
    print "WormBase version : ${WS_version}\n";
    print "\n";
    print "======================================================================\n";
    print " Write .ace files\n";
    print "  from database $db\n"                           if ($db);
    print "======================================================================\n";
    print "\n";
    print "Starting make_acefiles.pl .. \n\n";
    print "writing .acefiled to '$wormbasedir'\n";
  }


}



#################################################################################
# Erases old acefiles and make new ones                                         #
#################################################################################

sub mknewacefiles {
    
  my ($dbname,$dbdir,$targetdir,$exe,$object,$criteria,$criterianoasterisk,$follow);
  my ($filename,$extrafile,$filepath,$command,$outfile,$tag,@deletes,$include);
  
  local (*CONFIG);
  
  open (CONFIG, "<$autoace_config");
  while (<CONFIG>) {
    
    # some formating, remove leading and trailing whitespace
    s/^\s+//;
    s/\s+$//;    
    
    # next if comment line
    next if (/^\#/  || /^$/);
    
    # single database mode, next unless line contains db name
    if ($db) {
      next unless (/$db/);
    }
    # parse database information
    if (/^P\s+(\S+)\s+(\S+)$/) {
      $dbname    = $1;
      $dbdir     = $2;
      print LOG "\n\n",&runtime, ": Processing $dbname information in autoace_config\n";
      # need to change dbpath if in test mode
      $dbdir     =~ s/\/wormsrv2/$basedir/ if ($test);
      $targetdir = "$wormbasedir/$dbname";
      $exe       = "$tace $dbdir";
      next;
    }
    
    next if ($dbname =~ /misc/);
    
    # nuts and bolts of acefile generation
    # format:  database object criteria
    
    # clear out the @deletes array eachtime
    undef (@deletes);
    
    # parse filename
    ($filename) = (/^\S+\s+(\S+)/);
    $filepath   = "$wormbasedir/$dbname/$filename";
    $extrafile  = "$wormbasedir/$dbname/$filename.extra";

    print "Noting filename:$filename\n" if ($debug);
    
    # deal with queries requiring tag deletion
    if (/\[(\S+.+)\]$/) { 
      @deletes = split (/\s/,$1);
      print LOG "Adding " . scalar (@deletes) . " tags to delete array\n\n";
      my $report = join ' + ', @deletes;
      
      if (/^\S+\s+\S+\s+(\S+)\s+\[/) {  
	$object = $1; ($criteria,$criterianoasterisk) = ""; 
	
	$command = "nosave\nquery find $object\n";
	foreach my $delete (@deletes) {
	  $command .= "eedit -D $delete\n";
	}
	$command .= "show -a -T -f $filepath\nquit\n";
      } 
      elsif (/^\S+\s+\S+\s+(\S+)\s+(\S+.+)\s+\[/) {
	$object   = $1; ($criteria,$criterianoasterisk) = $2; 
	chop($criterianoasterisk) if ($criteria =~ /\*$/);
	
	$command ="nosave\nquery find $object where ($criteria)\n";
	foreach my $delete (@deletes) {
	  $command .= "eedit -D $delete\n";
	}
	$command .= "show -a -T -f $filepath\nquit\n";
      }
    }
    # queries requiring tag inclusion
    elsif (/\{(\S+)\}$/) { 
      $include = $1;
      # embedded follow statement
      if ($include =~ /\+$/) {
	  chop $include;
	  $follow = $include;
      }
      print LOG "add $include tags\n";
      if ($follow) { print LOG "... and follow $follow\n";}
      
      if (/^\S+\s+\S+\s+(\S+)\s+\{/) {  
	  $object = $1; 
	  ($criteria,$criterianoasterisk) = ""; 
	  
	  unless ($follow) {
	      $command = "nosave\nquery find $object\n";
	      $command .= "show -a -T $include -f $filepath\nquit\n";
	  } 
	  else {
	      $command = "nosave\nquery find $object\n";
	      $command .= "show -a -T $include -f $filepath\n";
	      $command .= "follow $include\n";
	      $command .= "show -a -T -f $extrafile\nquit\n";
	  }
	  
      } 
      elsif (/^\S+\s+\S+\s+(\S+)\s+(\S+.+)\s+\{/) {
	  $object   = $1; 
	  ($criteria,$criterianoasterisk) = $2; 
	  chop($criterianoasterisk) if ($criteria =~ /\*$/);
	  
	  unless ($follow) {
	      $command ="nosave\nquery find $object where ($criteria)\n";
	      $command .= "show -a -T $include -f $filepath\nquit\n";
	  }
	  else {
	      $command ="nosave\nquery find $object where ($criteria)\n";
	      $command .= "show -a -T $include -f $filepath\n";
	      $command .= "follow $include\n";
	      $command .= "show -a -T -f $extrafile\nquit\n";
	  }
      }
  }
    # simple query
    else { 
      if (/^\S+\s+\S+\s+(\S+)$/) {
	$object=$1; 
	($criteria,$criterianoasterisk) = ""; 
	if ($object eq "DNA") {
	  $command   = "nosave\nquery find Sequence\nfollow DNA\n";
	  $command  .= "show -a -T -f $filepath\nquit\n";
	}
	else {
	  $command   = "nosave\nquery find $object\n";
	  $command  .= "show -a -T -f $filepath\nquit\n";
	}
      }
      elsif (/^\S+\s+\S+\s+(\S+)\s+(\S+.+)$/) {
	$object   = $1; 
	$criteria = $2;
	($criteria,$criterianoasterisk) = $2;
	chop($criterianoasterisk) if ($criteria =~ /\*$/);
	
	if ($object eq "DNA") {
	  $command   = "nosave\nquery find Sequence where ($criteria)\nfollow DNA\n";
	  $command  .= "show -a -T -f $filepath\nquit\n";
	}
	else {
	  $command   ="nosave\nquery find $object where ($criteria)\n";
	  $command  .= "show -a -T -f $filepath\nquit\n";
	}
      }
    }
    
    # dump out from ACEDB
    print "\nFilename: $filepath\n";
    print "Command : 'find $object $criteria' in $dbname\n";

    print LOG &runtime, ": Dumping class $object from database\n";

    open (TACE,"| $exe");
    print TACE $command;
    close TACE;


    # Check file was made (if class is empty no output file will be made from tace, and downstream
    # checks/processing can be skipped)
    if(-e $filepath){
      print "Made file $filepath\n";

      # cat extra data into file if needed and clear up
      if ($follow) {
	system ("cat $extrafile >> $filepath");
	system ("rm -f $extrafile");
	undef ($follow);
      }

    # process database dumps to add database names into time stamps
      &process_ace_file($filepath,$dbname);
    }     
    
    next;
  }
  close CONFIG;
}

########################################################################
# Process ace files to change 'wormpub' in timestamps to database name #
########################################################################

sub process_ace_file{
    my $filename = shift;
    my $database = shift;

    print LOG &runtime, ": Changing timestamp information to contain $database\n";
    
    open (INPUT, "<$filename")            || die "Couldn't open $filename for reading\n";
    open (OUTPUT, ">${filename}.tmpfile") || die "Couldn't write to tmpfile\n";

    
    while (<INPUT>) {

      if (/\"\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:.*?_.*?\"/){
	s/(\"\d{4}-\d{2}-\d{2}_\d{2}:\d{2}:.*?)_.*?\"/$1_$database\"/g;
	print OUTPUT "$_";
      }
      else{
	print OUTPUT "$_";
      }
    }
    close(INPUT);
    close(OUTPUT);
    my $status = move("${filename}.tmpfile", "$filename");
    print "ERROR: Couldn't move file: $!\n" if ($status == 0);
  }


#######################################################################
# Help and error trap outputs                                         #
#######################################################################
 
sub usage {
    my $error = shift;
    
    if ($error eq "Bad_database_name") {
	# Single database mode: database name not recognised
	print "\nNo database of this name exists ('$db')\n";
	print "Check database names in config file\n\n";
	exit(0);
    }
    elsif ($error eq "Help") {
	# Normal help menu
	system ('perldoc',$0);
	exit (0);
    }
    elsif ($error eq "Debug") {
	# No debug bod named
	print "You haven't supplied your name\nI won't run in debug mode until i know who you are\n";
        exit (0);
    }

}

__END__

=pod

=head2   NAME - make_acefiles.pl

=head1 USAGE

=over 4

=item make_acefiles.pl [-options]

=back

make_acefiles.pl is the first step in a WS build.

to write a number of .acefiles which are data dumps from the respective
ACEDB databases which make up the full C.elegans release (i.e. camace,
stlace, geneace, brigace, citace & cshace)

make_acefiles.pl MANDATORY arguments:

=over 4

=item none, 

=back

make_acefiles.pl OPTIONAL arguments:

=over 4

=item -db text, single database mode. Only dumps acefiles for the named database

=item -debug, verbose report (e-mails dl1 only)

=item -test, builds new acefiles into test build area in ~wormpub/TEST_BUILD

=item -help, this help 

=back

make_acefiles.pl DEFAULT behaviour:

make_acefiles.pl will write .acefiles for each of the primary databases to the
/wormsrv2/wormbase/ directories camace, stlace, briggsae, caltech, csh,
emblace, and geneace.

=head1 REQUIREMENTS

=back

This script must run on a machine which can see the /wormsrv2 disk unless in
test mode.

This script requires the following data files:

=over 4

=item autoace_config

=back

This script is called by the following scripts:

=over 4

=item autoace_minder

=item dev_minder

=back

=head1 AUTHOR

=over 4

=item Dan Lawson (dl1@sanger.ac.uk)

=back

=cut


