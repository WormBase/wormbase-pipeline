#!/usr/local/bin/perl5.8.0 -w
#
# EMBLUpdate.pl
#
# by Dan Lawson, Allesandro, and Steve Jones
#
# [010619 dl] Horrible hack for date stamp - needs rewriting !!!!
# [ag3 0006] Re-writing of the original Steve Jones script
#
# This program will:-
#	
#	1) Check to see if the new version is different to before
#
#	2) If different, will write the new version into the cosmid/embl
#	 e.g. ZK1058.951111.embl and symbolically link this to ZK1058.current.embl
#		 	
#       3) Copy the version to ~wormpub/analysis/TO_SUBMIT.  
#
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2003-12-04 14:07:00 $      


use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use File::Copy;

##############
# variables  #
##############

my $maintainers = "All";
my $rundate     = &rundate;
my $name;
my $log;
my $first;
my $second,

my ($test ,$debug ,$help, $file);
GetOptions (
	    "debug=s"   => \$debug,
	    "file=s"    => \$file,
	    "help"      => \$help,
	    "test"      => \$test);

# help page
&usage("Help") if ($help);

# no debug name
&usage("Debug") if ((defined $debug) && ($debug eq ""));

# valid file specified?
die "$file does not exist\n\n" if (! -e $file);

# assign $maintainers if $debug set
($maintainers = $debug . '\@sanger.ac.uk') if ($debug);


# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir   = "/wormsrv2";
$basedir      = glob("~wormpub")."/TEST_BUILD" if ($test); 


&create_log_files;


print LOG "EMBLupdate.pl start at ",&runtime,"\n";


my ($sec,$min,$hour,$mday,$mon,$year);
($sec,$min,$hour,$mday,$mon,$year)=localtime(time);
$mon++;
if ($mon < 10) {$mon="0".$mon;}
if ($mday < 10) {$mday="0".$mday;}
$year = $year-100;

my $datestamp = "0" . $year . $mon . $mday;

my($cosmid, $cosnum);

open(IN,"<$file") || die "Could not open $file\n";

while (<IN>) {
  # if detect the beginning of an entry start a new file
  if (/^ID\s+CE(\S+)/) {
    $cosmid=$1;
    open(TEMP,">temp$$");
    $cosnum++;
  }
  # for some clones with long names the clone name cannot be got from the ID
  if (/^DE\s+Caenorhabditis elegans\s+\S+\s(\S+)$/) { $cosmid=$1; }

  # printout the whole file
  print TEMP $_;
  
  # if detect the end of an entry, close the temp file and begin
  if (/^\/\/$/) { 
    close TEMP;

    #------------------------------------
    # Get the current date for the cosmid
    #
    open(CURRENT, "grep ^$cosmid\/  /nfs/disk100/wormpub/analysis/cosmids/current.versions  |"); 
    my $date;
    while (<CURRENT>) {
      if (/^\S+\/(\d+)/) { 
	$date=$1; 
      } 
      else  { 
	print "** WARNING - NO CURRENT VERSIONS ENTRY FOR $cosmid\n"; 
	next;
      }
    }
    close CURRENT;
    
    #------------------------------------------------------------------------
    # If there never has been a cosmid entry so far then may wish to make one
    #
    
    if (!-e "/nfs/disk100/wormpub/analysis/cosmids/$cosmid/$date/embl/$cosmid.embl") {
      print "This cosmid has not been submitted before\n";
      print "Shall I make a $cosmid.embl file and make this the $cosmid.current.file and submit? (y or n)\n";
      my $answer=<STDIN>;
      if ($answer eq "y\n") {
	my $status = copy("temp$$", "/nfs/disk100/wormpub/analysis/cosmids/$cosmid/$date/embl/$cosmid.embl");
	print LOG "ERROR: Couldn't copy file: $!\n" if ($status == 0);
	system "ln -s /nfs/disk100/wormpub/analysis/cosmids/$cosmid/$date/embl/$cosmid.embl /nfs/disk100/wormpub/analysis/cosmids/$cosmid/$date/embl/$cosmid.current.embl";
	copy("temp$$", "/nfs/disk100/wormpub/analysis/TO_SUBMIT/$cosmid.$datestamp.embl");	 
	print LOG "ERROR: Couldn't copy file: $!\n" if ($status == 0);
      } 
      else {
	print "** WARNING - COSMID $cosmid DOES NOT HAVE A $cosmid.embl FILE AND HAS NOT BEEN SUBMITTED\n";
	next;
      }
    }
    
    
    #------------------------------------------------------------------
    # Is there a current.embl file for this sequence in the appropriate 
    # place - if not quit now and sort this out 	
    #

    if (-e "/nfs/disk100/wormpub/analysis/cosmids/$cosmid/$date/embl/$cosmid.current.embl"){ 
      print "Current file for $cosmid found (number $cosnum)\n";
    } 
    else {
      print "NO CURRENT FILE FOR COSMID $cosmid\n";
      print "Program has stopped. Cosmids which need resubmission up until now have been placed\n"; 
      print "in ~wormpub/analysis/TO_SUBMIT\n"; 
      print "You may need to just run make_current.embl to sort this out\n";
      unlink "temp$$";
      exit;
    }


    #------------------------------------------------
    # Is this different from the current.embl version 
    #
    
    my $compare = &detectdifference($cosmid,$date);
    if ($compare==1) {print "$cosmid differs\n";} 
    
    #--- -----------------------------------------------------------------
    # If there is a difference then copy the temp file to the approp. dir
    # and make a symbolic link to cosmid.current.embl
    # Report the change and place in TO_SUBMIT directory. 

    if ($compare==1) {
      my $status copy("temp$$", "/nfs/disk100/wormpub/analysis/cosmids/$cosmid/$date/embl/$cosmid.$datestamp.embl");
      print LOG "ERROR: Couldn't copy file: $!\n" if ($status == 0);
      system "ln -fs /nfs/disk100/wormpub/analysis/cosmids/$cosmid/$date/embl/$cosmid.$datestamp.embl /nfs/disk100/wormpub/analysis/cosmids/$cosmid/$date/embl/$cosmid.current.embl";
      $status = move("temp$$", "/nfs/disk100/wormpub/analysis/TO_SUBMIT/$cosmid.$datestamp.embl");
      print LOG "ERROR: Couldn't move file: $!\n" if ($status == 0);
      print "$cosmid should be resubmitted and is in TO_SUBMIT.\n";
    } 
    else {
      print "Cosmid $cosmid is up to date\n";
      unlink "temp$$";
    } 
    next;
  }
  
}
close(IN);
print LOG "$0 finished at ",&runtime, "\n";
close LOG;
&mail_maintainer("EMBLupdate",$maintainers,"$log");

exit(0);

#----------------------------
# Detect whether the new file 
# is different from the old one
# using diff
#
sub detectdifference {		 
  my $difference = 0;
  open (DIFFOUT,"diff temp$$ /nfs/disk100/wormpub/analysis/cosmids/$_[0]/$_[1]/embl/$cosmid.current.embl |");
  while (<DIFFOUT>) {
    if ($_ ne "") {
      $difference=1;
    }
  }
  return($difference);					   
  close(DIFFOUT);
}

sub usage {
  my $error = shift;
  
  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
  elsif ($error eq "Debug") {
    # No debug bod named
    print "You haven't supplied your name\nI won't run in debug mode
         until i know who you are\n";
    exit (0);
  }
}


######################################################################################


sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch $basedir/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate = &rundate;
  $log        = "$basedir/logs/$script_name.$rundate.$$";

  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name\n";
  print LOG "started at ",&rundate,"\n";
  print LOG "=============================================\n";
  print LOG "\n";

}






