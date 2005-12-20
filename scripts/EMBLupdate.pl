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
# Last updated by: $Author: pad $     
# Last updated on: $Date: 2005-12-20 14:19:51 $      


use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use File::Copy;
use Storable;

##############
# variables  #
##############

my $name;
my $first;
my $second;
my $wormbase;
my $store;

#command line options
my ($test ,$debug ,$help, $file);
GetOptions (
	    "debug=s"   => \$debug,
	    "file=s"    => \$file,
	    "help"      => \$help,
	    "test"      => \$test,
	    "store:s"  =>  \$store,
	   );


#Build wormbase storable
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}


# establish log file.
my $log = Log_files->make_build_log($wormbase);

#who will receive log file?
my $maintainers = "All";

#Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

###############################
# Variables using Wormbase.pm #
###############################

# help page
&usage("Help") if ($help);

# no debug name
&usage("Debug") if ((defined $debug) && ($debug eq ""));

# valid file specified?
die "$file does not exist\n\n" if (! -e $file);

# Set up top level base directory which is different if in test mode
my $basedir   = $wormbase->wormpub;

#&create_log_files;

$log->write_to("\nList of ERRORS\n-----------------\n");


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
    open(CURRENT, "grep ^$cosmid\/  $basedir/analysis/cosmids/current.versions  |"); 
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
    
    if (!-e "$basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.embl") {
      print "This cosmid has not been submitted before\n";
      print "Shall I make a $cosmid.embl file and make this the $cosmid.current.file and submit? (y or n)\n";
      my $answer=<STDIN>;
      if ($answer eq "y\n") {
	my $status = copy("temp$$", "$basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.embl");
	$log->write_to("ERROR: Couldn't copy file: $!\n") if ($status == 0);
	system "ln -s $basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.embl $basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.current.embl";
	copy("temp$$", "$basedir/analysis/TO_SUBMIT/$cosmid.$datestamp.embl");	 
	$log->write_to("ERROR: Couldn't copy file: $!\n") if ($status == 0);
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

    if (-e "$basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.current.embl"){ 
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
      my $status = copy("temp$$", "$basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.$datestamp.embl");
      $log->write_to("ERROR: Couldn't copy file: $!\n") if ($status == 0);
      system "ln -fs $basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.$datestamp.embl $basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.current.embl";
      $status = move("temp$$", "$basedir/analysis/TO_SUBMIT/$cosmid.$datestamp.embl");
      $log->write_to("ERROR: Couldn't move file: $!\n") if ($status == 0);
      print "$cosmid should be resubmitted and is in TO_SUBMIT.\n";
    } 
    else {
      print "Cosmid $cosmid is up to date\n";
      unlink "temp$$";
    } 
    next;
  }
  
}
$log->write_to("\nOUTPUT: .embl files can be found  in ~wormpub/analysis/TO_SUBMIT\n\n");
close(IN);
$log->mail("$maintainers");

exit(0);

#----------------------------
# Detect whether the new file 
# is different from the old one
# using diff
#
sub detectdifference {		 
  my $difference = 0;
  open (DIFFOUT,"diff temp$$ $basedir/analysis/cosmids/$_[0]/$_[1]/embl/$cosmid.current.embl |");
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

