#!/usr/local/bin/perl5.8.0 -w
#
# EMBLUpdate.pl
#
# This is a script to dump out all genomic canonical clones in embl format.
# Compare them to those submitted to embl in the past and if necessary:
#     * Update the record on file
#     * Place the new records in ~wormpub/analysis/TO_SUBMIT ready for submission.
#     and
#     * Update all necessary files on disk.
#
# Last updated by: $Author: pad $     
# Last updated on: $Date: 2006-07-17 14:39:19 $      

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use File::Copy;
use Storable;

#command line options
my ($test ,$debug ,$help, $file, $suball, $store);
GetOptions (
	    "debug=s"   => \$debug,
	    "file=s"    => \$file,
	    "help"      => \$help,
	    "test"      => \$test,
	    "store:s"   => \$store,
	    "suball"    => \$suball,
	   );

#Build wormbase storable
my $wormbase;
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

&usage("Help") if ($help); # help page
&usage("Debug") if ((defined $debug) && ($debug eq "")); # no debug name
$log->write_to("You have selected to re-submit all Genomic clones\n-------------------------------------------------\n\n") if ($suball);

$log->log_and_die("$file does not exist\n") if (! -e $file); # valid file specified?

# Set up top level base directory which is different if in test mode
my $basedir   = $wormbase->wormpub;

# define global variables
my ($compare,$cosmid, $cosnum, $date, $cossub);

# Timeing [010619 dl] Horrible hack for date stamp - needs rewriting !!!!
my ($sec,$min,$hour,$mday,$mon,$year);
($sec,$min,$hour,$mday,$mon,$year)=localtime(time);
$mon++;
if ($mon < 10) {$mon="0".$mon;}
if ($mday < 10) {$mday="0".$mday;}
$year = $year-100;
my $datestamp = "0" . $year . $mon . $mday;

###################################################################################
#                               Main body of script                               #
###################################################################################

open(IN,"<$file") || die "Could not open $file\n";
while (<IN>) {
  # if detect the beginning of an entry start a new file
  if (/^ID\s+\w+\;\s+XXX;/) {
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

    &getCurrentDate_dir;       # Extract the current date dir for each clone.

    if (!$suball) {
      &CreateCosmidDirStructure; # If there is a new genomic clone, create dir structure.
      &Check_current_version;    # Is there a current.embl file for this sequence.

      my $compare = &detectdifference($cosmid,$date);

      if ($compare==1) {
	print "$cosmid differs\n";
	$log->write_to("$cosmid differs\n");
	$cossub++;
	&PrepareCosmidSubmission;
      }
      else {
	$log->write_to("Cosmid $cosmid is up to date\n\n");
	unlink "temp$$";
	next;
      }
    }
    elsif ($suball) {
      &suball;
    }
  }
}
close(IN);
$log->write_to("\nThere were $cosnum record(s) checked in this round of EMBL updating\n\n$cossub need re-submitting\n\n\n\n");
$log->mail();
exit(0);

##################################################################################
#                                   Subroutines                                  #
##################################################################################

sub getCurrentDate_dir   {
  # Get the current datedir for the cosmid
  open(CURRENT, "grep ^$cosmid\/  $basedir/analysis/cosmids/current.versions  |"); 
  while (<CURRENT>) {
    if (/^\S+\/(\d+)/) {
      $date=$1;
      $log->write_to("The current date dir for cosmid $cosmid id $date.\n");
    }
    else  { 
      $log->write_to ("** WARNING - NO CURRENT VERSIONS ENTRY FOR $cosmid\n"); 
      next;
    }
  }
  close CURRENT;
}

sub CreateCosmidDirStructure {
  # If there never has been a cosmid entry so far then may wish to make one#
  if (!-e "$basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.embl") {
    print "This cosmid has not been submitted before\n";
    print "Shall I make a $cosmid.embl file and make this the $cosmid.current.file and submit? (y or n)\n";
    my $answer=<STDIN>;
    if ($answer eq "y\n") {
      $log->write_to("\nCreating $cosmid.embl and $cosmid.current.file\n");
      my $status = copy("temp$$", "$basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.embl");
      $log->write_to("ERROR: Couldn't copy file: $!\n") if ($status == 0);
      system "ln -s $basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.embl $basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.current.embl";
      copy("temp$$", "$basedir/analysis/TO_SUBMIT/$cosmid.$datestamp.embl");	 
      $log->write_to("ERROR: Couldn't copy file: $!\n") if ($status == 0);
    }
    else {
      print "** WARNING - COSMID $cosmid DOES NOT HAVE A $cosmid.embl FILE AND HAS NOT BEEN SUBMITTED\n" if ($debug);
      $log->write_to ("** WARNING - COSMID $cosmid DOES NOT HAVE A $cosmid.embl FILE AND HAS NOT BEEN SUBMITTED\n\n");
      next;
    }
  }
}

sub Check_current_version {
  # Check to see if there are differences between the dumped record and the record on file.
  if (-e "$basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.current.embl"){ 
    print "Current file for $cosmid found (number $cosnum)\n" if ($debug);
    $log->write_to("Current file for $cosmid found (number $cosnum)\n");
  } 
  else {
    $log->write_to("\nNO CURRENT FILE FOR COSMID $cosmid\n");
    $log->write_to("Program has stopped. Cosmids which need resubmission up until now have been placed\n"); 
    $log->write_to("in ~wormpub/analysis/TO_SUBMIT\n"); 
    $log->write_to("You may need to just run make_current.embl to sort this out\n");
    $log->log_and_die("There has been an error processing $cosmid.....please investigate!\n");
    unlink "temp$$";
    exit;
  }
}

sub PrepareCosmidSubmission {
  # If there is a difference then copy the temp file to the appropriate dir.
  my $status = copy("temp$$", "$basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.$datestamp.embl");
  $log->write_to("ERROR: Couldn't copy file: $!\n") if ($status == 0);
  # Make a symbolic link to cosmid.current.embl
  system "ln -fs $basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.$datestamp.embl $basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.current.embl";
  # Report the change and copy to the TO_SUBMIT directory.
  $status = move("temp$$", "$basedir/analysis/TO_SUBMIT/$cosmid.$datestamp.embl");
  $log->write_to("ERROR: Couldn't move file: $!\n") if ($status == 0);
  $log->write_to("$cosmid should be resubmitted and is in TO_SUBMIT.\n");
  $log->write_to("\nOUTPUT: .embl files can be found  in ~wormpub/analysis/TO_SUBMIT\n\n");
}

sub detectdifference {
  # Detect whether the new file is different from the old one using diff
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

sub suball {
  $log->write_to("Copying $basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.current.embl\n\n");
  my $status = copy("$basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.current.embl", "$basedir/analysis/TO_SUBMIT/$cosmid.$datestamp.embl");
  $cossub++;
  $log->write_to("ERROR: Couldn't copy file: $!\n") if ($status == 0);
  #system "cp $basedir/analysis/cosmids/$cosmid/$date/embl/$cosmid.current.embl", "$basedir/analysis/TO_SUBMIT/$cosmid.$datestamp.embl";
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

__END__
