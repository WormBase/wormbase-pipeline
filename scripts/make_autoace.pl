#!/usr/local/bin/perl5.8.0 -w
#
# make_autoace
# dl1/ag3 & others
#
# Usage : make_autoace [-options]
#
# This makes the autoace database from its composite sources.
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2004-11-24 15:06:56 $

use strict;
use lib  $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use POSIX qw(:signal_h :errno_h :sys_wait_h);
use Getopt::Long;
use Cwd;
use File::Copy;
use File::Path;
use Log_files;


#################################
# Command-line options          #
#################################

our ($help, $debug, $test);

GetOptions ("help"         => \$help,
            "debug=s"      => \$debug,
	    "test"         => \$test);


# Display help if required
&usage("Help") if ($help);


###################################
# Check command-line arguments    #
###################################


# Exit if no database specified
# database/file paths and locations
my $basedir     = glob("~wormpub");
$basedir       .= "/TEST_BUILD" if ($test); 

my $wormbasedir = "$basedir/wormbase";
my $autoacedir  = "$basedir/autoace";
my $stlacedir   = "$basedir/stlace";
my $camacedir   = "$basedir/camace";
my $configfile  = "$basedir/autoace_config/autoace.config";
my $dbpath = $autoacedir;
my $CWD = cwd;

# make log files
my $log = Log_files->make_build_log($debug);

##################
# misc variables # 
##################

my $WS_current;
my $WS_version;

# need to change this if in test mode
if($test){
  $WS_current       = "666";
  $WS_version       = "WS666";
}
else{
  $WS_current    = &get_wormbase_version;
  $WS_version    = &get_wormbase_version_name;
}
my @filenames; # for storing contents of autoace_config

my $tace   = &tace;
my $giface = &giface;
my $errors = 0; # for tracking system call related errors

# start doing stuff
#&mail_reminder;

#Set up correct database structure if it doesn't exist
#&createdirs;	

# Parse config file
#&parseconfig;

# Re-initialize the database
# Re-initializing database is more complex as the .acefiles will be spread over a number of
# directories - use the config file to find them all
#reinitdb();

# remove temp genes
#&rmtempgene();

# Read in the physical map and make all maps
#&physical_map_stuff();

# Make the chromosomal links
#&makechromlink();

&check_make_autoace;

$log->mail;
exit (0);



##############################################################
#
# Subroutines
#
##############################################################


###################################################
# Subroutine for writing to a given database      

sub DbWrite {
    my ($command,$exec,$dir,$name)=@_;
    open (WRITEDB,"| $exec $dir >> STDERR") or $log->log_and_die("$name DbWrite failed\n");
    print WRITEDB $command;
    close WRITEDB;
}


###################################################
# Get time coordinates

sub GetTime {
    my ($SECS,$MINS,$HOURS,$DAY,$MONTH,$YEAR)=(localtime)[0,1,2,3,4,5];
    if ($MINS=~/^\d{1,1}$/) {
	$MINS="0"."$MINS";
    }
    my $REALMONTH=$MONTH+1;
    my $REALYEAR=$YEAR+1900;
    my $NOW = "$DAY-$REALMONTH-$REALYEAR:$HOURS:$MINS";
    return $NOW;
} 


################################################
# Remove temp_gene sequences from stlace, camace 

sub rmtempgene {
  
  $log->write_to( &runtime." : starting rmtempgene subroutine\n");
  my $camace  = "$camacedir";
  my $stlace  = "$stlacedir";
  my $command = "query find elegans_CDS method = hand_built\nkill\nsave\nquit\n";
  &DbWrite($command,$tace,$camace,"CamAce");
  &DbWrite($command,$tace,$stlace,"StlAce");
  my $command2 = "query find elegans_CDS temp*\nkill\nsave\nquit\n";
  &DbWrite($command2,$tace,$camace,"CamAce");
  &DbWrite($command2,$tace,$stlace,"StlAce");
  $log->write_to(&runtime." : Finished.\n\n");
}


############################################
# Create directories

sub createdirs {

  $log->write_to( &runtime. ": starting createdirs subroutine\n");

  my $chromes  = "$dbpath/CHROMOSOMES";
  my $db       = "$dbpath/database";		
  my $new      = "$db/new";
  my $touch    = "$db/touched";
  my $ace      = "$dbpath/acefiles";
  my $rel      = "$dbpath/release";
  my $wspec    = "$dbpath/wspec";
  my @dirarray    = ("$dbpath","$chromes","$db","$new","$touch","$ace","$rel","$wspec");
  
  foreach (@dirarray) {	
    my $present_dir = $_;
    if (-d $present_dir) {
      $log->write_to( "\t** $present_dir - already present\n");
      print "** Skipping $present_dir - already present\n";
      next;
    }
    else {
      $log->write_to("making $present_dir\n");
      mkpath($present_dir);
    }			
  }
  $log->write_to( &runtime, ": Finished\n\n");
}


###################################################
# Parses the lists from the config files     

sub parseconfig {

  $log->write_to( &runtime. ": starting parseconfig subroutine\n");
  my ($filename,$dbname);
  open(CONFIG,"$configfile");
  while(<CONFIG>) {

    # some formating, remove leading and trailing whitespace
    s/^\s+//;
    s/\s+$//;   
    
    # next if comment line
    next if (/^\#/ || /^$/);
    
    # parse database information
    if (/^P\s+(\S+)\s+(\S+)$/) {
      $dbname = $1;
      #	    $dbdir  = $2;
      #	    $targetdir="$wormbasedir"."/$dbname";
      next;
    }
    
    # parse file name
    if (/^\S+\s+(\S+)/) {
      $filename = $1;
    }
    
    # next if no filename parsed
    if (!defined $filename) {
      $log->write_to( "ERROR: Failed to parse filename ..\n");
      next;
    }
    
    # check that file exists before adding to array and is not zero bytes
    if (-e "$wormbasedir"."/$dbname/"."$filename") {
      if (-z "$wormbasedir"."/$dbname/"."$filename") {
	$log->write_to( "ERROR: file $wormbasedir/$dbname/$filename is zero bytes !\n");
	$errors++;
      }
      else{
	push (@filenames,"$wormbasedir"."/$dbname/"."$filename");
	$log->write_to( "* Parse config file : file $wormbasedir/$dbname/$filename noted ..\n");
      }
    } 
    else {
      $log->write_to( "ERROR: file $wormbasedir/$dbname/$filename is not existent !\n");
      $errors++;
      next;
    }
    
  }
  close(CONFIG);
  $log->write_to( &runtime. ": Finished\n\n");
}




###################################################
# Re-initialize the database
# cleans and re-initalizes the ACEDB residing in $dbpath
# then parses .ace files in @filenames

sub reinitdb {

  $log->write_to( &runtime. ": Starting reinitdb subroutine\n");

  if( -e "$dbpath/wspec/models.wrm" ) {
    $log->write_to("$dbpath/wspec/model already exists . . \n");
    if (-e "$dbpath/database/lock.wrm") {
      $log->log_and_die( "*Reinitdb error - lock.wrm file present..\n");
    }

    &delete_files_from("$dbpath/database",".\.wrm","-");
  }
  else {
    $log->write_to("Create a new database from scratch at $dbpath");

    # CVS checkout the wspec dir in to there
    $log->write_to("running CVS checkout of wspec dir\n");
    &run_command("cd $dbpath;cvs checkout -d wspec wormbase/wspec");
    $log->log_and_die("cvs checkout failed\n") unless (-e "$dbpath/wspec/models.wrm");
  }

  my $command = "y\n";
  $log->write_to( &runtime. ": reinitializing the database\n");
  &DbWrite($command,$tace,$dbpath,"ReInitDB");


  foreach my $filename (@filenames) {
    my $command = "pparse $filename\nsave\nquit\n";
    if (-e $filename) {
      my $runtime = &runtime;
      $log->write_to( "* Reinitdb: started parsing $filename at $runtime\n");
      my ($tsuser) = $filename =~ (/^\S+\/(\S+)\_/);
      &DbWrite($command,"$tace -tsuser $tsuser",$dbpath,"ParseFile");
    }
    else {
      $log->write_to( "* Reinitdb: $filename is not existent - skipping ..\n");
      next;
    }
  }
  $log->write_to( &runtime. ": Finished.\n\n");

}


###################################################
# Alan Coulson maintains a ContigC
# database in ~cemap for the physical map.
# This is dumped in file ~cemap/cen2hs.ace 
# (makecen2hs.pl nightly cron job on rathbin)

sub physical_map_stuff{

  $log->write_to( &runtime. ": starting physical_map_stuff subroutine\n");

# This data is permanently in geneace now, the previous file was no longer
# being updated once Alan left
  
  # now make the maps
  my $command = "gif makemaps -all\nsave\ngif makemaps -seqclonemap $dbpath/acefiles/seqclonemap.ace\n";
  $command .= "pparse $dbpath/acefiles/seqclonemap.ace\nsave\nquit\n";
  &DbWrite($command,$giface,$dbpath,"MakeMaps");

  $log->write_to( &runtime. ": finished\n\n");
}


###################################################
# Set the date correctly in displays.wrm

sub setdate {
  my @t   = localtime ; while ($t[5] >= 100) { $t[5] -= 100 ; }
  my $dat = sprintf "%02d\/%02d\/%02d", $t[3], $t[4]+1, $t[5] ;
  my $status = move("$dbpath/wspec/displays.wrm", "$dbpath/wspec/displays.old");
  print "ERROR: Couldn't move file: $!\n" if ($status == 0);


  open(FILE,"$dbpath/wspec/displays.old") or do { $log->write_to( "failed to open $dbpath/wspec/displays.old\n"); return 1;};
  open(NEWFILE,">$dbpath/wspec/displays.wrm") or do { $log->write_to( "failed to open $dbpath/wspec/displays.wrm\n"); return 1;};
  while (<FILE>) {
    if (/^_DDtMain/) {
      print NEWFILE "_DDtMain -g TEXT_FIT -t \"C.elegans database $dat\"  -w .43 -height .23 -help acedb\n";
    } 
    else {
      print NEWFILE $_;
    }
  }
  close(FILE);
  close(NEWFILE);
  unlink "$dbpath/wspec/displays.old" or $log->write_to( "ERROR: Couldn't unlink file: $!\n");
}





#############################
# Make chromosomal links    #
#############################

sub makechromlink {

  $log->write_to( &runtime. ": starting makechromlink subroutine\n");
  my $chromlink_file = "$dbpath/acefiles/chromlinks.ace";

  unlink "$chromlink_file" or $log->write_to( "ERROR: Couldn't unlink file $chromlink_file: $!\n");
  my $command = "makeChromLinks.pl -out $chromlink_file";
  $command = "makeChromLinks.pl -test -out $chromlink_file"  if ($test);
  &run_command("$command", $log);

  if (-z "$chromlink_file") {
    $log->write_to( "*Makechromlink: $chromlink_file has ZERO size\n");
    return;
  } 
  else {
    my $command = "pparse $chromlink_file\nsave\nquit\n";
    &DbWrite($command,"$tace -tsuser make_autoace",$dbpath,"MakeChromLinks");
  }
  
  $log->write_to( &runtime. ": Finished.\n\n");
}



#########################################
#  Help and usage
#########################################
sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

sub check_make_autoace {
  local (*BUILDLOG);

  $log->write_to( &runtime. ": Entering check_make_autoace subroutine\n");

  my $log_file = $log->get_file;

  print "Looking at log file: $log_file\n";  
  print "Open log file $log_file\n";
  
  my ($parsefile,$parsefilename);
  my $builderrors = 0;

  open (BUILDLOG, "<$log_file") || $log->log_and_die("Couldn't open $log_file out\n");
  while (<BUILDLOG>) {
    if (/^\* Reinitdb: started parsing (\S+)/) {
      $parsefile = $1;
    }
    if ((/^\/\/ objects processed: (\d+) found, (\d+) parsed ok, (\d+) parse failed/) && ($parsefile ne "")) {
      my $object_count = $1;
      my $error_count  = $3;
      (printf "%6s parse failures of %6s objects from file: $parsefile\n", $error_count,$object_count);
      if ($error_count > 0) {
	$parsefilename = $parsefile;
	$parsefilename =~ s/$basedir//;
	$log->write_to(sprintf("%6s parse failures of %6s objects from file: $parsefilename\n", $error_count,$object_count));
	$builderrors++;
      }
      ($parsefile,$parsefilename) = "";
    }
  }
  close BUILDLOG;


  # look for objects with no Gene tag
  $log->write_to( "\n". &runtime. ": Looking for CDSs, Transcripts, Pseudogenes with no Gene tag\n");
  my $db = Ace->connect(-path=>$autoacedir, -program =>$tace) || $log->log_and_die("Connection failure: ". Ace->error);

  my @genes= $db->fetch(-query=>'find worm_genes NOT Gene');
  if(@genes){
    foreach (@genes){
      $log->write_to( "ERROR: $_ has no Gene tag, please add valid Gene ID from geneace\n");
      $builderrors++;
    }
  }
  $db->close;

  $log->write_to( &runtime. ": Finished subroutine\n\n");

  return ($builderrors);
}

sub allcmid
  {
    # make the allcmid file needed for the farm
    my $command = "query find genome_sequence\nDNA -f /wormsrv2/autoace/allcmid\nquit\n";
    open (WRITEDB, "| $tace $basedir/autoace ") || $log->log_and_die("Couldn't open pipe to autoace\n");
    print WRITEDB $command;
    close (WRITEDB);
  }

sub mail_reminder 
  {
    open (EMAIL,  "|/bin/mailx -s \"WormBase build reminder\" \"wormbase\@sanger.ac.uk\" ");
    print EMAIL "Dear builder,\n\n";
    print EMAIL "You have just run autoace_minder.pl -build.  This will probably take 5-6 hours\n";
    print EMAIL "to run.  You should therefore start work on the blast pipeline. So put down that\n";
    print EMAIL "coffee and do some work.\n\n";
    print EMAIL "Yours sincerely,\nOtto\n";
    close (EMAIL);
  }
__END__

=pod

=head1 NAME - make_autoace.pl

=head2 USAGE

make_autoace.pl  makes the autoace database in a given directory.

make_autoace.pl  arguments:

=over 4

=item *

-database path_of_database => location of the new autoace directory

=item *

-buildautoace => rebuild the database from source databases (geneace, camace etc)

=item *

-buildrelease => creates only the distribution (release) files

=item *

-test => runs in test mode and uses test environment in ~wormpub/TEST_BUILD
=back
 

=cut
















