#!/usr/local/bin/perl5.8.0 -w
#
# make_autoace
# dl1/ag3 & others
#
# Usage : make_autoace [-options]
#
# This makes the autoace database from its composite sources.
#
# Last edited by: $Author: krb $
# Last edited on: $Date: 2004-03-02 15:09:36 $

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use POSIX qw(:signal_h :errno_h :sys_wait_h);
use Getopt::Long;
use Cwd;
use File::Copy;


#################################
# Command-line options          #
#################################

our ($help, $debug, $database, $buildautoace, $buildrelease, $log, $test);

GetOptions ("help"         => \$help,
            "debug=s"      => \$debug,
	    "database=s"   => \$database,
	    "buildautoace" => \$buildautoace,
	    "buildrelease" => \$buildrelease,
	    "test"         => \$test);


# Display help if required
&usage("Help") if ($help);

my $maintainers   = "All";

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}


###################################
# Check command-line arguments    #
###################################

&usage("Help") if ((!$buildautoace)&&(!$buildrelease));

# Exit if no database specified
my $dbpath;
my $CWD = cwd;

if(!$database){
  &usage("Database");
}
else{
  if($database =~ /^(\~\w+)\//){ # do we need to expand path if ~path specified?
    $dbpath = glob("$1");
    $dbpath =~ s/\/tmp_mnt//;
    my $filename = "$'";
    $dbpath = "$dbpath"."/"."$filename";
  } 
  elsif($database =~ /^(\w+)/) { # for incomplete paths, expand using CWD
    $dbpath="$CWD"."/"."$database";
  }  
  elsif ($database =~ /\/\w+/) { # else assume path is ok
    $dbpath=$database;
  } 
  else {
    &usage("Help");
  }
}


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

# database/file paths and locations
my $basedir     = "/wormsrv2";
$basedir        = glob("~wormpub")."/TEST_BUILD" if ($test); 

my $wormbasedir = "$basedir/wormbase";
my $autoacedir  = "$basedir/autoace";
my $stlacedir   = "$basedir/stlace";
my $camacedir   = "$basedir/camace";
my $configfile  = "$basedir/autoace_config/autoace.config";

my $tace   = &tace;
my $giface = &giface;
my $errors = 0; # for tracking system call related errors


# Open logfile                                   
&create_log_files;


################################################
# Avoid filling process table with zombies     #
################################################

$SIG{CHLD} = \&REAPER;
sub REAPER {
  my $pid;
  $pid=waitpid(-1,&WNOHANG);
  $SIG{CHLD}=\&REAPER;
}


##################################################	
# Create the directory structure                 #
##################################################

&buildautoace if($buildautoace);
&buildrelease if($buildrelease);


################################################
# Finish and tidy up                           #
################################################

print LOG &runtime, " make_autoace.pl finished\n\n";
close (LOG);

# warn about errors in subject line if there were any
if($errors == 0){
  &mail_maintainer("BUILD REPORT: make_autoace.pl",$maintainers,$log);
}
elsif ($errors ==1){
  &mail_maintainer("BUILD REPORT: make_autoace.pl : $errors ERROR!",$maintainers,$log);
}
else{
  &mail_maintainer("BUILD REPORT: make_autoace.pl : $errors ERRORS!!!",$maintainers,$log);
}

exit (0);



##############################################################
#
# Subroutines
#
##############################################################




################################################
# Build new autoace database
################################################

sub buildautoace{

  #Set up correct database structure if it doesn't exist
  &createdirs;	


  # Parse config file                            
  &parseconfig;
  

  # Re-initialize the database                      
  # Re-initializing database is more complex as the .acefiles will be spread over a number of
  # directories - use the config file to find them all
  &reinitdb();

  
  # remove temp genes
  &rmtempgene();
  
  

  # Read in the physical map and make all maps   
  &physical_map_stuff();


  # Make the chromosomal links                  
  &makechromlink();
    
}


###################################################################


sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; 
  system("touch $basedir/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate = &rundate;
  $log        = "$basedir/logs/$script_name.$WS_version.$rundate.$$";

  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name\n";
  print LOG "started at ",`date`,"\n";
  print LOG "=============================================\n";
  print LOG "\n";

}




###################################################
# Subroutine for writing to a given database      

sub DbWrite {
    my ($command,$exec,$dir,$name)=@_;
    open (WRITEDB,"| $exec $dir >> $log") or do {print LOG "$name DbWrite failed\n";close LOG; die();};
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
  
  print LOG &runtime, ": starting rmtempgene subroutine\n";
  my $camace  = "$camacedir";
  my $stlace  = "$stlacedir";
  my $command = "query find elegans_CDS method = hand_built\nkill\nsave\nquit\n";
  &DbWrite($command,$tace,$camace,"CamAce");
  &DbWrite($command,$tace,$stlace,"StlAce");
  my $command2 = "query find elegans_CDS temp*\nkill\nsave\nquit\n";
  &DbWrite($command2,$tace,$camace,"CamAce");
  &DbWrite($command2,$tace,$stlace,"StlAce");
  print LOG &runtime, ": Finished.\n\n";
}


############################################
# Create directories

sub createdirs {

  print LOG &runtime, ": starting createdirs subroutine\n";

  my $chromes  = "$dbpath/CHROMOSOMES";
  my $db       = "$dbpath/database";		
  my $new      = "$db/new";
  my $touch    = "$db/touched";
  my $ace      = "$dbpath/acefiles";
  my $rel      = "$dbpath/release";
  my $wspec    = "$dbpath/wspec";
  my $pictures = "$dbpath/pictures";
  my @dirarray    = ("$dbpath","$chromes","$db","$new","$touch","$ace","$rel","$wspec");
  my @args1       = ();
  my @args2       = ("/bin/mkdir");
  
  foreach (@dirarray) {	
    my $present_dir = $_;
    if (-d $present_dir) {
      print "** Skipping $present_dir - already present\n";
      next;
    }
    else {
      push (@args1,"$present_dir");
    }			
  }
  my $argsize = scalar (@args1);
  if ($argsize == 0) {
    print "** No new directories to create .. end mkdirectories\n";
    return;
  }
  my $command = "@args2 @args1";
  &run_command ($command);
  foreach my $made_dir (@args1) {
    if ($made_dir =~ /wspec/) {
      print "** Copying wspec from $autoacedir .. \n";
      my $status = copy("$autoacedir/wspec/*", "$wspec/.");
      print LOG "ERROR: Couldn't copy file: $!\n" if ($status == 0);
    }
    if (!-d $made_dir) {
      print " ** mkdir for $made_dir failed ** \n\n";
      print LOG "ERROR: mkdir command failed\n";
      $errors++;
      die(0);
    } 
  }
  &run_command("/bin/ln -s $basedir/geneace/pictures $pictures");

  print LOG &runtime, ": Finished.\n\n";
  return;
}


###################################################
# Parses the lists from the config files     

sub parseconfig {

  print LOG &runtime, ": starting parseconfig subroutine\n";
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
      print LOG "ERROR: Failed to parse filename ..\n";
      next;
    }
    
    # check that file exists before adding to array and is not zero bytes
    if (-e "$wormbasedir"."/$dbname/"."$filename") {
      if (-z "$wormbasedir"."/$dbname/"."$filename") {
	print LOG "ERROR: file $wormbasedir/$dbname/$filename is zero bytes !\n";
	$errors++;
      }
      else{
	push (@filenames,"$wormbasedir"."/$dbname/"."$filename");
      }
    } 
    else {
      print LOG "ERROR: file $wormbasedir/$dbname/$filename is not existent !\n";
      $errors++;
      next;
    }
    
  }
  close(CONFIG);
  print LOG &runtime, ": Finished.\n\n";
}




###################################################
# Re-initialize the database                      
# cleans and re-initalizes the ACEDB residing in $dbpath
# then parses .ace files in @filenames
#
# 011016 : dl  : Added '-f' option to the rm lines. This will ensure that non 
#                wormpub owned files within group worm are deleted without the
#                need for an interactive prompt
# 021025 : dl  : Added single Dbwrite command which parses the first section of
#                a filename to assign the timestamp user string. i.e. a file
#                camace_Sequence.ace will be loaded as camace.

sub reinitdb {

  print LOG &runtime, ": Starting reinitdb subroutine\n";

  &delete_files_from("$dbpath/database/new","*","-") or print LOG "ERROR: Problems removing files from $dbpath/database/new\n";
  &delete_files_from("$dbpath/database/touched","*","-") or print LOG "ERROR: Problems removing files from $dbpath/database/touched\n";


  if (-e "$dbpath/database/lock.wrm") {
    print LOG "*Reinitdb error - lock.wrm file present..\n";
    close LOG;
    die();
  }

  &delete_files_from("$dbpath/database",".\.wrm","-") or print LOG "ERROR: Problems removing files from $dbpath/database\n";
  
  my $command = "y\n";
  print LOG &runtime, ": reinitializing the database\n";
  &DbWrite($command,$tace,$dbpath,"ReInitDB");
 
  foreach my $filename (@filenames) {
    my $command = "pparse $filename\nsave\nquit\n";
    if (-e $filename) {
      my $runtime = &runtime;
      print LOG &runtime, ": parsing $filename\n";
      LOG->autoflush();
      my ($tsuser) = $filename =~ (/^\S+\/(\S+)\_/);
      &DbWrite($command,"$tace -tsuser $tsuser",$dbpath,"ParseFile");
      $runtime = &runtime;
      LOG->autoflush();
    }
    else {
      print LOG "* Reinitdb: $filename is not existent - skipping ..\n";
      next;
    }
  }
  print LOG &runtime, ": Finished.\n\n";

}


###################################################
# Alan Coulson maintains a ContigC                
# database in ~cemap for the physical map.
# This is dumped in file ~cemap/cen2hs.ace 
# (makecen2hs.pl nightly cron job on rathbin)

sub physical_map_stuff{

  print LOG &runtime, ": starting physical_map_stuff subroutine\n";

  
  my $command = "find clone\nedit -D pMap\nedit -D Fingerprint\nedit -D Contig9\nedit -D Remark\n";
  $command .= "pparse $autoacedir/physical_map/cen2hs.ace\nsave\nquit\n";
  &DbWrite($command,"$tace -tsuser Coulson",$dbpath,"ContigC");
  
  # now make the maps
  $command = "gif makemaps -all\nsave\ngif makemaps -seqclonemap $dbpath/acefiles/seqclonemap.ace\n";
  $command .= "pparse $dbpath/acefiles/seqclonemap.ace\nsave\nquit\n";
  &DbWrite($command,$giface,$dbpath,"MakeMaps");

  print LOG &runtime, ": finished\n\n";
}


###################################################
# Set the date correctly in displays.wrm          

sub setdate {
  my @t   = localtime ; while ($t[5] >= 100) { $t[5] -= 100 ; }
  my $dat = sprintf "%02d\/%02d\/%02d", $t[3], $t[4]+1, $t[5] ;
  my $status = move("$dbpath/wspec/displays.wrm", "$dbpath/wspec/displays.old");
  print "ERROR: Couldn't move file: $!\n" if ($status == 0);


  open(FILE,"$dbpath/wspec/displays.old") or do { print LOG "failed to open $dbpath/wspec/displays.old\n"; return 1;};
  open(NEWFILE,">$dbpath/wspec/displays.wrm") or do { print LOG "failed to open $dbpath/wspec/displays.wrm\n"; return 1;};
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
  unlink "$dbpath/wspec/displays.old" or print LOG "ERROR: Couldn't unlink file: $!\n";
}





#############################
# Make chromosomal links    #
#############################
#
# 010605 : rd  : removed extra "wormbase/" from paths following $wormbasedir
# 001218 : dl  : Altered path of file to fit wormbase designations
# 001307 : dl  : Added call to remove old file and state 'autoace' as the db to query

sub makechromlink {

  print LOG &runtime, ": starting makechromlink subroutine\n";
  unlink "$wormbasedir/misc/misc_chromlinks.ace" or print LOG "ERROR: Couldn't unlink file: $!\n";
  my $command = "$basedir/scripts/makeChromLinks.pl > $wormbasedir/misc/misc_chromlinks.ace";
  $command = "$basedir/scripts/makeChromLinks.pl -test > $wormbasedir/misc/misc_chromlinks.ace" if ($test);
  &run_command("$command"); 

  if (-z "$wormbasedir/misc/misc_chromlinks.ace") {
    print LOG "*Makechromlink: chromlinks.ace has ZERO size\n";  
    return;
  } 
  else {
    my $command = "pparse $wormbasedir/misc/misc_chromlinks.ace\nsave\nquit\n";
    &DbWrite($command,"$tace -tsuser make_autoace",$dbpath,"MakeChromLinks");
  }

  print LOG &runtime, ": Finished.\n\n";
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

  elsif ($error eq "Database") {
    # Normal help menu
    print "\n\nYou need to specify a name of a database to build (-buildautoace option) or a name\n";
    print "of an existing database where you can dump release files (-buildrelease option\n\n";
    exit (0);
  }
}


################################################
# Make the db files for distribution           *	
################################################

sub buildrelease{	

  print LOG &runtime, "Starting to build release files\n\n";

  # Remove old release files if present
  my $WS_previous   = $WS_current - 1;

  if (-e "$dbpath/release/database.WS"."$WS_previous".".4-0.tar.gz"){
    print LOG "Older WS version files exist, removing them\n";
    &delete_files_from("$dbpath/release","*WS"."$WS_previous"."*") or print LOG "ERROR: Problems removing files from $dbpath/release: $!\n";
  }
  
  my $dbname;
  
  open (DBWRM,"<$dbpath/wspec/database.wrm");
  while (<DBWRM>) {
    if ($_ =~ /^NAME\s+(\S+)/) {
      $dbname=$1;
    }
  }
  close DBWRM;
  print LOG "makedistr: dbname $dbname\n\n";
  
  &run_command("/bin/touch $dbpath/release/files_in_tar");
  &run_command("/bin/touch $dbpath/release/md5sum.${dbname}");
  
  my @tarfiles;
  $tarfiles[0] = "wspec/cachesize.wrm  wspec/constraints.wrm wspec/copyright wspec/database.wrm wspec/displays.wrm wspec/help.wrm wspec/layout.wrm wspec/models.wrm wspec/options.wrm wspec/passwd.wrm wspec/psfonts.wrm wspec/subclasses.wrm wspec/xfonts.wrm wgf wquery wscripts  pictures database/log.wrm database/database.map database/ACEDB.wrm" ;
  
  for (my $i = 1 ; -e "$dbpath/database/block$i.wrm" ; ++$i) {
    $tarfiles[($i+4)/5] .= " database/block$i.wrm" ;
  }
  print LOG "* Makedistr: beginning tar ..\n";
  &delete_files_from("$dbpath/release","database\.$dbname\..*\.tar","-") or print LOG "ERROR removing files from $dbpath/release\n";

  for (my $i = 0; $i < @tarfiles; ++$i) {
    &run_command("cd $dbpath; tar -hcf $dbpath/release/database.$dbname.4-$i.tar $tarfiles[$i]\"");
    
    # list files in the tar archive
    &run_command("tar -tf $dbpath/release/database.$dbname.4-$i.tar >> $dbpath/release/files_in_tar");
    
    # gzip the tar archive
    &run_command("/bin/gzip $dbpath/release/database.$dbname.4-$i.tar"); 
    
    # check consistency of gzip file
    &run_command("/bin/gzip -t $dbpath/release/database.$dbname.4-$i.tar.gz >> $dbpath/release/files_in_tar");
    
    # calculate md5sum for the gzip file
    &run_command("/nfs/disk100/wormpub/bin.ALPHA/md5sum $dbpath/release/database.$dbname.4-$i.tar.gz >> $dbpath/release/md5sum.$dbname");
  }
}





##################################################################################
#
# Simple routine which will run commands via system calls but also check the 
# return status of a system call and complain if non-zero, increments error check 
# count, and prints a log file error
#
##################################################################################

sub run_command{
  my $command = shift;
  print LOG &runtime, ": started running $command\n";
  my $status = system($command);
 if(($status >>8) != 0){
    $errors++;
    print LOG "ERROR: $command failed. \$\? = $status\n";
  }
  print LOG &runtime, ": finished running\n\n";

  # for optional further testing by calling subroutine
  return($status);
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
















