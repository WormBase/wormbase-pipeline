#!/usr/local/bin/perl5.8.0 -w
#
# make_acefiles.pl 
#
# dl1
#
# Generates the .acefiles from the primary databases as a prelim for building
# autoace.
#
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2006-11-01 14:26:18 $

#################################################################################
# Variables                                                                     #
#################################################################################

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use IPC::Open2;
use IO::Handle;
use Getopt::Long;
use File::Copy;
use File::Spec;
use File::Path;
use Log_files;
use Storable;

##############################
# command-line options       #
##############################

our $help;       # Help perdoc
our $debug;      # Debug mode, verbose output to runner only
my $test;        # If set, script will use TEST_BUILD directory under ~wormpub
my $db;
my $database;
my $basedir;
my $store;

GetOptions ("debug=s"    => \$debug,
	    		"help"       => \$help,
	   	 	"database=s" => \$database,
	    		"db:s"       => \$db,
	    		"test"       => \$test,
	    		"store:s"    => \$store
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


if( $database ) {
  my ($volume,$directories,$fname) = File::Spec->splitpath( $database );
  $basedir = $directories;
  chop $basedir; # remove the last /
}
##############################
# Script variables (run)     #
##############################

my $rundate     = $wormbase->rundate;
my $runtime     = $wormbase->runtime;
my $tace        = $wormbase->tace;
my $log         = Log_files->make_build_log( $wormbase );

$basedir     = $wormbase->basedir unless $basedir;
my $wormbasedir = $wormbase->autoace."/acefiles/primary";
mkpath( $wormbasedir ) unless -e ( $wormbasedir );
my $miscdir     = $wormbase->misc_static;

my $autoace_config = "$basedir/autoace_config/autoace.config";

# set WS version depending on whether in test mode
my  $WS_version = $wormbase->get_wormbase_version_name;

# help page
&usage("Help") if ($help);

# single database mode
if ($db) {
    &usage("Bad_database_name") unless (($db eq "camace") || ($db eq "stlace") || ($db eq "brigace") || ($db eq "cshace") || ($db eq "caltech") || ($db eq "geneace"));
}



# Main Loop
&mknewacefiles;


##############################
# tidy up                    #
##############################


close(STDOUT);
close(STDERR);

$log->mail;

exit (0);





#################################################################################
### Subroutines                                                               ###
#################################################################################



#################################################################################
# Erases old acefiles and make new ones                                         #
#################################################################################

sub mknewacefiles {
    
  my ($dbname,$dbdir,$targetdir,$exe,$object,$criteria,$follow);
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
      ($dbname)  = $1;
      $dbdir     = $wormbase->primary("$dbname");
      my $msg = "\n\n". $wormbase->runtime . " : Processing $dbname information in autoace_config\n";
      $log->write_to("$msg");

      # need to change dbpath if in test mode
      $targetdir = $wormbase->acefiles."/primaries/$dbname";
      mkpath( $targetdir ) unless ( -e "$targetdir" );
      $exe       = "$tace $dbdir";
      next;
    }
    
    # nuts and bolts of acefile generation
    # format:  database object criteria
    
    # clear out the @deletes array eachtime
    undef (@deletes);
    
    # parse filename
    ($filename) = (/^\S+\s+(\S+)/);
    $filepath   = "$targetdir/$filename";
    $extrafile  = "$targetdir/$filename.extra";

    print "Noting filename:$filename\n" if ($wormbase->debug);

    # if misc_static file copy from source dir under ~wormpub
    if ($dbname =~ /misc_static/) {
      my $misc_stat_dir = $wormbase->misc_static;
      system("cp -f $misc_stat_dir/$filename $filepath") and $log->write_to("ERROR : couldnt copy $misc_stat_dir/$filename : $!\n");
      next;
    }

    # deal with queries requiring tag deletion
    if (/\[(\S+.+)\]$/) { 
      @deletes = split (/\s/,$1);
      my $msg = $wormbase->runtime." : Adding " . scalar (@deletes) . " tags to delete array\n";
      $log->write_to("$msg");
      my $report = join ' + ', @deletes;

      if (/^\S+\s+\S+\s+(\S+)\s+\[/) {  
	$object = $1; 
	$criteria = ""; 
	
	$command = "nosave\nquery find $object\n";
	foreach my $delete (@deletes) {
	  $command .= "eedit -D $delete\n";
	}
	$command .= "show -a -T -f $filepath\nquit\n";
      } 
      elsif (/^\S+\s+\S+\s+(\S+)\s+(\S+.+)\s+\[/) {
	$object   = $1; 
	$criteria = $2; 
		
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
      $log->write_to("Adding $include tags\n");
      if ($follow) { $log->write_to("\t... and follow $follow\n");}
      
      if (/^\S+\s+\S+\s+(\S+)\s+\{/) {  
	  $object = $1; 
	  $criteria = ""; 
	  
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
	  $criteria = $2; 
	  
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
	$criteria = ""; 
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

    $log->write_to(": Dumping $object class\n");

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


