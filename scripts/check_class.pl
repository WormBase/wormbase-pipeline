#!/usr/local/bin/perl5.8.0 -w
#
# dbcomp.pl
#
# Counts the number of objects in an ACEDB database for each Class stated in the config file
# Compares this number to those from a second database.
#
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2007-02-28 14:30:01 $


use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use IO::Handle;


$|=1;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($database, $database1, $database2, $classes);
our ($db_1, $db_2, $dbname_1, $dbname_2);

GetOptions (
	    "help"          => \$help,
            "debug=s"       => \$debug,
	    "test"          => \$test,
	    "verbose"       => \$verbose,
	    "store:s"         => \$store,
	    "database=s"    => \$database,
	    "database1=s"    => \$database1,
	    "database2=s"   => \$database2,
	    "classes=s"        => \$classes,
	    );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


############################################
# Get paths and versions
############################################

my $WS_current  = $wormbase->get_wormbase_version;
my $WS_previous = $WS_current - 1;
my $exec        = $wormbase->tace;

#################################################################
# Compare autoace to previous build unless -database specified
#################################################################

$dbname_1    = "WS${WS_previous}";
$db_1        = $wormbase->database("WS${WS_previous}");

$dbname_2    = "WS${WS_current}";
$db_2        = $wormbase->autoace;

$database = $database1 unless $database;

# First alternative database specified?
if ($database) {
    $dbname_1  = "$database";
    $db_1      = "$database"; 
}

# Second alternative database specified?
if ($database2) {
    $dbname_2  = "$database2";
    $db_2      = "$database2"; 
}

#########################################################################
# Main part of script
#########################################################################

my @classes = split /[\s,;]/, $classes;

$log->write_to("Checking $dbname_1 vs $dbname_2 for classes:\n@classes\n\n");

my ($class_count_1, $class_count_2) = &count_classes($db_1, $db_2, @classes);

$log->write_to(sprintf("%-12s\t%-7s\t%-7s\t%-7s\n", "CLASS",$dbname_1,$dbname_2,"Difference"));

my $count = 0;
foreach my $class (@classes) {

  ##################################################
  # Calculate difference between databases         #
  ##################################################
  
  my $count1 = $$class_count_1[$count];
  my $count2 = $$class_count_2[$count];
  my $diff = $count2 - $count1;
  my $err = "";
  if ($count2 == 0 || 
      $count2 < $count1 * 0.9 || 
      $count2 > $count1 * 1.1) {
    $err = "***** POSSIBLE ERROR *****";
    $log->error;
  }

  $log->write_to(sprintf("%-12s\t%7d\t%7d\t%7d\t%s\n", $class,$count1,$count2,$diff,$err));
  $count++;
}


# Email log file
$log->mail();

print "Finished.\n" if ($verbose);
exit (0);




##############################################################
#
# Subroutines
#
##############################################################



##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################

sub count_classes {

  my ($db_1, $db_2, @classes)  = @_;

  my @class_count1;
  my @class_count2;

  # Formulate query
  my $command;
  foreach my $class (@classes) {
    $command .= "query find '$class'\n";
  }
  $command .= "quit\n";


  ####################################
  # Count objects in first database
  ####################################
    
  # open tace connection and count how many objects in each class
  open (TACE, "echo '$command' | $exec $db_1 | ");
  while (<TACE>) {
    (push @class_count1, $1) if (/^\/\/ (\d+) Active Objects/);
  }
  close (TACE);
    
    
  ####################################
  # Count objects in second database
  ####################################
    
  # open tace connection and count how many objects in each class
  open (TACE, "echo '$command' | $exec $db_2 | ");
  while (<TACE>) {
    (push @class_count2, $1) if (/^\/\/ (\d+) Active Objects/);
  }
  close (TACE);

  return (\@class_count1, \@class_count2);
    
}


###############################################

=pod

=head2   NAME - check_class.pl

=head1 USAGE

=over 4

=item check_class.pl [-options]

=back

This script will (by default) perform a class by class comparison of autoace and 
the previous WormBase release database.  It will calculate what database this
will be but you can force the choice of the second database by using the -database
option.  Classes to be compared are specified within the script and should be amended
as necessary when classes are added/removed to the database.


=over 4

=item MANDATORY arguments: none

=back

=over 4

=item OPTIONAL arguments: -debug, -help, -database


-debug and -help are standard Wormbase script options.

-classes are a comma-delimited list of classes to check

-database allows you to specify the path of a database which will
then be compared to the BUILD/autoace directory.

-database2 allows you to specify a second database to compare to that specified
by -database (rather than compare with previous build)


=back

=head1 AUTHOR - Gary Williams



=cut

