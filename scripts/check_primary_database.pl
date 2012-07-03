#!/software/bin/perl -w
#
# check_primary_database.pl
# 
# by pad                         
#
# This script checks all of the species primary sequence databases....and geneace 
# so that gene discrepancies can be identified early in the build.
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2012-07-03 11:40:01 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
#use Ace;
#use Sequence_extract;
#use Coords_converter;

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $species,);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "organism=s" => \$species, # Specify an organism
	    "store:s"    => \$store,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species,
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

#################################
# Set up some useful paths      #
#################################

my $logs_dir        = $wormbase->logs;        # AUTOACE LOGS
my $xaceinstances;
my @xaceinstances;
# some database paths
my $geneace   = $wormbase->primary('geneace');
print "$geneace\n" if ($debug);
my $camace    = $wormbase->primary('camace');
print "$camace\n" if ($debug);
my $stlace    = $wormbase->primary('stlace');
print "$stlace\n" if ($debug);
my $citace    = $wormbase->primary('citace');
print "$citace\n" if ($debug);
my $cshace    = $wormbase->primary('cshace');
print "$cshace\n" if ($debug);
my $brigace   = $wormbase->primary('brigace');
print "$brigace\n" if ($debug);
my $remace    = $wormbase->primary('remace');
print "$remace\n" if ($debug);
my $brenace   = $wormbase->primary('brenace');
print "$brenace\n" if ($debug);
my $japace    = $wormbase->primary('japace');
print "$japace\n" if ($debug);

##########################
# MAIN BODY OF SCRIPT
##########################

# Load the array with databases to check based on the species.
if ((-e $camace) && ($species eq 'elegans'))  {push(@xaceinstances,"$camace");}
if ((-e $stlace) && ($species eq 'elegans')) {push(@xaceinstances,"$stlace");}
if ((-e $citace) && ($species eq 'elegans')) {push(@xaceinstances,"$citace");}
if ((-e $cshace) && ($species eq 'elegans')) {push(@xaceinstances,"$cshace");}
if ((-e $brigace) && ($species eq 'briggsae')) {push(@xaceinstances,"$brigace");}
if ((-e $remace) && ($species eq 'remanei')) {push(@xaceinstances,"$remace");}
if ((-e $brenace) && ($species eq 'brenneri')) {push(@xaceinstances,"$brenace");}
if ((-e $japace) && ($species eq 'japonica')) {push(@xaceinstances,"$japace");}
# example of running anther script

foreach $xaceinstances (@xaceinstances) {
  print "running camace_nameDB_comm.pl -database $xaceinstances\n";
  $wormbase->run_script("NAMEDB/camace_nameDB_comm.pl -database $xaceinstances", $log);
  
  # Run additional checks for camace and stlace
  if (($xaceinstances eq "/nfs/wormpub/BUILD/PRIMARIES/camace") || ($xaceinstances eq "/nfs/wormpub/BUILD/PRIMARIES/stlace")){
    print "Checking $xaceinstances for gene curation errors\n\nrunning check_predicted_genes.pl -database $xaceinstances\n";
    $wormbase->run_script("check_predicted_genes.pl -database $xaceinstances -basic", $log);
  }
}

# Checking primary gene database
if ((-e $geneace) && ($species eq 'elegans'))  {
$log->write_to("running geneace_nameDB_comm.pl\n");
  $wormbase->run_script("NAMEDB/geneace_nameDB_comm.pl", $log);
}

# Close log files and exit
$log->write_to("\n\nFinished\n");
$log->write_to("----------\n\n");
#$log->write_to("Put some statistics here.\n");

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






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




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - check_primary_database.pl

=head1 USAGE

=over 4

=item check_primary_database.pl  [-options]

=back

This script does...blah blah blah

script_template.pl MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item pad (pad@sanger.ac.uk)

=back

=cut
