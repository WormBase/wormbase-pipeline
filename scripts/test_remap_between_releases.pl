#!/usr/local/bin/perl5.8.0 -w 
#
#   test_remap__between_releases.pl                 
# 
# by Gary Williams                         
#
# This tests to see if we need to remap
# sequence positions between two releases.
#
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2007-05-15 15:30:26 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

use Modules::Remap_Sequence_Change;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($release1, $release2, $version);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "release1=i" => \$release1,
	    "release2=i" => \$release2,
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


if (! defined $release1 || ! defined $release2) {
  die "Specify the release numbers to use\n";
}


##########################
# read in the mapping data
##########################

my @mapping_data = Remap_Sequence_Change::read_mapping_data($release1, $release2);

if (Remap_Sequence_Change::remap_test($release1, $release2, @mapping_data)) {
  # there are changes
  $log->write_to("WARNING: There have been genomic sequence changes.\nThe remapping programs will therefore be run.\nThis may take some time.\n");

  $log->mail();
  exit(1);

} else {
  # there are no changes
  $log->write_to("There are no genomic sequence changes.\nThe remapping programs will not be run.\n");

  $log->mail();
  exit(0);
}





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

=head2 NAME - test_remap__between_releases.pl

=head1 USAGE

=over 4

=item test_remap__between_releases.pl [options]

=back

This script tests to see if there have been genomic sequemce changes
bewteen two releases and exits with status 1 if there have been
changes. If there are no changes it exits with status 0.

script_template.pl MANDATORY arguments:

=over 4

=item -release1 The first (earlier) database to convert from e.g. 140

=back

=item -release2 The second (later) database to convert to e.g. 155

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

=item Gary Williams

=back

=cut
