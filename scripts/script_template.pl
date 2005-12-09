#!/usr/local/bin/perl5.8.0 -w
#
# script_template.pl                           
# 
# by Keith Bradnam                         
#
# Script to find candidate genes for splitting
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2005-12-09 16:50:52 $      

use strict;                                      
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
#use Ace;
#use Sequence_extract;
#use Coords_converter;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose);
my $maintainers = "All";

my $database;


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "database=s" => \$database,
	    );

my $log = Log_files->make_build_log($debug);

# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if ($debug) {
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

# in test mode?
if ($test) {
  print "In test mode\n";

}



##########################
# MAIN BODY OF SCRIPT
##########################

# main stuff goes here




# Close log files and exit
$log->write_to("\n\nStatistics\n");
$log->write_to("----------\n\n");

$log->write_to("put some statistics here");

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
# This should expand on your brief description above and add details of any options
# that can be used with the program.  Such documentation can be viewed using the perldoc
# command.


__END__

=pod

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

=back

This script does...blah blah blah

script_template.pl MANDATORY arguments:

=over 4

=item none

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Verbose/Debug mode
 
=back

=over 4

=item -test, Test mode, generate the acefile but do not upload themrun the script, but don't change anything

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back
                                                                                             

=head1 REQUIREMENTS

=over 4

=item This script must run on a machine which can see the /wormsrv2 disk.

=back

=head1 AUTHOR

=over 4

=item Keith Bradnam (krb@sanger.ac.uk)

=back

=cut
