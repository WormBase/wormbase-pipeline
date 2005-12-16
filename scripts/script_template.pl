#!/usr/local/bin/perl5.8.0 -w
#
# script_template.pl                           
# 
# by Keith Bradnam                         
#
<<<<<<< script_template.pl
# Script to find candidate genes for splitting
=======
# This is a example of a good script template
>>>>>>> 1.8.6.3
#
# Last updated by: $Author: ar2 $     
# Last updated on: $Date: 2005-12-16 11:18:55 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
<<<<<<< script_template.pl
use Log_files;
#use Ace;
#use Sequence_extract;
#use Coords_converter;

=======
use Log_files;
use Storable;
#use Ace;
#use Sequence_extract;
#use Coords_converter;

>>>>>>> 1.8.6.3

######################################
# variables and command-line options # 
######################################

<<<<<<< script_template.pl
my ($help, $debug, $test, $verbose);
my $maintainers = "All";

my $database;

=======
my ($help, $debug, $test, $verbose, $store, $wormbase);

my $database;

>>>>>>> 1.8.6.3

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
<<<<<<< script_template.pl
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "database=s" => \$database,
	    );

my $log = Log_files->make_build_log($debug);
=======
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "database=s" => \$database,
	    "store"      => \$store,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}
>>>>>>> 1.8.6.3

# Display help if required
&usage("Help") if ($help);

<<<<<<< script_template.pl
# Use debug mode?
if ($debug) {
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
=======
# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);
# $debug = `whoami` unless ($debug);


>>>>>>> 1.8.6.3
}

<<<<<<< script_template.pl
# in test mode?
if ($test) {
  print "In test mode\n";

}
=======
# establish log file.
my $log = Log_files->make_build_log($wormbase);

#################################
# Set up some useful paths      #
#################################

# Set up top level base directories (these are different if in test mode)
my $ace_dir    = $wormbase->autoace; # AUTOACE DATABASE DIR
my $ftp_upload = $wormbase->ftp_upload;	# "/nfs/ftp_uploads/wormbase"
my $wormpep    = $wormbase->wormpep; # CURRENT WORMPEP
my $wormrna    = $wormbase->wormrna; # CURRENT WORMRNA
my $data_dir   = $wormbase->common_data; # AUTOACE COMMON_DATA
my $chromosomes = $wormbase->chromosomes; # AUTOACE CHROMSOMES
my $reports    = $wormbase->reports; # AUTOACE REPORTS
my $gff        = $wormbase->gff;# AUTOACE GFF
my $gff_splits = $wormbase->gff_splits; # AUTOACE GFF SPLIT FILES
my $tace       = $wormbase->tace; # TACE PATH
my $giface     = $wormbase->giface; # GIFACE PATH



>>>>>>> 1.8.6.3



##########################
# MAIN BODY OF SCRIPT
##########################

<<<<<<< script_template.pl
# main stuff goes here

=======
# main stuff goes here

# example of running anther script
$wormbase->run_script("other_script -options", $log);
>>>>>>> 1.8.6.3



# Close log files and exit
$log->write_to("\n\nStatistics\n");
$log->write_to("----------\n\n");

$log->write_to("put some statistics here");

<<<<<<< script_template.pl
$log->mail();
=======
# Close log files and exit
$log->write_to("\n\nStatistics\n");
$log->write_to("----------\n\n");
>>>>>>> 1.8.6.3

<<<<<<< script_template.pl
print "Finished.\n" if ($verbose);
=======
$log->write_to("Put some statistics here\n");

$log->mail();
print "Finished.\n" if ($verbose);
>>>>>>> 1.8.6.3
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

=head2 NAME - script_template.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

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

<<<<<<< script_template.pl
=over 4
 
=item -debug, Verbose/Debug mode
 
=back

=over 4

=item -test, Test mode, generate the acefile but do not upload themrun the script, but don't change anything

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back
                                                                                             

=======
=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back



>>>>>>> 1.8.6.3
=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Keith Bradnam (krb@sanger.ac.uk)

=back

=cut
