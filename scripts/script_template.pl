#!/usr/local/bin/perl5.8.0 -w
#
# script_template.pl                           
# 
# by Keith Bradnam                         
#
# This is a example of a good script template
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2006-01-10 14:23:32 $      

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

my ($help, $debug, $test, $verbose, $store, $wormbase);


GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store"      => \$store,
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

#################################
# Set up some useful paths      #
#################################

# Set up top level base directories (these are different if in test mode)
my $basedir         = $wormbase->basedir;     # BASE DIR
my $ace_dir         = $wormbase->autoace;     # AUTOACE DATABASE DIR
my $wormpep_dir     = $wormbase->wormpep;     # CURRENT WORMPEP
my $wormrna_dir     = $wormbase->wormrna;     # CURRENT WORMRNA
my $common_data_dir = $wormbase->common_data; # AUTOACE COMMON_DATA
my $chromosomes_dir = $wormbase->chromosomes; # AUTOACE CHROMSOMES
my $reports_dir     = $wormbase->reports;     # AUTOACE REPORTS
my $gff_dir         = $wormbase->gff;         # AUTOACE GFF
my $gff_splits_dir  = $wormbase->gff_splits;  # AUTOACE GFF SPLIT
my $logs_dir        = $wormbase->logs;        # AUTOACE LOGS

# some database paths
my $geneace   = $wormbase->database('geneace');
my $camace    = $wormbase->database('camace');
my $currentdb = $wormbase->database('current');
my $stlace    = $wormbase->database('stlace');
my $citace    = $wormbase->database('citace');
my $cshace    = $wormbase->database('cshace');
my $brigace   = $wormbase->database('brigace');

# other paths
my $ftp_upload_dir  = $wormbase->ftp_upload;  # "/nfs/ftp_uploads/wormbase"
my $tace            = $wormbase->tace;        # TACE PATH
my $giface          = $wormbase->giface;      # GIFACE PATH



##########################
# MAIN BODY OF SCRIPT
##########################

# main stuff goes here

# example of running anther script
$wormbase->run_script("other_script -options", $log);


# Close log files and exit
$log->write_to("\n\nStatistics\n");
$log->write_to("----------\n\n");
$log->write_to("Put some statistics here.");

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

=item Keith Bradnam (krb@sanger.ac.uk)

=back

=cut
