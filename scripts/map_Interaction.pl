#!/usr/local/bin/perl5.8.0 -w
#
# map_Interaction.pl
#
# Add information to Interaction objects via aceperl follows....
#
# by Dan Lawson
#
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-05-02 20:18:23 $

use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Ace;

#####################################
# variables and command-line options #
######################################

my $maintainers = "All";
my $output;     # output ace file
my $help;       # Help perldoc
my $test;       # Test mode
my $debug;      # Debug mode, verbose output to user running script
my $verbose;    # Verbose mode
my $load;       # for loading into autoace
my $store;      # specify a frozen configuration file

GetOptions(
    'debug=s'   => \$debug,
    'verbose'   => \$verbose,
    'test'      => \$test,
    'help'      => \$help,
    'load'      => \$load,
    'acefile=s' => \$output,
    'store=s'   => \$store
);

# Display help if required
&usage('Help') if ($help);

############################
# recreate configuration   #
############################
my $wb;
if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug, -test => $test, ) }

###########################################
# Variables Part II (depending on $wb)    #
###########################################
$test  = $wb->test  if $wb->test;     # Test mode
$debug = $wb->debug if $wb->debug;    # Debug mode, output only goes to one user
my $tace  = $wb->tace;                # tace executable path
my $dbdir = $wb->autoace;             # Database path

# Use debug mode?
if ($debug) {
    print "DEBUG = \"$debug\"\n\n";
    ( $maintainers = $debug . '\@sanger.ac.uk' );
}

$output = $output ? $output : "$dbdir/acefiles/Interaction_connections.ace";    # output file path
print "// Test mode:\n// searching against $dbdir\n// output written to $output\n\n" if $test;
my $log = Log_files->make_build_log($wb);

open( OUTPUT, ">$output" ) || die "Can't open output file $output\n";

#####################
# open a connection #
#####################

my $db = Ace->connect( -path    => "$dbdir", -program => $tace)
  || do { print "Connection failure: ", Ace->error; die(); };

my @interactions = $db->fetch( -class => 'Interaction', -name  => '*');

# Loop through each YH object
foreach my $interaction (@interactions) {

  print "// Interaction : \"$interaction\"\n" if ($verbose);
  
  if ( defined( $interaction->PCR_interactor ) ) {
    my $pcr_name = $interaction->PCR_interactor;
    
    my $pcr = $db->fetch( PCR_product => $pcr_name );
    
    if ( defined( $pcr->Overlaps_CDS ) ) {
      my @cds_names = $pcr->Overlaps_CDS;
      print OUTPUT "\nInteraction : \"$interaction\"\n";
      
      foreach my $cds_name (@cds_names) {
        print OUTPUT "Interactor_overlapping_CDS $cds_name\n";
        
        my $cds = $db->fetch( CDS => $cds_name );
        my $gene = $cds->Gene;
        print OUTPUT "Interactor_overlapping_Gene $gene\n";
        $cds->DESTROY();
      }
      print OUTPUT "\n";
      
      $pcr->DESTROY();
    }
  }
  
  if ( defined( $interaction->Sequence_interactor ) ) {
    my $target = $interaction->Sequence_interactor;

    my $seq = $db->fetch( Sequence => $target );

    if ( defined( $seq->Matching_CDS ) ) {
      my @cds_target = $seq->Matching_CDS;
      print OUTPUT "\nInteraction : \"$interaction\"\n";

      foreach my $cds_name (@cds_target) {
        print OUTPUT "Interactor_overlapping_CDS $cds_name\n";

        my $cds = $db->fetch( CDS => $cds_name );
        my $gene_target = $cds->Gene;
        print OUTPUT "Interactor_overlapping_Gene $gene_target\n";
        $cds->DESTROY();
      }
      print OUTPUT "\n";
      $seq->DESTROY();
    }
  }
  $interaction->DESTROY();
}

$db->close;

close(OUTPUT);    # close the output filehandle

###############
# hasta luego #
###############
if ($load) {
  $log->write_to("Loading file to autoace\n");
  $wb->load_to_database( $wb->autoace, $output, 'map_y2h_script', $log );
}

$log->mail();
exit(0);

###############################
# Prints help and disappears  #
###############################

##############################################################
#
# Subroutines
#
##############################################################

##########################################

sub usage {
    my $error = shift;

    if ( $error eq 'Help' ) {

        # Normal help menu
        system( 'perldoc', $0 );
        exit(0);
    }
}

############################################

__END__

=pod

=head2 NAME - map_Y2H.pl

=head1 USAGE

=over 4

=item map_Y2H.pl [-options]

=back

map_Y2H.pl connects Yeast2Hybruid objects to the Gene and CDS objects
via the PCR_bait and Sequence_target data stored within the database.

map_Y2H.pl mandatory arguments:

=over 4

=item none

=back

map_Y2H.pl optional arguments:

=over 4

=item -debug, Debug mode

=item -verbose, Verbose mode

=item -test, Test mode, generate the acefile but do not upload them 

=item -load, loads file to autoace

=item -help, Help pages

=item -store specify a configuration file

=item -outfile specify location for the output Ace (if you want to use it outside of the build)

=back

=cut
