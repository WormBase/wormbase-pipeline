#!/usr/local/bin/perl5.8.0 -w
#
# map_Y2H.pl
#
# Add information to Y2H objects via aceperl follows....
#
# by Dan Lawson
#
<<<<<<< map_Y2H.pl
# Last updated by: $Author: ar2 $                      
# Last updated on: $Date: 2005-12-16 11:18:55 $        
=======
# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2005-12-16 11:18:55 $
>>>>>>> 1.3.4.2

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Ace;

#####################################
# variables and command-line options #
######################################

<<<<<<< map_Y2H.pl
my $tace        = &tace;                                                       # tace executable path
my $dbdir       = "/wormsrv2/autoace";                                         # Database path
my $output      = "$dbdir/acefiles/Y2H_connections.ace";  # output file path
my $db_version  = &get_wormbase_version_name;                                  # WS version name

my %output      = (); # for Y2H

=======
>>>>>>> 1.3.4.2
my $maintainers = "All";
my $output;          # output ace file
my $help;            # Help perldoc
my $test;            # Test mode
my $debug;           # Debug mode, verbose output to user running script
my $verbose;         # Verbose mode
my $load;            # for loading into autoace
my $store;           # specify a frozen configuration file

GetOptions(
    'debug=s'   => \$debug,
    'verbose'   => \$verbose,
    'test'      => \$test,
    'help'      => \$help,
    'load'      => \$load,
    'acefile=s' => \$output,
    'store=s'   => \$store
);

<<<<<<< map_Y2H.pl
my $rundate = &rundate;
my $runtime = &runtime;

my $help;       # Help perldoc
my $test;       # Test mode
my $debug;      # Debug mode, verbose output to user running script
my $verbose;    # Verbose mode
my $load;       # for loading into autoace

GetOptions ("debug=s"   => \$debug,
	    "verbose"   => \$verbose,
 	    "test"      => \$test,
            "help"      => \$help,
	    "load"      => \$load);
=======
# Display help if required
&usage('Help') if ($help);
>>>>>>> 1.3.4.2

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
my $tace = $wb->tace;    # tace executable path
my $dbdir = $wb->autoace;               # Database path

# Use debug mode?
<<<<<<< map_Y2H.pl
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

my $log = Log_files->make_build_log($debug);

if ($test) {
    $dbdir       = glob("~wormpub/TEST_BUILD/autoace");                    # Database path
    $output      = "$dbdir/acefiles/Y2H_connections.ace";               # output file path
    print "// Test mode:\n// searching against $dbdir\n// output written to $output\n\n";

=======
if ($debug) {
    print "DEBUG = \"$debug\"\n\n";
    ( $maintainers = $debug . '\@sanger.ac.uk' );
>>>>>>> 1.3.4.2
}

$output = $output ? $output : "$dbdir/acefiles/Y2H_connections.ace";    # output file path
print "// Test mode:\n// searching against $dbdir\n// output written to $output\n\n" if $test;
my $log = Log_files->make_build_log($debug);

##############################
# Main                       #
##############################
my $bait;
my @CDS_bait;
my $PCR;
my $target;
my $seq;
my @CDS_target;
my $gene_bait;
my $gene_target;
my $CDS;
my $i;

#######################
# MAIN BODY OF SCRIPT #
#######################

open( OUTPUT, ">$output" ) || die "Can't open output file $output\n";

#####################
# open a connection #
#####################

my $db = Ace->connect(
    -path    => "$dbdir",
    -program => $tace
  )
  || do { print "Connection failure: ", Ace->error; die(); };

my @Y2H = $db->fetch(
    -class => 'Y2H',
    -name  => '*'
);

# Loop through each Y2H object

foreach my $Y2H (@Y2H) {

  print "// Y2H : \"$Y2H\"\n" if ($verbose);

<<<<<<< map_Y2H.pl
  if (defined ($Y2H->PCR_bait)) { 
    $bait = $Y2H->PCR_bait;
    
    $PCR = $db->fetch(PCR_product => $bait);
    
    if (defined ($PCR->Overlaps_CDS)) { 
      @CDS_bait =  $PCR->Overlaps_CDS;
      print OUTPUT "\nY2H : \"$Y2H\"\n";
      
      foreach  $i (@CDS_bait) {
	print OUTPUT "Bait_overlapping_CDS $i\n";
	
	$CDS = $db->fetch(CDS => $i);
	$gene_bait = $CDS->Gene;
	print OUTPUT "Bait_overlapping_Gene $gene_bait\n";
      }
      print OUTPUT "\n";
      
      $CDS->DESTROY();
      $PCR->DESTROY();
    }    
  }
  
  if (defined ($Y2H->Sequence_target)) { 
    $target = $Y2H->Sequence_target;
    
    $seq = $db->fetch(Sequence => $target);
    
    if (defined ($seq->Matching_CDS)) { 
      @CDS_target =  $seq->Matching_CDS;
      print OUTPUT "\nY2H : \"$Y2H\"\n";
      
      foreach $i (@CDS_target) {
	print OUTPUT "Target_overlapping_CDS $i\n";
	
	$CDS = $db->fetch(CDS => $i);
	$gene_target = $CDS->Gene;
	print OUTPUT "Target_overlapping_Gene $gene_target\n";
      }
      print OUTPUT "\n";
      
      $seq->DESTROY();
      
    }    
  }
  
  $Y2H->DESTROY();
   
=======
    if ( defined( $Y2H->PCR_bait ) ) {
        $bait = $Y2H->PCR_bait;

        $PCR = $db->fetch( PCR_product => $bait );

        if ( defined( $PCR->Overlaps_CDS ) ) {
            @CDS_bait = $PCR->Overlaps_CDS;
            print OUTPUT "\nY2H : \"$Y2H\"\n";

            foreach $i (@CDS_bait) {
                print OUTPUT "Bait_overlapping_CDS $i\n";

                $CDS = $db->fetch( CDS => $i );
                $gene_bait = $CDS->Gene;
                print OUTPUT "Bait_overlapping_Gene $gene_bait\n";
            }
            print OUTPUT "\n";

            $CDS->DESTROY();
            $PCR->DESTROY();
        }
    }

    if ( defined( $Y2H->Sequence_target ) ) {
        $target = $Y2H->Sequence_target;

        $seq = $db->fetch( Sequence => $target );

        if ( defined( $seq->Matching_CDS ) ) {
            @CDS_target = $seq->Matching_CDS;
            print OUTPUT "\nY2H : \"$Y2H\"\n";

            foreach $i (@CDS_target) {
                print OUTPUT "Target_overlapping_CDS $i\n";

                $CDS = $db->fetch( CDS => $i );
                $gene_target = $CDS->Gene;
                print OUTPUT "Target_overlapping_Gene $gene_target\n";
            }
            print OUTPUT "\n";

            $seq->DESTROY();

        }
    }

    $Y2H->DESTROY();

>>>>>>> 1.3.4.2
}

$db->close;

close(OUTPUT);    # close the output filehandle

###############
# hasta luego #
###############
<<<<<<< map_Y2H.pl
if($load){
  $log->write_to("Loading file to autoace\n");
  my $command = "autoace_minder.pl -load $dbdir/acefiles/Y2H_connections.ace -tsuser interpro_motifs";

  my $status = system($command);
  if(($status >>8) != 0){
    $log->write_to("ERROR: Loading Y2H_connections.ace file failed \$\? = $status\n");
  }
}
=======
if ($load) {
    $log->write_to("Loading file to autoace\n");
    my $command = "autoace_minder.pl -load $output -tsuser interpro_motifs";

    my $status = system($command);
    if ( ( $status >> 8 ) != 0 ) {
        $log->write_to("ERROR: Loading $output file failed \$\? = $status\n");
    }
}
>>>>>>> 1.3.4.2

<<<<<<< map_Y2H.pl
$log->mail("$maintainers","BUILD REPORT: $0");
=======
$log->mail( "$maintainers", "BUILD REPORT: $0" );
>>>>>>> 1.3.4.2

exit(0);

###############################
# Prints help and disappears  #
###############################

##############################################################
#
# Subroutines
#
##############################################################

<<<<<<< map_Y2H.pl

=======
>>>>>>> 1.3.4.2
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
