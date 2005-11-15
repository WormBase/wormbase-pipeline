#!/usr/local/bin/perl5.8.0 -w
#
# map_Y2H.pl
#
# Add information to Y2H objects via aceperl follows....
#
# by Dan Lawson
#
# Last updated by: $Author: mh6 $                      
# Last updated on: $Date: 2005-11-15 14:55:53 $        

use strict;
use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Getopt::Long;
use Cwd;
use Ace;

#####################################
# variables and command-line options #
######################################

my $tace        = &tace;                                                       # tace executable path
my $dbdir       = "/wormsrv2/autoace";                                         # Database path
my $output      = "$dbdir/acefiles/Y2H_connections.ace";  # output file path
my $db_version  = &get_wormbase_version_name;                                  # WS version name

my %output      = (); # for Y2H

my $maintainers = "All";

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


# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = $debug . '\@sanger.ac.uk');
}

my $log = Log_files->make_build_log($debug);

if ($test) {
    $dbdir       = glob("~wormpub/TEST_BUILD/autoace");                    # Database path
    $output      = "$dbdir/acefiles/Y2H_connections.ace";               # output file path
    print "// Test mode:\n// searching against $dbdir\n// output written to $output\n\n";

}


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

open (OUTPUT, ">$output") || die "Can't open output file $output\n";

#####################
# open a connection # 
#####################

my $db = Ace->connect(-path  => "$dbdir",
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};
                                                                                                                       
my @Y2H = $db->fetch(-class => 'Y2H',
                               -name  => '*');

# Loop through each Y2H object

foreach my $Y2H (@Y2H) {

  print "// Y2H : \"$Y2H\"\n" if ($verbose);

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
   
}

$db->close;

close (OUTPUT);      # close the output filehandle

###############
# hasta luego #
###############
if($load){
  $log->write_to("Loading file to autoace\n");
  my $command = "autoace_minder.pl -load $dbdir/acefiles/Y2H_connections.ace -tsuser interpro_motifs";

  my $status = system($command);
  if(($status >>8) != 0){
    $log->write_to("ERROR: Loading Y2H_connections.ace file failed \$\? = $status\n");
  }
}

$log->mail("$maintainers","BUILD REPORT: $0");

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

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
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

=back

=cut
