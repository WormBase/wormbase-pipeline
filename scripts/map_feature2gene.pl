#!/usr/local/bin/perl
#
# map_feature2gene.pl
#
# make the connections between TSL features and CDS/Gene
#
# Dan Lawson
#
# Usage : map_feature2gene.pl [-options]
#
# Last edited by: $Author: dl1 $
# Last edited on: $Date: 2005-04-26 13:45:15 $
 
#################################################################################
# Initialise variables                                                          #
#################################################################################

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Getopt::Long;
use vars;
use Ace;
 
#####################################
# variables and command-line options #
######################################
                                                                                                 
my $tace        = &tace;                                                       # tace executable path
my $dbdir       = "/wormsrv2/autoace";                                         # Database path
my $output      = "$dbdir/acefiles/TSL_feature_connections.ace";               # output file path
my $db_version  = &get_wormbase_version_name;                      

my $test;               # In test build mode
my $help;               # Help/Usage page
my $verbose;            # Verbose mode
my $load;               # upload acefile after run

GetOptions (
            "help"           => \$help,
	    "load"           => \$load,
	    "verbose"        => \$verbose,
            "test"           => \$test,
);

# Help pod if needed
&usage(0) if ($help);

#######################
# MAIN BODY OF SCRIPT #
#######################

open (OUTPUT, ">$output") || die "Can't open output file $output\n";

#####################
# open a connection #
#####################

print "// Connect to database $dbdir\n" if ($verbose);
my $db = Ace->connect(-path  => "$dbdir",
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};
print "// done\n" if ($verbose);

my @TSL = $db->fetch(-class => 'Feature',
                               -name  => '*');

# Loop through each Feature object
  
my $method;
my @sequence_name;
my $sequence;
my @CDS;
my $seq;
my $i;
my $cds;
my %match;
my $sum;

foreach my $feature (@TSL) {
    
    $method = $feature->Method(1);
    
    next unless ( ($method eq "SL1") || ($method eq "SL2") );
    
    @sequence_name = $feature->Defined_by_sequence(1);

    undef %match ;           # reset hash
    
    foreach $i (@sequence_name) {
	print "// defined by sequence $i\n" if ($verbose);

	$sequence = $db->fetch(Sequence => $i);
	@CDS = $sequence->Matching_CDS(1);
	
	foreach $cds (@CDS) {
	    next if ($cds eq "");
	    $match{$cds} = 1;
	    print "// matches CDS '$cds'\n" if ($verbose);
	}
    }
    
    # write output and error tracking
    unless (defined %match) {
	print "\n// No connection found for Feature : $feature\n" ;
    }
    else {
	print OUTPUT "Feature : \"$feature\"\n";
	foreach my $j (keys %match) {
	    print OUTPUT "Associated_with_CDS $j\n";
	}
    }
    print OUTPUT "\n";

    $sequence->DESTROY();
    $feature->DESTROY();
    
}


close OUTPUT;

# Upload file to autoace (if you have been asked to)

if ($load) {
    
    my $command = "autoace_minder.pl -load $output -tsuser TSL_CDS_connect";
    
    my $status = system($command);
    if ( ($status >> 8) != 0 ) {
	print "ERROR: Loading $output file failed \$\? = $status\n";
    }
}

exit(0);



