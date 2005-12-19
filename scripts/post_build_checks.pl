#!/usr/local/bin/perl5.8.0 -w
#
# wrapper to call introncheck, estcheck and overlapcheck.pl on gff files, 
# recommended after database rebuilt :o)
#
# by Kerstin Jekosch
# 13/07/01
#
# N.B. Previously called gffcheck
#
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2005-12-19 12:31:38 $

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use IO::Handle;
$|=1;


###########
# options #
###########

use vars qw / $opt_a $opt_o $opt_e $opt_i $opt_t $opt_h /;

$opt_a="";   # does all the following
$opt_o="";   # performs overlapcheck.pl
$opt_e="";   # performs estcheck 
$opt_i="";   # performs introncheck
$opt_t = ""; # performs TSL checks
$opt_h="";   # Help/Usage page


my ($help, $debug, $test, $verbose, $store, $wormbase);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store"      => \$store,
	    "a"          => \$opt_a,
	    "e"          => \$opt_e,
	    "i"          => \$opt_i,
	    "o"          => \$opt_o,
	    "t"          => \$opt_t,
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


#################
# set variables #
#################

my $rundate = $wormbase->rundate;
my $runtime = $wormbase->runtime;
my $WS_version = $wormbase->get_wormbase_version;


#########################
# tiny little main loop #
#########################

&runintroncheck  if ($opt_i);

&runestcheck     if ($opt_e);

&runoverlapcheck if ($opt_o);

&runTSLcheck     if ($opt_t);

##############################
# Tidy up                    #
##############################

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);




##############################################################
#
# Subroutines
#
##############################################################


##########################################


sub runestcheck {
  $log->write_to( $wormbase->runtime, ": Starting estcheck\n");
  $wormbase->run_script("estcheck", $log) && die "Can't run estcheck\n";
  $log->write_to( $wormbase->runtime, ": Finished running\n\n");
}

sub runTSLcheck {
  $log->write_to( $wormbase->runtime, ": Starting TSLcheck.pl\n");
  $wormbase->run_script("TSLcheck.pl", $log) && die "Can't run TSLcheck.pl\n";
  $log->write_to( $wormbase->runtime, ": Finished running\n\n");
}

sub runintroncheck {
  $log->write_to( $wormbase->runtime, ": Starting introncheck\n");
  $wormbase->run_script("introncheck", $log) && die "cant run introncheck\n";
  $log->write_to( $wormbase->runtime, ": Finished running\n\n");
}

sub runoverlapcheck {
  $log->write_to( $wormbase->runtime, ": Starting overlapcheck.pl\n");
  $wormbase->run_script("overlapcheck.pl", $log) && die "Can't run overlapcheck.pl\n";
  $log->write_to( $wormbase->runtime, ": Finished running\n\n");
}
   

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


__END__

=pod

=head2   NAME - post_build_checks.pl

=head1 USAGE

=over 4

=item post_build_checks.pl [-options]

=back

post_build_checks.pl is a wrapper to drive scripts to check the gff files for confirmed introns (introncheck), 
inconsistencies in EST assignments (estcheck), overlapping genes, ESTs matching introns and repeats 
within genes (overlapcheck.pl).  

post_build_checks.pl mandatory arguments:

=over 4

=item none, (but it won't do anything)

=back

post_build_checks.pl OPTIONAL arguments:

=over 4

=item -a, executes all of the following -eiowl

=item -e, executes estcheck

=item -t, executes TSLcheck

=item -i, executes introncheck 

=item -o, execute overlapcheck.pl

=back

=cut
