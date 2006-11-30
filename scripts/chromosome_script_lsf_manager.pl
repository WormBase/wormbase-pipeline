#!/nfs/team71/worm/mh6/bin/perl
#####################################################
# Put a job on the queue for each chromosome and wait for them to exit.
# The command to run on each chromosome is expected to understand the
# command-line option -chrom <chromosome> e.g. -chrom IV
####################################################
use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;

##############################
# Script variables (run)     #
##############################

my ( $help, $debug, $test, $store );
my $command;			# the command to run for each chromosome
my $prefix;			# if set then chromosome has a 'CHROMSOME_' prefix
my $mito;			# if set then the mitochondrion is included in the set of chromosomes
my $verbose;    # for toggling extra output

##############################
# command-line options       #
##############################

GetOptions(
    "help"      => \$help,
    "debug=s"   => \$debug,
    "test"      => \$test,
    "store:s"   => \$store,
    "command:s" => \$command,      # command to execute for each chromosome
    "prefix"    => \$prefix,       # if set then chromosome has a 'CHROMSOME_' prefix
    "mitochondrion"    => \$mito,  # if set then the mitochondrion is included in the set of chromosomes
);

# recreate configuration ##########
my $wormbase;
my $flags = "";
if ($store) {
    $wormbase = Storable::retrieve($store) or croak("cant restore wormbase from $store\n");
    $flags = "-store $store";
}
else {
    $wormbase = Wormbase->new( -debug => $debug, -test => $debug );
    $flags .= "-debug $debug " if $debug;
    $flags .= "-test "         if $test;
}

# Variables Part II (depending on $wormbase)
$debug = $wormbase->debug if $wormbase->debug;    # Debug mode, output only goes to one user

my $log = Log_files->make_build_log($wormbase);

####################################

my $m      = LSF::JobManager->new();

my @chroms = $wormbase->get_chromosome_names(-mito =>$mito, -prefix => $prefix);

foreach my $chrom ( @chroms ) {
  my $mother = $m->submit("$command -chrom $chrom $flags");
}

$m->wait_all_children( history => 1 );
print "All children have completed!\n";

############################
for my $job ( $m->jobs ) {    # much quicker if history is pre-cached
    $log->write_to("$job exited non zero\n") if $job->history->exit_status != 0;
}
$m->clear;                    # clear out the job manager to reuse.

#################

$log->mail();
