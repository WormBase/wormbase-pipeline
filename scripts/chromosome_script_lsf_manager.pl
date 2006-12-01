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


# Display help if required
&usage("Help") if ($help);

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
  my $mother = $m->submit("$command -chromosome $chrom $flags");
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

=head2 NAME - chromosome_script_lsf_manager.pl

=head1 USAGE

=over 4

=item chromosome_script_lsf_manager.pl [-options]

=back

This script runs a command once for each chromosome by appending '
-chromosome ' and then the name of the chromosome on a command-line
which is submitted at an LSF job to the normal queue.

chromosome_script_lsf_manager.pl MANDATORY arguments:

=over 4

=item -command, Command to run.

=back

This is the script that is run under LSF, once for each chromosome. It is expected that this script understands the argument '-chromosome' on its command-line.

chromosome_script_lsf_manager.pl  OPTIONAL arguments:

=over 4

=item -prefix

=back

If this is set then 'CHROMOSOME_' is prefixed to the name of the
chromosome when it is placed on the command-line of the script to be
run

=over 4

=item -mitochondrion

=back

If this is set then the mitochondrion is included in the list of chromsomes.

=over 4

=item -h, Help

=back

This help message.

=over 4
 
=item -debug, Debug mode
 
=back

Set this to the username who should receive the emailed log
messages. The default is that everyone in the group receives them.

=over 4

=item -test

=back

Test mode, run the script, but don't change anything.

=over 4
    
=item -verbose

=back

Output lots of chatty test messages.

=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Gary Williams (gw3@sanger.ac.uk)

=back

=cut
