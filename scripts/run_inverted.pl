#!/usr/local/bin/perl5.8.0
#
# run_inverted.pl
#
# dl 
#
# Usage : run_inverted.pk [-options] <sequence_file>

#################################################################################
# Initialise variables                                                          #
#################################################################################
 
use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Storable;
use Log_files;

##############################
# command-line options       #
##############################
                                                                                                                                     
my $debug;              # debug mode
my $test;
my $store;
my $help;               # Help/Usage page
my $sequence;           # Sequence file handle
my $all;                # Compute for all clones


GetOptions (
	    "all"         => \$all,
	    "sequence:s"  => \$sequence,
            "debug:s"     => \$debug,
            "help"        => \$help,
	    "test"        => \$test,
	    "store:s"     => \$store
	    );


# Help menu

if ($help) {
    exec ('perldoc',$0);
}

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}


# make log files
my $log = Log_files->make_build_log($wormbase);

# Crash out if no sequence filename supplied
unless ($all||$sequence) {
  $log->log_and_die("You must choose an option, [-all or -sequence <seqname>]. Ciao\n\n");
}

my %clone2sequence  = $wormbase->FetchData('clone2sequence');      # Genomic_canonical => DNA sequence
my %clonesize       = $wormbase->FetchData('clonesize');           # Genomic_canonical => length in bp


if (($sequence) && (!defined $clone2sequence{$sequence})) {
  $log->log_and_die("Sequence name '$sequence' is not recognised in the COMMON_DATA hash\n\n");
}

# deal with array of clones to process

my @clones2process = sort keys %clone2sequence;

if ($sequence) {
    @clones2process = "";
    push (@clones2process, $sequence);
}

# open output acefile

my $outputfilename = $wormbase->acefiles."inverted_repeats.ace";

open (ACE, ">$outputfilename") || die "Failed to open output acefile: '$outputfilename'\n";

# Loop through all clone which need to be dealt with

my ($score, $percent,$gaps);
my ($loop_1_start,$loop_1_stop);
my ($loop_2_start,$loop_2_stop);
my ($loop_len);
my $tag;
my @output;

foreach my $clone (@clones2process) {

    next if ($clone eq "");

    # generate sequence file for the query sequence    
    open (INPUT, ">/tmp/inverted_temp.$$") || die "Failed to open file for dumping sequence to\n";
    print INPUT ">$clone\n";
    print INPUT $clone2sequence{$clone};
    close INPUT;

    ################################
    # Run inverted on the sequence #
    ################################
    
    # inverted output format: 
    #
    # Score 52: 20/22 ( 90%) matches, 0 gaps
    #      3810 ataaaaactcgaattcaaaaaa 3831
    #                                 
    #      4164 tatttgtgagcttaattttttt 4143

    # Final ace output
    #
    # Feature Inverted        3810    3831    90      "loop 353"

    open (INV, "inverted /tmp/inverted_temp.$$ | ");
    while (<INV>) {
	
	chomp;
	next if ($_ eq "");     # Shortcut empty lines
	
	# Main data 
	if (/^Score (\d+)\: \S+ \( (\d+)\%\) matches\, (\d+) gaps/) {
	    ($score, $percent,$gaps) = ($1,$2,$3);
	    $tag = 1;
	    next;
	}
	elsif (/^Score (\d+)\: \S+ \(100\%\) matches\, (\d+) gaps/) {
	    ($score, $percent,$gaps) = ($1,"100",$3);
	    $tag = 1;
	    next;
	}
	
	# start loop
	if (($tag == 1) && (/^\s+(\d+) \S+ (\d+)/)) {
	    ($loop_1_start,$loop_1_stop) = ($1,$2);
	    $tag++;
	    next;
	} 
	
	# end loop
	if (($tag == 2) && (/^\s+(\d+) \S+ (\d+)/)) {
	    ($loop_2_start,$loop_2_stop) = ($1,$2);
	    
	    $loop_len = $loop_2_start - $loop_1_start - 1;
	    
	    # output ace format for both stem structures at the same time
	    
	    if ($gaps > 1) {
		push (@output, "Feature Inverted $loop_1_start $loop_1_stop $percent \"loop $loop_len, $gaps gaps\"\n");
		push (@output, "Feature Inverted $loop_2_start $loop_2_stop $percent \"loop $loop_len, $gaps gaps\"\n");
	    }
	    elsif ($gaps == 1) {
		push (@output, "Feature Inverted $loop_1_start $loop_1_stop $percent \"loop $loop_len, 1 gap\"\n");
		push (@output, "Feature Inverted $loop_2_start $loop_2_stop $percent \"loop $loop_len, 1 gap\"\n");
	    }
	    elsif ($gaps == 0) {
		push (@output, "Feature Inverted $loop_1_start $loop_1_stop $percent \"loop $loop_len\"\n");
		push (@output, "Feature Inverted $loop_2_start $loop_2_stop $percent \"loop $loop_len\"\n");
	    }
	    
	    $tag == 0;
	}
	
    }
    close INV;

    # tidy up and remove temp file
    system ("rm -f /tmp/inverted_temp.$$");

    if (scalar @output < 1) {
	next;
    }

    print ACE "Sequence : \"$clone\"\n";
    print ACE "Feature_data \"${clone}:inverted\" 1 $clonesize{$clone}\n\n";
    
    print ACE "Feature_data : \"${clone}:inverted\"\n";
    foreach (@output) {
	next if ($_ eq "");
	print ACE;
    }
    
    print ACE "\n";
    
    @output = "";

} #_ loop over each clone (@clones2process)


close ACE;
close STDOUT;

$log->mail();

# Hasta luego
exit (0);


__END__
                                                                                                                                     
=pod
                                                                                                                                     
=head2   NAME - run_inverted.pl

=head1 USAGE
                                                                                                                                     
=over 4
                                                                                                                                     
=item run_inverted.pl [-options]
                                                                                                                                     
=back
                                                                                                                                     
run_inverted.pl is a wrapper to run Richard Durbin's inverted script over the
genome sequences. You can run the whole data set or specify a single sequence
name. 

All data is taken from the COMMON_DATA hashes hosted on wormsrv2.

Output files are written to /wormsrv2/autoace/acefiles

=back

=head1 run_inverted.pl MANDATORY arguments (one of the two):

=over 4 

=item -all, Processes all genome sequences

=back 

=item -sequence <sequence_name>, Process a single genome sequence

=back

=head1 run_inverted.pl OPTIONAL arguements:

=item -help, these help pages

=back
 
=cut
