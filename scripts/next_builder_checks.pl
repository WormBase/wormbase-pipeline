#!/usr/local/bin/perl5.8.0 -w
# next_builder_checks.pl
#
# by Keith Bradnam
#
# A simple script to send a check list to the person who will be performing the next
# build to check the current build
#
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2005-12-19 14:12:39 $
use strict;
use warnings;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;

my $store;                                          # to specify saved commandline arguments
my $maintainers = "All";
my ($help,$debug,$user);

GetOptions(
    "help"    => \$help,
    "debug=s" => \$debug,
    "user=s"  => \$user,
    'store=s' => \$store
);

# Display help if required
&usage("Help") if ($help);

############################
# recreate configuration   #
############################
my $wb;
if ($store) { $wb = Storable::retrieve($store) or croak("cant restore wormbase from $store\n") }
else { $wb = Wormbase->new( -debug => $debug, -test => $debug, ) }

##########################################
# Variables Part II (depending on $wb)    #
###########################################
$debug = $wb->debug if $wb->debug;    # Debug mode, output only goes to one user
my $WS_current = $wb->get_wormbase_version;


# Use debug mode?
if ($debug) {
    print "DEBUG = \"$debug\"\n\n";
    $maintainers = $debug.'\@sanger.ac.uk';
}

# Use -user setting to send email
if ($user) {
    $maintainers = $user.'\@sanger.ac.uk';
}

my $log = Log_files->make_build_log($wb);

#######################
#
# Main part of script
#
#######################

$log->write_to("\nGreetings $user.  You are meant to be building the next release of WormBase, so\n");
$log->write_to("you should take some time to check the current release.\n");
$log->write_to("====================================================================================\n\n");
$log->write_to("THINGS TO CHECK:\n\n");

my $log_msg= <<'LOGMSG';
1) The following 12 clones are representative of the whole genome in that
they include one Sanger and one St. Louis clone for each chromosome.  Check
each clone to ensure that it contains BLAT data (EST and mRNA), BLAST data,
waba data, gene models, UTRs etc.  Also check for presence of tandem and inverted
repeats which have gone missing in the past

i)    C25A1
ii)   F56A3
iii)  C04H5
iv)   B0432
v)    C07A9
vi)   F30H5
vii)  C10C6
viii) B0545
ix)   C12D8
x)    K04F1
xi)   C02C6
xii)  AH9

2) If the genome sequence has changed, you should inspect the clones containing
those changes to see if there are any strange errors (e.g. duplicate sets of data
which are slightly out of sync.)

3) Check /wormsrv2/autoace/CHROMOSOMES/composition.all - are there any non-ATCGN
characters

4a) Check that the latest WormPep proteins have proper protein and motif homologies
This has been a problem in some builds where all new WormPep proteins don't get any
BLAST analysis.  Pick a few random Wormpep proteins and especially check that all of
the various blastp homologies are there (human, fly, worm, yeast etc.) and try to
check at least one protein from the /wormsrv2/WORMPEP/wormpepXXX/new_entries.WSXXX file

4b) Now that we have a curated set of brigpep, should do this periodically for
C. briggase protein objects too...these now have their own set of blastp hits

5) Open a keyset of all genome sequences and then query them for:
 'Transcript AND NEXT AND NOT NEXT'
I.e. tRNAs not attached properly to parent sequence.  This has happened before and you
should notify the responsible group (HX or RW) to fix them for next build

6) Check all Protein objects have a Species tag set

7) Check all PFAM Motif objects have a title tag

8) Run composition *.dna in the CHROMOSOMES directory.  Make sure this is the same as it was
at the start of the build.  Bad Homol objects can lead to errors esp when chromosome length has been reduced

9) Confirm that GFFmunger.pl has done its job ie 
 grep landmark /wormsrv2/autoace/CHROMOSOMES/CHROMOSOME_III.gff | head
 grep five_prime /wormsrv2/autoace/CHROMOSOMES/CHROMOSOME_III.gff
 grep Partially /wormsrv2/autoace/CHROMOSOMES/CHROMOSOME_III.gff  and check these are overloaded CDS lines.
 These should have about 10 columns if overloaded and only about 5 or 6 if there are problems.

That's all...for now!  If you are satisfied the build is ok, please inform the person
building the database. Please continue to add to this list as appropriate

========================================================================================
LOGMSG

$log->write_to($log_msg);

$log->mail( "$maintainers", "Please check the ongoing build of WS${WS_current}");

exit(0);

##############################################################
#
# Subroutines
#
##############################################################

sub usage {
    my $error = shift;

    if ( $error eq "Help" ) {

        # Normal help menu
        system( 'perldoc', $0 );
        exit(0);
    }
}

##################################################################

__END__

=pod

=head1 NAME - next_builder_checks.pl

=head1 USAGE

=over 4

=item next_builder_checks.pl --user <user>

=back

This script simply sends a list of check items to the next person doing the build.
They should 'sign off' on each build and hand back to the main person when they
are happy all is ok.

=item MANDATORY arguments: -user <valid unix username>

Needed to send email

=back

=over 4

=item OPTIONAL arguments: -help, -debug <user>


 
=head1 AUTHOR - Keith Bradnam

Email krb@sanger.ac.uk


=cut
