#!/usr/local/bin/perl5.8.0 -w
#
# make_keysets.pl
#
# dl 021025
#
# Usage : make_keysets.pl [-options]
#
# Last edited by: $Author: dl1 $
# Last edited on: $Date: 2004-10-08 15:58:40 $

#################################################################################
# variables                                                                     #
#################################################################################

use IO::Handle;
use Getopt::Long;
use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Socket;


##############################
# command-line options       #
##############################


use vars qw/ $all $rnai $expr $history $pcr $go $touched $noload $debug $help $test/;
GetOptions (
	    "all"       => \$all,
	    "rnai"      => \$rnai,
	    "expr"      => \$expr,
	    "pcr"       => \$pcr,
	    "go"        => \$go,
	    "touched"   => \$touched,
	    "history=s" => \$history,
	    "noload"    => \$noload,
	    "test"      => \$test,
	    "debug"     => \$debug,
	    "help"      => \$help
	    );

&usage('help') if ($help);

# status report to stdout?

##############################
# Paths etc                  #
##############################

our $tace     = &tace;                                    # tace executable path
our $dbpath   = "/wormsrv2/autoace";                      # Database path (live runs)
our $outpath  = "/wormsrv2/tmp";                          # Output directory

# Use current_DB if you are testing 
$dbpath   = "/nfs/disk100/wormpub/DATABASES/current_DB" if ($test);        # Database path (test runs)

# only tell Dan if running debug mode
my $maintainers = "dl1\@sanger.ac.uk" if ($debug);

#################################################
# series of if loops for each keyset to be made #
#################################################

# CDS with RNAi
if (($rnai) || ($all)) {
    print "Tace query for CDS_with_RNAi\t" if ($debug);
    my $command ="nosave\nquery find elegans_CDS where RNAi_result\nlist -a\nquit\n";
    &tace_it($command,'CDS_with_RNAi');
    print "....load into db\t" if ($debug);
    &load_it('CDS_with_RNAi','keyset_lists') unless ($noload);
    print "...hasta luego\n\n" if ($debug);
}

# CDS with Expr_pattern
if (($expr) || ($all)) {
    print "Tace query for CDS_with_Expr_pattern\t" if ($debug);
    my $command ="nosave\nquery find elegans_CDS where Expr_pattern\nlist -a\nquit\n";
    &tace_it($command,'CDS_with_Expr_pattern');
    print "....load into db\t" if ($debug);
    &load_it('CDS_with_Expr_pattern','keyset_lists') unless ($noload);
    print "...hasta luego\n\n" if ($debug);
}

# CDS with PCR_product
if (($pcr) || ($all)) {
    print "Tace query for CDS_with_PCR_product\t" if ($debug);
    my $command ="nosave\nquery find elegans_CDS where Corresponding_PCR_product\nlist -a\nquit\n";
    &tace_it($command,'CDS_with_PCR_product');
    print "....load into db\t" if ($debug);
    &load_it('CDS_with_PCR_product','keyset_lists') unless ($noload);
    print "...hasta luego\n\n" if ($debug);
}

# CDS with GO_term
if (($go) || ($all)) {
    print "Tace query for CDS_with_GO_term\t" if ($debug);
    my $command ="nosave\nquery find elegans_CDS where GO_term\nlist -a\nquit\n";
    &tace_it($command,'CDS_with_GO_term');
    print "....load into db\t" if ( ($debug) && (!$noload) );
    &load_it('CDS_with_GO_term','keyset_lists') unless ($noload);
    print "...hasta luego\n\n" if ($debug);
}

# CDS with transcript_evidence
if (($touched) || ($all)) {
    print "Tace query for CDS_with_transcript_evidence\t" if ($debug);
    my $command ="nosave\nquery find elegans_CDS where Matching_cDNA\nlist -a\nquit\n";
    &tace_it($command,'CDS_with_transcript_evidence');
    print "....load into db\t" if (($debug) && (!$noload));
    &load_it('CDS_with_transcript_evidence','keyset_lists') unless ($noload);
    print "...hasta luego\n\n" if ($debug);
}

# wormpep histories
if (($history) || ($all)) {
    print "Tace query for wormpep_mods_since_WSnn\n" if ($debug);

    my @releases = (77,100,110,120,130);
    push (@releases,$history-1);

    foreach my $i (@releases) {
	my $wsname = "wormpep_mods_since_WS" . $i;
	(print "Tace query for " . $wsname . "\t") if ($debug);
	my $command ="nosave\nquery find wormpep where history AND NEXT > $i\nlist -a\nquit\n";
	&tace_it($command,$wsname);
	print "....load into db\t" if ($debug);
	&load_it($wsname,'keyset_lists') unless ($noload);
 	print "...hasta luego\n\n" if ($debug);
    } 
    print "calculated keysets of changed wormpep entries\n\n";
} 

print "le fin\n";

exit(0);

#################
## SUBROUTINES ##
#################

sub tace_it {
    my ($command,$name) = @_;
    my @results if ($debug);
    
    open (OUTPUT, ">$outpath/$name.ace") || die "Can't open the output file\n";
    open (TACE,"echo '$command' | $tace $dbpath |");
    while (<TACE>) { 
	next if /\/\//;
	next if /acedb\>/;
	next if ($_ eq "");
	s/KeySet : Answer_1/Keyset : $name/g;
	push (@results,$_) if ($debug);
	print OUTPUT;
    }
    close TACE;
    close OUTPUT;

    return (@results) if ($debug);
}
#_ end of tace_it _#


sub load_it {
    my ($name,$user) = @_;

    my $command = "pparse $outpath/$name.ace\nsave\nquit\n";

    open (TACE,"echo '$command' | $tace $dbpath |") or die "Can't make db connection\n";
    while (<TACE>) {
	print if ($debug);
    }
    close TACE;
}
#_ end of load_it _#

#######################################################################
# Help and error trap outputs                                         #
#######################################################################

sub usage {
    my $error = shift;

    if ($error eq "No_WormBase_release_number") {
        # No WormBase release number file
        print "The WormBase release number cannot be parsed\n";
        print "Check File: /wormsrv2/autoace/wspec/database.wrm\n\n";
        exit(0);
    }
    elsif ($error eq "help") {
        # Normal help menu
        exec ('perldoc',$0);
    }
}




  
__END__

=pod

=head2   NAME - make_keysets.pl

=head1 USAGE

=over 4

=item make_keysets.pl [-options]

=back

make_keysets.pl is a wee perl script to drive acedb queries via tace which
result in a keyset of objects. The keyset is then loaded back into the
database. These keysets act as a extendable version of the subclass methodology
which can be delicately configured.

make_keysets.pl  mandatory arguments:

=over 4

=item none, (but it won't do anything)

=back

make_keysets.pl OPTIONAL arguments:

=over 4

=item -expr,

=back 

CDS sequences with expression patterns

=item B<-all,> Make all of the following keysets

=item B<-rnai,> CDS sequences with RNAi experimental data

=item B<-pcr,> CDS sequences which are amplified within PCR_product objects

=item B<-go,> CDS sequences with inferred GO_terms

=item B<-touched,> CDS sequences with partial/full transcript evidence

=item B<-history XX,> wormpep accessions touched since release of wormpepXX

=item B<-noload,> create keysets but do not load them into the database

=item B<-debug,> verbose debug mode

=item B<-help,> these help pages

=back

=head1 AUTHOR (& person to blame)

=over 4

=item Dan Lawson dl1@sanger.ac.uk

=back

