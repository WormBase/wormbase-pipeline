#!/usr/local/bin/perl5.8.0                   
#
# map_features.pl
#
# by Dan Lawson
#
# This maps features to the genome based on their flanking sequence.
# Uses Ant's Feature_mapper.pm module
#
#
# Last updated by: $Author: dl1 $                      # These lines will get filled in by cvs and helps us
# Last updated on: $Date: 2004-08-12 12:23:30 $        # quickly see when script was last changed and by whom


$|=1;
use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Feature_mapper;
use Wormbase;
use Ace;
use Getopt::Long;


my ($feature, $clone, $flanking_left, $flanking_right, $coords, $span);

my $help;                    # Help menu
my $debug;                   # Debug mode 
my $verbose;                 # Verbose mode
my $all;                     # Do all the following features:
my $SL1;                     #  SL1 trans-splice leader acceptors
my $SL2;                     #  SL2 trans-splice leader acceptors
my $polyA_site;              #  polyA_site
my $polyA_signal;            #  polyA_signal
my $adhoc;                   # Run against a file, output to screen
my $start;
my $stop;


GetOptions (
	    "all"          => \$all,
	    "SL1"          => \$SL1,
	    "SL2"          => \$SL2,
	    "polyA_site"   => \$polyA_site,
	    "polyA_signal" => \$polyA_signal,
	    "adhoc=s"      => \$adhoc,
            "debug=s"      => \$debug,
            "verbose"      => \$verbose,
	    "help"         => \$help
           );

# Help pod if needed
&usage(0) if ($help);

#######################
# ACEDB and databases #
#######################

our $tace   = &tace; 
our $outdir = "/wormsrv2/autoace/FEATURES";
our $dbdir  = "/wormsrv2/autoace";
$dbdir = "/nfs/disk100/wormpub/DATABASES/current_DB" if ($debug);


# WS version for output files

our ($WS_version) = &get_wormbase_version_name;

# coordinates for Feature_mapper.pm module

my $mapper      = Feature_mapper->new($dbdir);

# sanity checks for the length of feature types

my %sanity = (
	      'SL1'          => 0,
	      'SL2'          => 0,
	      'polyA_site'   => 0,
	      'polyA_signal' => 6
	      );

# queue which Feature types you want to map

my @features2map;
push (@features2map, "SL1")           if (($SL1) || ($all));
push (@features2map, "SL2")           if (($SL2) || ($all));
push (@features2map, "polyA_signal")  if (($polyA_signal) || ($all));
push (@features2map, "polyA_site")    if (($polyA_site) || ($all));

#############
# main loop #
#############

foreach my $query (@features2map) {

    print "// mapping $query features\n" if ($verbose);
    
    # open output files

    open (OUTPUT, ">$outdir/${WS_version}_feature_${query}.ace") or die "Failed to open output file\n" unless ($adhoc);
    open (ERRORS, ">$outdir/${WS_version}_feature_${query}.err") or die "Failed to open error file\n" unless ($adhoc);

    # start tace session for input data (or find file for adhoc run)
    if ($adhoc) {
	open (TACE, "<$adhoc") or die "Failed to open input file: $adhoc\n";
	print "// Opening a file for input: $adhoc\n" if ($verbose);
    }
    else {
	open (TACE, "echo 'Table-maker -p $dbdir/wquery/feature_${query}.def\nquit\n' | $tace $dbdir | ");
    }
    while (<TACE>) {

	# when it finds a good line
	if (/^\"(\S+)\"\s+\"(\S+)\"\s+\"(\S+)\"\s+\"(\S+)\"/) {
	    ($feature,$clone,$flanking_left,$flanking_right) = ($1,$2,$3,$4);

	    my @coords = $mapper->map_feature($clone,$flanking_left,$flanking_right);

	    $start = $coords[1];
	    $stop  = $coords[2];

	    # munge returned coordinates to get the span of the mapped feature

	    # Deal with polyA_signal features
	    if ($polyA_signal) {
		if ($start < $stop) {
		    $start++;
		    $stop--;
		    $span = $stop - $start + 1;
		}
		else {
		    $start--;
		    $stop++;
		    $span = $start - $stop + 1;
		}
	    }
	    # else deal with butt-ended features (e.g. SL1, SL2 & polyA_site)
	    elsif ($start > $stop) {
		$span = $start - $stop - 1;
	    }
	    else {
		$span = $stop - $start - 1;
	    }

	    # check feature span is sane
	    
	    if ($span == $sanity{$query}) {
		
		if ($adhoc) {
		    print "$feature maps to $clone $start -> $stop, feature span is $span bp\n";
		}
		else {
		    print OUTPUT "//$feature maps to $clone $start -> $stop, feature span is $span bp\n";
		    print OUTPUT "\nSequence : \"$clone\"\n";
		    print OUTPUT "Feature_object $feature $start  $stop\n\n";
		}
	    }
	    else {
		print ERRORS "// $feature maps to $clone $start -> $stop, feature span is $span bp\n" unless ($adhoc);
	    }
	} #_ if match line
    }
    close TACE;
}

close OUTPUT unless ($adhoc);
close ERRORS unless ($adhoc);

###############
# hasta luego #
###############

exit(0);


#######################################################################
# Help and error trap outputs                                         #
#######################################################################

sub usage {
    my $error = shift;

    if ($error == 1) {
    }
    elsif ($error == 0) {
        # Normal help menu
        exec ('perldoc',$0);
    }
}



__END__

=pod

=head2 NAME - map_features.pl

=head1 USAGE

=over 4

=item map_features.pl [-options]

=back

map_features.pl mandatory arguments:

=over 4

=item none

=back

map_features.pl optional arguments:

=over 4

=item -all 

map all of the following feature types:

=item -SL1

map SL1 trans-splice leader acceptor sites (2 bp feature)

=item -SL2

map SL2 trans-splice leader acceptor sites (2 bp feature)

=item -polyA_site

map polyA attachement sites (2 bp feature)

=item -polyA_signal

map polyA signal sequence sites (6 bp feature)

=item -debug <user> 

Queries current_DB rather than autoace

=item -adhoc <file> 

Queries a flatfile for the Feature data. Flatfile format is
"<feature_name>"  "<clone>" "<flanking_sequence_left>" "<flanking_sequence_right>" 

=item -help      

This help page

=cut
