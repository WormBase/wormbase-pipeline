#!/usr/local/bin/perl5.6.1                   
#
# map_features.pl
#
# by Dan Lawson
#
# This maps features to the genome based on their flanking sequence
#
# Last updated by: $Author: dl1 $                      # These lines will get filled in by cvs and helps us
# Last updated on: $Date: 2004-04-15 09:20:20 $        # quickly see when script was last changed and by whom


$|=1;
use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Feature_mapper;
use Wormbase;
use Ace;
use Getopt::Long;


my $feature;
my $clone;
my $flanking_left;
my $flanking_right;
my $coords;
my $span;

my $help;
my $debug;
my $verbose;
my $all;
my $SL1;
my $SL2;
my $polyA_site;
my $polyA_signal;
my $adhoc;

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



# ACEDB and databases
my $tace  = &tace; 
my $dbdir = "/wormsrv2/autoace";
$dbdir = "/nfs/disk100/wormpub/DATABASES/current_DB" if ($debug);

our ($WS_version) = &get_wormbase_version_name;

# assign genomic sequences to hash
#print "// Reading sequence hash\n" if ($verbose);
#my %clone2seq = &FetchData('clone2seq');

# assign Tablemaker defs for each Feature type
my %command;
$command{SL1}          = "Table-maker -p $dbdir/wquery/feature_SL1.def\nquit\n";
$command{SL2}          = "Table-maker -p $dbdir/wquery/feature_SL2.def\nquit\n"; 
$command{polyA_site}   = "Table-maker -p $dbdir/wquery/feature_polyA_site.def\nquit\n";
$command{polyA_signal} = "Table-maker -p $dbdir/wquery/feature_polyA_signal_sequence.def\nquit\n";

my %sanity;
$sanity{SL1} = 0;
$sanity{SL2} = 0;
$sanity{polyA_site} = 0;
$sanity{polyA_signal} = 6;

# queue which Feature types you want to map
my @features2map;
push (@features2map, "SL1")           if (($SL1) || ($all));
push (@features2map, "SL2")           if (($SL2) || ($all));
push (@features2map, "polyA_signal")  if (($polyA_signal) || ($all));
push (@features2map, "polyA_site")    if (($polyA_site) || ($all));

my $mapper      = Feature_mapper->new($dbdir);

# main loop

foreach my $query (@features2map) {

    print "// mapping $query features\n" if ($verbose);

    open (OUTPUT, ">/wormsrv2/autoace/FEATURES/${WS_version}_feature_${query}.ace") or die "Failed to open output file\n" unless ($adhoc);
    open (ERRORS, ">/wormsrv2/autoace/FEATURES/${WS_version}_feature_${query}.err") or die "Failed to open error file\n" unless ($adhoc);

    if ($adhoc) {
	open (TACE, "<$adhoc") or die "Failed to open input file: $adhoc\n";
	print "// Opening a file for input: $adhoc\n";
    }
    else {
	open (TACE, "echo '$command{$query}' | $tace $dbdir | ");
    }
    while (<TACE>) {
	# when it finds a good line
	if (/^\"(\S+)\"\s+\"(\S+)\"\s+\"(\S+)\"\s+\"(\S+)\"/) {
	    ($feature,$clone,$flanking_left,$flanking_right) = ($1,$2,$3,$4);

#	    $coords = &Map_feature($clone2seq{$clone},$flanking_left,$flanking_right);
	    my @coords = $mapper->map_feature($clone,$flanking_left,$flanking_right);

	    if ($coords[1] eq "") {
		print "// do it again\n";
		my $rev_left     = &DNA_string_reverse2($flanking_left);
		my $rev_right    = &DNA_string_reverse2($flanking_right);
		@coords = $mapper->map_feature($clone,$rev_left,$rev_right);
	    }


	    if ($coords[1] > $coords[2]) {
		$span = $coords[1] - $coords[2] -1 if ($SL1 || $SL2 || $polyA_site);  # +1 for TSL tight junction
		$span = $coords[1] - $coords[2] +1 if ($polyA_signal);                # -1 for polyA_signal

	    }
	    else {
		$span = $coords[2] - $coords[1] -1 if ($SL1 || $SL2 || $polyA_site) ;  # +1 for TSL tight junction
		$span = $coords[2] - $coords[1] +1 if ($polyA_signal);                 # -1 for polyA_signal
	    }

	    print "// $feature\n[$flanking_left]\n[$flanking_right]\n" if ($adhoc);
	    print "//$feature maps to $clone $coords[1] -> $coords[2], feature span is $span bp\n";

#	    my @coords = &MapFeature2($feature,$clone,$flanking_left,$flanking_right,$clone2seq{$clone});
#	    
#	    print  "Mapping $feature [$flanking_left | $flanking_right] to $clone: $coords[0] -> $coords[1], span $coords[2]\n";
#	    
	    if ($span == $sanity{$query}) {
		print OUTPUT "//$feature maps to $clone $coords[1] -> $coords[2], feature span is $span bp\n" if ($verbose);
		print OUTPUT "\nSequence : \"$clone\"\n" unless ($adhoc);
		print OUTPUT "Feature_object $feature $coords[1]  $coords[2]\n\n" unless ($adhoc);
		print "$feature maps to $clone $coords[1] -> $coords[2], feature span is $span bp\n" if ($adhoc);
	    }
	    else {
		print ERRORS "//$feature maps to $clone $coords[1] -> $coords[2], feature span is $span bp\n" unless ($adhoc);
	    }
	} #_ if match line
    }
    close TACE;
}

close OUTPUT unless ($adhoc);
close ERRORS unless ($adhoc);

exit;


sub DNA_string_reverse2 {
    my $revseq = reverse shift;
    $revseq =~ tr/A/x/;
    $revseq =~ tr/T/A/;
    $revseq =~ tr/x/T/;
    $revseq =~ tr/G/x/;
    $revseq =~ tr/C/G/;
    $revseq =~ tr/x/C/;
    return ($revseq);
}    
