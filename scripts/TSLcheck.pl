#!/usr/local/bin/perl5.8.0 -w
#
# TSLcheck.pl
#
# by Dan
#
# checks whether TSL acceptor sites overlap with exons.
# Will discard matches to isoforms (Warning).
#
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2005-12-19 13:32:45 $


use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;

#######################################
# command-line options                #
#######################################
 
my $database;       # specify release to read GFF split files from
my $debug;          # debug mode,
my $help;           # help mode, show perldoc
my $verbose;        # verbose mode, extra output to screen
my $test;
my $store;
my $wormbase;
 
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "database=s" => \$database,
	    "store"      => \$store,
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





# database options
my $datadir  = $wormbase->gff_splits;  # AUTOACE GFF SPLIT DIR
my $ace_dir  = $wormbase->autoace;     # AUTOACE DATABASE DIR


if ($database) {
    $datadir = $ace_dir . "/GFF_SPLITS/$database"; 
}

# output options
my $outdir      = $ace_dir . "/CHECKS";


# get clone2 centre data
my %clone2centre = $wormbase->FetchData('clone2centre');       # Clone => (HX|RW) centre designation

# vars
my @chromosomes  = qw(I II III IV V X);
my %overlap;                                         # No of overlaps per chromosome       
my %isoform;                                         # No of isoform overlaps per chromosome
my %active;                                          # No of overlaps needing to be look at
our $chromosome;
my $cds; 
my $line;


#############
# Main Loop #
#############

foreach $chromosome (@chromosomes) {
 
    # open output files
    open (CAM, ">$outdir/CHROMOSOME_${chromosome}.overlapping_TSL_cam");
    open (STL, ">$outdir/CHROMOSOME_${chromosome}.overlapping_TSL_stl");
    
    # input lines are processed GFF (see below)
    # CHROMOSOME_X    curated exon    2089999 2090139 .       +       .       CDS "F49H12.1" intersect(CHROMOSOME_X.TSL_site.gff)=(21)

    # check that the input files exist
    &usage("No TSL file")  unless (-e "$datadir/CHROMOSOME_${chromosome}.TSL_site.gff");
    &usage("No exon file") unless (-e "$datadir/CHROMOSOME_${chromosome}.exon.gff");

    # cycle through overlap file and assign to centre
    open (IN, "gff_overlap -minfrac1 1 $datadir/CHROMOSOME_${chromosome}.TSL_site.gff $datadir/CHROMOSOME_${chromosome}.CDS.gff |");
    while (<IN>) {
	
	$line = $_;
	print if ($debug);

	# get CDS name
	($cds) = $line =~ (/CDS \"(\S+)\./);
	$overlap{$chromosome}++;                                          # increment overlap count

	# discard isoforms 
	if ($line =~ /[a-z]\" /) {
	    print "// discard as an isoform match\n"if ($verbose);
	    $isoform{$chromosome}++;                                      # increment isoform count
	    $log->write_to("WARNING: Isoform match $line") if ($verbose);
	    next;
	}

	print "// processing $cds [$clone2centre{$cds}]\n" if ($verbose);
	$active{$chromosome}++;                                           # increment active count
	
	if ($clone2centre{$cds} eq "HX") {
	    print CAM $line;
	}
	elsif ($clone2centre{$cds} eq "RW") {
	    print STL $line;
	}
	
    }
    close IN;
    
    # close output filehandles
    close CAM;
    close STL;
    
    # write output summary to mail
    $log->write_to("CHROMOSOME $chromosome\t: $overlap{$chromosome} overlaps of which $isoform{$chromosome} are to isoforms  \t($active{$chromosome})\n");

}

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);    
 

sub usage {
    my $error = shift;
    
    if ($error eq "Help") {
        # Help menu
	exec ('perldoc',$0);
    }
    elsif ($error eq "No TSL file") {
	print "No TSL_site GFF split file exists for chromosome $chromosome\n";
	exit;
    }
    elsif ($error eq "No exon file") {
	print "No exon GFF split file exists for chromosome $chromosome\n";
	exit;
    }
    elsif ($error == 0) {
	# Normal help menu
	exec ('perldoc',$0);
    } 
}



__END__
 
=pod
 
=head1 NAME - TSLcheck.pl

=head2 USAGE 

TSLcheck.pl maps the overlaps of the TSL features to curated exons based on the GFF files. Currently, any overlap to exons from alternate
isoforms are discarded (although counted). The output is 12 GFF files (two per chromosome split based on the centre).

TSLcheck.pl arguements:

=over 4

=item -database, process a named GFF_SPLIT directory rather than the default autoace current build (e.g. -database WS100) 

=item -debug <user>, Debug mode - prints all GFF lines
 
=item -verbose, Verbose mode toggle on extra command line output
 
=item -help, these help pages
 
=back

=head1 AUTHOR

Dan Lawson (dl1@sanger.ac.uk)

=cut
