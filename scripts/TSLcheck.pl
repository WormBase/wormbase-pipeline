#!/usr/local/bin/perl5.8.0 -w
#
# TSLcheck.pl
#
# by Dan
#
# checks whether TSL acceptor sites overlap with exons.
# Will discard matches to isoforms (Warning).
#
# Last updated by: $Author: dl1 $
# Last updated on: $Date: 2004-09-14 15:47:51 $

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts"  : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
 

#######################################
# command-line options                #
#######################################
 
my $database;       # specify release to read GFF split files from
my $debug;          # debug mode,
my $help;           # help mode, show perldoc
my $verbose;        # verbose mode, extra output to screen
 
GetOptions( "debug=s"    => \$debug,
            "database=s" => \$database,
            "help"       => \$help,
            "verbose"    => \$verbose
	    );
                                              
# check command line options

# Help pod if needed
&usage("Help") if ($help);

# database options
my $datadir     = "/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS"; 

if ($database) {
    $datadir = "/wormsrv2/autoace/GFF_SPLITS/$database"; 
}

# output options
my $outdir      = "/wormsrv2/autoace/CHECKS";


# get clone2 centre data
my %clone2centre = &FetchData('clone2centre');       # Clone => (HX|RW) centre designation

# vars
my @chromosomes  = qw(I II III IV V X);
my %overlap;                                         # No of overlaps per chromosome       
my %isoform;                                         # No of isoform overlaps per chromosome
my %active;                                          # No of overlaps needing to be look at
our $chromosome;
my $cds; 
my $line;

# create log file, open output file handles
my $log = Log_files->make_build_log();


#############
# Main Loop #
#############

foreach $chromosome (@chromosomes) {
 
    # open output files
    open (CAM, ">$outdir/CHROMOSOME_${chromosome}.overlapping_TSL_cam.gff");
    open (STL, ">$outdir/CHROMOSOME_${chromosome}.overlapping_TSL_stl.gff");
    
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

$log->mail("dl1\@sanger.ac.uk", "BUILD REPORT: TSLcheck.pl");

# hasta luego
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
