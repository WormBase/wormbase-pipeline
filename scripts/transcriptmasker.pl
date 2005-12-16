#!/usr/local/bin/perl5.8.0 -w
#
# transcriptmasker.pl
#
# masks out ?Feature data spans in mRNA/ESTs prior to the BLAT analysis
# (essentially to remove TSL and polyA sequences)

# 031023 dl1

# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2005-12-16 11:18:55 $

#################################################################################
# Initialise variables                                                          #
#################################################################################

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Ace;
use Getopt::Long;
use Carp;
use File::Path;
use Storable;

$|=1;

##############################
# command-line options       #
##############################

my $debug;              # debug mode
<<<<<<< transcriptmasker.pl
my $report;             # simple report mode
=======
my $database;
my $test;
my $store;
my $verbose;            # verbose mode
>>>>>>> 1.9.4.1
my $help;               # Help/Usage page
my $mrna;               # mRNA data
my $ncrna;              # ncRNA data
my $est;                # EST data
my $ost;                # OST data
my $all;                # all of the above

GetOptions (
	    "all"            => \$all,
	    "mrna"           => \$mrna,
	    "ncrna"          => \$ncrna,
	    "est"            => \$est,
	    "ost"            => \$ost,
            "debug:s"        => \$debug,
<<<<<<< transcriptmasker.pl
	    "report:s"       => \$report,
            "help"           => \$help,
=======
	    "verbose"        => \$verbose,
            "help|h"         => \$help,
	    "database:s"     => \$database,
	    "test"           => \$test,
	    "store:s"        => \$store
>>>>>>> 1.9.4.1
	    );

# Help pod if needed
&usage("Help") if ($help);

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}
if ($report) {
    print "// REPORT : \"$report\"\n\n";
    ($maintainers = $report . '\@sanger.ac.uk');
}

my $log = Log_files->make_build_log($wormbase);


# datafiles for input
my $EST_dir = $wormbase->wormpub."/analysis/ESTs";
our %datafiles = (
		  "mrna"  => "elegans_mRNAs",
		  "ncrna" => "elegans_ncRNAs",
		  "est"   => "elegans_ESTs",
		  "ost"   => "elegans_OSTs"
		  );

# valid Feature_data methods
our @valid_methods = (
		      "TSL",                                          # TSL data
		      "polyA",                                        # polyA tail sequences
		      "poly_nucleotide",                              # poly nucleotide sequences (e.g. 5' AAAAAAA)
		      "BLAT_discrepency",                             # General method for BLAT problems
		      "vector"                                        # Vector sequences
		      );

our $valid_methods = join(' ', @valid_methods);

# transcript accessions to names from a hash in common data

<<<<<<< transcriptmasker.pl
print "// Reading EST_names.dat hash\n\n" if ($debug);
our %EST_name = &FetchData('NDBaccession2est');
print "// Finished reading EST_names.dat hash\n\n" if ($debug);
=======
print "// Reading EST_names.dat hash\n\n" if ($verbose);
our %EST_name = $wormbase->FetchData('NDBaccession2est');
print "// Finished reading EST_names.dat hash\n\n" if ($verbose);
>>>>>>> 1.9.4.1

# which database
<<<<<<< transcriptmasker.pl
my $dbdir = "/wormsrv2/autoace";
if ($debug) {$dbdir = "/nfs/disk100/wormpub/camace_dl1";}    # Use dan's split for the debug option
if ($report) {$dbdir = "/nfs/disk100/wormpub/camace_dl1";}    # Use dan's split for the debug option
my $tace  = &tace;                                    # tace executable path
=======
$database = $wormbase->autoace unless $database;

my $blat_dir = "$database/BLAT";
mkpath($blat_dir) unless (-e $blat_dir );

my $tace  = $wormbase->tace;                                    # tace executable path

>>>>>>> 1.9.4.1
my $acc;                                              # accession for the entry
my $id;                                               # id for the entry
my $seq;                                              # raw sequence for the entry
my %sequence;                                         # 
my @features;                                         # list of feature_data objects for the sequence
my $feature;                                          #
my $no_features;
my ($type,$method,$start,$stop,$length,$remark);      #
my ($cut_to,$cut_from,$cut_length,$newseq);           #
my $seqmasked;                                        #
my $seqlength;                                        #
my $masked;                                           # No of entries masked
my $ignored;                                          # No of entries ignored
my $ignore;

# which data file to parse
<<<<<<< transcriptmasker.pl
$masked = &MaskSequence($datafiles{mrna})  if ($mrna || $all);
print LOG &runtime, ": masked $masked mRNA sequences\n"   if ($mrna || $all);
print LOG &runtime, ": ignored $ignored mRNA sequences\n" if ($mrna || $all);
=======
$masked = &MaskSequence($datafiles{mrna}) if ($mrna || $all);
$log->write_to($wormbase->runtime." : masked $masked mRNA sequences\n") if ($mrna || $all);
>>>>>>> 1.9.4.1

$masked = &MaskSequence($datafiles{ncrna}) if ($ncrna || $all);
<<<<<<< transcriptmasker.pl
print LOG &runtime, ": masked $masked ncRNA sequences\n"   if ($ncrna || $all);
print LOG &runtime, ": ignored $ignored ncRNA sequences\n" if ($ncrna || $all);
=======
$log->write_to($wormbase->runtime." : masked $masked ncRNA sequences\n") if ($ncrna || $all);

$masked = &MaskSequence($datafiles{est})  if ($est || $all);
$log->write_to($wormbase->runtime." : masked $masked EST sequences\n") if ($est || $all);
>>>>>>> 1.9.4.1

<<<<<<< transcriptmasker.pl
$masked = &MaskSequence($datafiles{est})   if ($est || $all);
print LOG &runtime, ": masked $masked EST sequences\n"   if ($est || $all);
print LOG &runtime, ": ignored $ignored EST sequences\n" if ($est || $all);

$masked = &MaskSequence($datafiles{ost})   if ($ost || $all);
print LOG &runtime, ": masked $masked OST sequences\n"   if ($ost || $all);
print LOG &runtime, ": ignored $ignored OST sequences\n" if ($ost || $all);
=======
$masked = &MaskSequence($datafiles{ost})  if ($ost || $all);
$log->write_to($wormbase->runtime." : masked $masked OST sequences\n") if ($ost || $all);
>>>>>>> 1.9.4.1

$log->write_to("\n=============================================\n\n");

#########################################
# hasta luego                           #
#########################################

$log->mail;

exit(0);


############################################################
######################## Subroutines #######################
############################################################

#_ MaskSequence -#
# 
# pass type of transcript data to be masked (e.g. mRNA, EST etc)

sub MaskSequence {
    my $data   = shift;
    my $masked = 0;
<<<<<<< transcriptmasker.pl
    my $ignored = 0;

=======
    my $ignore ;

>>>>>>> 1.9.4.1
    # connect to database
    print  "Opening database for masking $data ..\n" if ($debug);
    my $db = Ace->connect(-path=>$database,
                          -program =>$tace) || $log->log_and_die("Connection failure: ".Ace->error."\n");

    # set input record seperator
    $/ = ">";

    # assign output file
<<<<<<< transcriptmasker.pl
    if ( ($debug) || ($report) ) {
	open (OUTPUT, ">${data}.testmasked") || die "ERROR: Can't open output file: '${data}.testmasked'";
    }
    else {
	open (OUTPUT, ">${data}.masked") || die "ERROR: Can't open output file: '${data}.masked'";
    }
=======
    open (OUTPUT, ">$blat_dir/$data") || $log->log_and_die("ERROR: Can't open output file: $blat_dir/$data");
>>>>>>> 1.9.4.1

    # input file loop structure
    my $skip = 1;
    
    open (INPUT, "<$EST_dir/$data")     || $log->log_and_die("ERROR: Can't open input file: $EST_dir/$data");
    while (<INPUT>) {
	chomp;
	next if ($_ eq "");                 # catch empty lines
	if ($skip == 1) {                   # skip first bad entry
	    $skip--;
	    next;
	}
	
	if (/^(\S+)\s+\S+.+\n/) {           # deal with accessions {$acc} and WormBase internal names {$id}
	    $acc = $1;
	    if (defined $EST_name{$acc}) {
		$id  = $EST_name{$acc};
	    }
	    else {
#		print "// ERROR: No accession-id connection for this sequence [$acc]\n" if ($debug);
		$EST_name{$acc} = $acc;
		$id = $acc;
	    }
	}
	$seq = "$'";                                  # assign the rest of the string to $seq
	$seq =~ s/[^gatcn]//g;                        # remove non-gatcn characters (i.e. newlines)
	$seqmasked = $seq;                            # copy sequence to masked file and
	$seqlength = length ($seq);                   # calculate the length of the sequence

	# fetch the sequence object from the database. push the feature_data objects to memory
	my $obj = $db->fetch(Sequence=>$id);
	if (!defined ($obj)) {
	    print "ERROR: Could not fetch sequence $id \n\n" if ($debug);
	    next;
	}
	
<<<<<<< transcriptmasker.pl
	# fetch the feature_data for the database object
	@features = $obj->Feature_data;                          # ?Method tag (i.e. SL1/SL2/polyA)
	$no_features = scalar (@features);                       # Number of features to parse
	
	if ($no_features == 0) {
	    print "WARNING: No Features to parse from $id\n\n" if ($debug);
=======
	@features = $obj->Feature_data(1);

	if ( scalar (@features) == 0) {
	    print "ERROR: No Features to parse \n\n" if ($debug);
>>>>>>> 1.9.4.1
	}
	else {
	    for (my $i=0; $i < $no_features; $i++) {             # loop through each attached ?Feature_data
		($type) = $features[$i] =~ (/^$id\:(\S+)/);
		
		next unless ($valid_methods =~ /$type/);         # only mask valid Feature_data methods

		$method = $features[$i]->Feature(1);
		$start  = $features[$i]->Feature(2);
		$stop   = $features[$i]->Feature(3);
		
<<<<<<< transcriptmasker.pl
		unless ((defined $start) && (defined $stop)) {
		    print LOG "// ERROR: No coordinates for the $method feature in $id\n";
		    next;
		}

		print "// Parse $features[$i]\t$start - $stop\n" if ($debug);

		$cut_to     = $start - 1;                        # manipulations for clipping 
		$cut_from   = $stop;
		$cut_length = $stop - $start + 1;
		    
		if ($cut_to < 0 ) {$cut_to = 0;}                 # fudge to ensure non-negative clipping coords
		    
		print "// $acc [$id] $method:" . substr($seq,$cut_to,$cut_length) . " [$start - $stop]\n\n" if ($debug);
		print "// $acc [$id] $method:" . substr($seq,$cut_to,$cut_length) . " [$start - $stop]\n" if ($report);
		$newseq = substr($seqmasked,0,$cut_to) . ('n' x $cut_length)  . substr($seqmasked,$cut_from);
		$seqmasked = $newseq;
=======
		$type  = $obj->Feature_data->Feature(1);         # Feature type (e.g. SL1,SL2,polyA)
		if (defined($type)) {
		    $start = $obj->Feature_data->Feature(2);         # start coord
		    $stop  = $obj->Feature_data->Feature(3);         # stop coord

		    $cut_to     = $start - 1;                        # manipulations for clipping 
		    $cut_from   = $stop;
		    $cut_length = $stop - $start + 1;

		    if ($cut_to < 0 ) {$cut_to = 0;}                 # fudge to ensure non-negative clipping coords

		    print "$acc [$id]: '$type' $start -> $stop [$cut_to : $cut_from ($cut_length)]\n" if ($debug);
		    print "// # $acc [$id] $type:" . (substr($seq,$cut_to,$cut_length)) . " [$start - $stop]\n\n" if ($verbose);
		    $newseq = (substr($seqmasked,0,$cut_to)) . ('n' x $cut_length)  . (substr($seqmasked,$cut_from));
		    $seqmasked = $newseq;
		}
>>>>>>> 1.9.4.1
	    }
	    
	    # increment count of sequences masked
	    $masked++;

	}
	
<<<<<<< transcriptmasker.pl
	# Is the Ignore tag set?

	$ignore = $obj->Ignore;        # Ignore tag set?
	unless (defined $ignore) {
	    # output masked sequence
	    print OUTPUT ">$acc $id\n$seqmasked\n\n";
	}
	else {
	    print "\n// Ignore tag set for $acc $id\n\n" if ($debug);
	    print "\n// Ignore tag set for $acc $id\n\n" if ($report);
	    $ignored++;
	}	
	undef $ignore;
=======
	# Is the Ignore tag set?

	my $ignore = $obj->Ignore;        # Ignore tag set?
	unless (defined $ignore) {
	    # output masked sequence
	    print OUTPUT ">$acc $id\n$seqmasked\n";
	}
	else {
	    print "\n// Ignore tag set for $acc $id\n\n" if ($verbose);
	}	
	undef $ignore;
>>>>>>> 1.9.4.1
	
	# close object
	$obj->DESTROY();
	
    }
    close INPUT;
    $/ = "\n";

    close OUTPUT;

    return ($masked,$ignored);
}
#_ end MaskSequence _#



###############################################################

sub usage {
    my $error = shift;
    
    if ($error eq "Help") {
	# Help menu
	 exec ('perldoc',$0);
     }
    elsif ($error == 0) {
        # Normal help menu
        exec ('perldoc',$0);
    }
}


__END__

=pod

=head2   NAME - transcriptmasker.pl

=head1 USAGE

=over 4

=item transcriptmasker.pl [-options]

=back

transcriptmasker.pl takes the fasta file generated using SRS and queries autoace
for ?Feature_data classes. It will mask out any annotated features (currently TSL
acceptor sites, and polyA+ tails) with N's. These masked files are then used by 
the BLAT scripts to map back to the genome.

transcriptmasker.pl mandatory arguments:

=over 4

=item none,

=back

transcriptmasker.pl OPTIONAL arguments:

=over 4

=item -all, process EST, mRNA & OST datafiles 

=item -est, process EST datafile 

=item -ost, process OST datafile 

=item -mrna, process mRNA datafile 

=item -debug <user>, Debug mode - log file only goes to named user

=item -help, these help pages

=back

=cut
