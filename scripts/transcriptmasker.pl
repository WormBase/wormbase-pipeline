#!/usr/local/bin/perl5.8.0 -w
#
# transcriptmasker.pl
#
# masks out ?Feature data spans in mRNA/ESTs prior to the BLAT analysis
# (essentially to remove TSL and polyA sequences)

# 031023 dl1

# Last edited by: $Author: krb $
# Last edited on: $Date: 2004-08-03 14:16:26 $

#################################################################################
# Initialise variables                                                          #
#################################################################################

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use IO::Handle;
use Ace;
use Getopt::Long;
use Carp;

$|=1;

##############################
# command-line options       #
##############################

my $maintainers = "All";
our $log;

my $debug;              # debug mode
my $verbose;            # verbose mode
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
	    "verbose"        => \$verbose,
            "help"           => \$help,
            "h"              => \$help
	    );

# Help pod if needed
&usage("Help") if ($help);

# Use debug mode?
if ($debug) {
    print "// DEBUG : \"$debug\"\n\n";
    ($maintainers = $debug . '\@sanger.ac.uk');
}

&create_log_files;


# datafiles for input
our %datafiles = (
		  "mrna"  => "/nfs/disk100/wormpub/analysis/ESTs/elegans_mRNAs",
		  "ncrna" => "/nfs/disk100/wormpub/analysis/ESTs/elegans_ncRNAs",
		  "est"   => "/nfs/disk100/wormpub/analysis/ESTs/elegans_ESTs",
		  "ost"   => "/nfs/disk100/wormpub/analysis/ESTs/elegans_OSTs"
		  );

# transcript accessions to names from a hash in common data

print "// Reading EST_names.dat hash\n\n" if ($verbose);
our %EST_name = &FetchData('NDBaccession2est');
print "// Finished reading EST_names.dat hash\n\n" if ($verbose);

# which database
my $dbdir = "/wormsrv2/autoace";
my $tace  = &tace;                                    # tace executable path

my $acc;                                              # accession for the entry
my $id;                                               # id for the entry
my $seq;                                              # raw sequence for the entry
my %sequence;                                         # 
my @features;                                         # list of feature_data objects for the sequence
my $feature;                                          #
my ($type,$start,$stop,$length,$remark);              #
my ($cut_to,$cut_from,$cut_length,$newseq);           #
my $seqmasked;                                        #
my $seqlength;                                        #
my $masked;                                           # No of entries masked

# which data file to parse
$masked = &MaskSequence($datafiles{mrna}) if ($mrna || $all);
print LOG &runtime, ": masked $masked mRNA sequences\n" if ($mrna || $all);

$masked = &MaskSequence($datafiles{ncrna}) if ($ncrna || $all);
print LOG &runtime, ": masked $masked ncRNA sequences\n" if ($ncrna || $all);

$masked = &MaskSequence($datafiles{est})  if ($est || $all);
print LOG &runtime, ": masked $masked EST sequences\n" if ($est || $all);

$masked = &MaskSequence($datafiles{ost})  if ($ost || $all);
print LOG &runtime, ": masked $masked OST sequences\n" if ($ost || $all);

print LOG "\n";
print LOG "=============================================\n";
print LOG "\n";
close LOG;

#########################################
# hasta luego                           #
#########################################

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
  
    # connect to database
    print  "Opening database for masking $data ..\n" if ($debug);
    my $db = Ace->connect(-path=>$dbdir,
                          -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

    # set input record seperator
    $/ = ">";

    # assign output file
    if ($debug) {
	open (OUTPUT, ">${data}.testmasked") || die "ERROR: Can't open output file: '${data}.testmasked'";
    }
    else {
	open (OUTPUT, ">${data}.masked") || die "ERROR: Can't open output file: '${data}.masked'";
    }

    # input file loop structure
    my $skip = 1;
    
    open (INPUT, "<$data")     || die "ERROR: Can't open input file: '$data'";
    while (<INPUT>) {
	chomp;
	next if ($_ eq "");                 # catch empty lines
	if ($skip == 1) {$skip--;next;}     # skip first bad entry
	
	if (/^(\S+)\s+\S+.+\n/) {           # deal with accessions {$acc} and WormBase internal names {$id}
	    $acc = $1;
	    if (defined $EST_name{$acc}) {
		$id  = $EST_name{$acc};
	    }
	    else {
		print "// ERROR: No accession-id connection for this sequence [$acc]\n" if ($verbose || $debug);
		$EST_name{$acc} = $acc;
		$id = $acc;
	    }
	}
	$seq = "$'";                        # assign the rest of the string to $seq
	$seq =~ s/[^gatcn]//g;              # remove non-gatcn characters (i.e. newlines)
	$seqmasked = $seq;                  # copy sequence to masked file and
	$seqlength = length ($seq);         # calculate the length of the sequence
	
	print "-> Parsing $acc [$id]\n" if ($verbose);	

	# fetch the sequence object from the database. push the feature_data objects to memory
	my $obj = $db->fetch(Sequence=>$id);
	if (!defined ($obj)) {
	    print "ERROR: Could not fetch sequence $id \n" if ($debug);
	    next;
	}
	
	@features = $obj->Feature_data(1);
	
	if ( scalar (@features) == 0) {
	    print "ERROR: No Features to parse \n" if ($debug);
	}
	else {
	    foreach $feature (@features) {
		print "\n// parse $feature\n" if ($verbose);
		
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
	    }

	    # increment count of sequences masked
	    $masked++;
	    
	}
	
	# output masked sequence
	print OUTPUT ">$acc $id\n$seqmasked\n";
	
	# close object
	$obj->DESTROY();
	
    }
    close INPUT;
    $/ = "\n";

    close OUTPUT;

    return ($masked);
}
#_ end MaskSequence _#

###############################################################

sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name    =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate     = `date +%y%m%d`; chomp $rundate;
  $log            = "/wormsrv2/logs/$script_name.$rundate.$$";

  open (LOG, ">$log") or die "cant open $log";
  print LOG "=============================================\n";
  print LOG "$script_name\n";
  print LOG "started at ",`date`,"\n";
  print LOG "=============================================\n";
  print LOG "\n";

}


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

=item -verbose, Verbose mode toggle on extra command line output

=item -help, these help pages

=back

=cut
