#!/usr/local/bin/perl5.8.0 -w

# transcriptmasker.pl
#
# masks out ?Feature data spans in mRNA/ESTs prior to the BLAT analysis
# (essentially to remove TSL and polyA sequences)

# 031023 dl1

# Last edited by: $Author: dl1 $
# Last edited on: $Date: 2003-10-31 16:29:45 $

#################################################################################
# Initialise variables                                                          #
#################################################################################

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Ace;
use Data::Dumper;
use Getopt::Long;
use Carp;

##############################
# command-line options       #
##############################

my $maintainers = "All";
our $log;

my $debug;              # debug mode
my $verbose;            # verbose mode
my $help;               # Help/Usage page
my $mrna;               # mRNA data
my $est;                # EST data
my $ost;                # OST data

GetOptions (
	    "mrna"           => \$mrna,
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
		  "mrna" => "/nfs/disk100/wormpub/analysis/ESTs/elegans_mRNAs",
		  "est"  => "/nfs/disk100/wormpub/analysis/ESTs/elegans_ESTs",
		  "ost"  => "/nfs/disk100/wormpub/analysis/ESTs/elegans_OSTs"
		  );

# transcript accessions to names from a hash
# 
# !! need to deal with making this .dat file etc 
#

our %EST_name;
&recreate_hashes("/nfs/disk100/wormpub/analysis/ESTs/EST_names.dat");

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


# which data file to parse
my $data;
$data = $datafiles{mrna} if ($mrna);
$data = $datafiles{est}  if ($est);
$data = $datafiles{ost}  if ($ost);

# connect to database

print  "Opening database ..\n" if ($debug);
my $db = Ace->connect(-path=>$dbdir,
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

# set input record seperator

$/ = ">";

# assign output file

open (OUTPUT, ">${data}.masked") || die "ERROR: Can't open output file: '${data}.masked'";

# input file loop structure

open (INPUT, "<$data")     || die "ERROR: Can't open input file: '$data'";
while (<INPUT>) {
    chomp;
    next if ($_ eq "");                 # catch empty lines
    
    if (/^(\S+)\s+\S+.+\n/) {           # deal with accessions {$acc} and WormBase internal names {$id}
	$acc = $1;
	$id  = $EST_name{$acc};
    }
    $seq = "$'";                        # assign the rest of the string to $seq
    $seq =~ s/[^gatcn]//g;              # remove non-gatcn characters (i.e. newlines)
    $seqmasked = $seq;                  # copy sequence to masked file and
    $seqlength = length ($seq);         # calculate the length of the sequence

    print "// -> Parsing $acc [$id]\n" if ($verbose);
    
    # ERROR if we haven't added this accession to the database yet.
    # Use the accession for the mean time. 
    
    if ($id eq "") {
	print "// ERROR: No accession-id connection for this sequence\n" if ($verbose || $debug);
	$EST_name{$acc} = $acc;
	$id = $acc;
    }
    
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
	    print "// parse $feature\n" if ($verbose);
	    
	    $type  = $obj->Feature_data->Feature(1);         # Feature type (e.g. SL1,SL2,polyA)
	    $start = $obj->Feature_data->Feature(2);         # start coord
	    $stop  = $obj->Feature_data->Feature(3);         # stop coord
	    
	    $cut_to     = $start - 1;                        # manipulations for clipping 
	    $cut_from   = $stop;
	    $cut_length = $stop - $start + 1;
	    
	    if ($cut_to < 0 ) {$cut_to = 0;}                 # fudge to ensure non-negative clipping coords
	    
	    print "$acc [$id]: '$type' $start -> $stop [$cut_to : $cut_from ($cut_length)]\n" if ($debug);
	    print "// # $acc [$id] $type:" . (substr($seq,$cut_to,$cut_length)) . " [$start - $stop]\n\n";
       		$newseq = (substr($seqmasked,0,$cut_to))  . ('n' x $cut_length)  . (substr($seqmasked,$cut_from));
	    $seqmasked = $newseq;
	    
	}
    }
    
    print OUTPUT ">$acc $id\n$seqmasked\n";
    
    # close object
    $obj->DESTROY();
    
}
close INPUT;
$/ = "\n";


close OUTPUT;
close LOG;

exit(0);


###############################################################

sub recreate_hashes {
    my $hash = shift;
    
    open (FH, "<$hash") or die "Can't open file $hash : $!";
    undef $/;
    $data = <FH>;
    eval $data;
    die if $@;
    $/ = "\n";
    close FH;
}


###############################################################

sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name
  my $rundate     = `date +%y%m%d`; chomp $rundate;
  $log        = "/wormsrv2/logs/$script_name.$rundate.$$";

  open (LOG, ">$log") or die "cant open $log";
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

=item -est, process EST datafile 

=item -ost, process OST datafile 

=item -mrna, process mRNA datafile 

=item -debug <user>, Debug mode - log file only goes to named user

=item -verbose, Verbose mode toggle on extra command line output

=item -help, these help pages

=back

=cut
