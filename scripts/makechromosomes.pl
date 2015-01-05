#!/software/bin/perl -w
#
# makechromosomes.pl
#
# Populates the chromosome objects Subsequences based on Overlap_right tags
# pad
#
# Last updated by: $Author: pad $
# Last updated on: $Date: 2015-01-05 16:39:22 $
 
$!=1;
use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace ;


##############################
# command-line options       #
##############################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($stlace, $camace);
 
# Database name for databases other than ~wormpub/DATABASES/camace
my $db;         
# output file for ace
my $acefile;		      

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
            "verbose"    => \$verbose,
            "store:s"    => \$store,
            "db=s"       => \$db,
	    "acefile=s"  => \$acefile,
	    "stlace"    => \$stlace,  # set if we are dealing with stlace
	    "camace"    => \$camace,  # the default
            );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             );
}
 
# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$camace = 1; # the default
if ($stlace)  {
  $camace = 0;
  $stlace = 1;
}


##############################
# Hardcoded LH seed clones   #
##############################

my %clone2super;

%clone2super = ("cTel33B" => "I",
		"cTel52S" => "II",
		"cTel54X" => "III",
		"cTel4X"  => "IV",
		"cTel3X"  => "V",
		"cTel7X"  => "X",
	       ) ;

##############################
# open output file
##############################

if (! $acefile) {die "-acefile is not specified\n";}
open (ACE, ">$acefile") || die "Can't open file $acefile\n";

##############################
# get data from camace       #
##############################

my $tace = $wormbase->tace;

my $campath = $wormbase->database('camace');
($campath = $db) if ($db);
# check database path
unless (-e "$campath/database/database.map") {$log->log_and_die("Failed to connect to $campath\n")};
warn "Acessing database $campath\n" if ($verbose);

# connect to database
my $camdb   = Ace->connect(-path=>"$campath",
			   -program => $tace) || die "failed to connect to database\n";


# first get lots of info about Genome_sequence objects

my (%isGenomeSequence,%length,%right,%rightOffset,%left,%isLinkCandidate,%isExternal);
my (%currSource,%currStart,%currEnd);
my (%CDSSource,%CDSStart,%CDSEnd);
my (%PseudoSource,%PseudoStart,%PseudoEnd);
my (%TransposonSource,%TransposonStart,%TransposonEnd);
my (%TranscriptSource,%TranscriptStart,%TranscriptEnd);
my ($it,$obj,$seq,$start,$end);

my $error = 0;			# error status to return from program

$it = $camdb->fetch_many(Genome_sequence => '*') ;
while ($obj = $it->next) {
  if (! defined $obj) {
    $log->write_to("Genome_sequence not defined\n");
    $error = 1;
  }
    $isGenomeSequence{$obj} = 1 ;
    $length{$obj}           = $obj->DNA(2) ;
    if (!$length{$obj})     { 
      $log->write_to("No length for $obj\n") ; 
      $error = 1;
    }
    $right{$obj}            = $obj->Overlap_right ;
    $rightOffset{$obj}      = $obj->Overlap_right(2) ;
    $left{$obj}             = $obj->Overlap_left ;

    $isLinkCandidate{$obj}  = (!$isExternal{$obj} && $length{$obj} > 0);
}

# then some info about current links

$it = $camdb->fetch_many(Sequence => 'CHROMOSOME*') ;
while ($obj = $it->next) {
    foreach $a ($obj->at('Structure.Subsequence')) {
      if (! defined $a) {
	$log->write_to("superlink Structure.Subsequence not defined\n"); 
	$error = 1;
      }
      ($seq, $start, $end) = $a->row;
      if (! defined $seq || ! defined $start || ! defined $end) {
	$log->write_to("Structure.Subsequence row not defined\n"); 
	$error = 1;
      }
      $currSource{$seq}    = $obj;
      $currStart{$seq}     = $start;
      $currEnd{$seq}       = $end;
	
#	print "// push $seq to subsequence hash\n";

    }
  }

warn "Stored data to hash\n" if ($verbose);

###########################################
# make links
###########################################

my ($lk,%start,%end,%link,$parent,$startright);

foreach $seq (keys %isGenomeSequence) {
    
    # only keep seeds
    next if (!$isLinkCandidate{$seq} ||
	     ($left{$seq} && $isLinkCandidate{$left{$seq}} && $rightOffset{$left{$seq}}) ||
	     !$rightOffset{$seq} ||
	     !$isLinkCandidate{$right{$seq}});

    # print LINK header

      #$lk = "SUPERLINK_CB_$clone2super{$seq}";
      $lk = "CHROMOSOME_$clone2super{$seq}";
      print ACE "\nSequence $lk\n";
      print ACE "From_laboratory HX\n";


    # loop over subsequences
    $startright = 1;
    while ($isLinkCandidate{$seq}) {
	$start{$seq}  = $startright; 
	$end{$seq}    = $startright + $length{$seq} - 1;
	print ACE "Subsequence $seq $start{$seq} $end{$seq}\n";
	$link{$seq}   = $lk;
	if (!$rightOffset{$seq}) {
	    warn "ending loop here because rightOffset{$seq} is not set\n" if ($verbose);
	    last;
	}     # POSS EXIT FROM LOOP
	$startright   = $startright + $rightOffset{$seq} - 1;
	$seq          = $right{$seq};
    }		
}


###########################################
# Do a quick check to make sure that everything 
# that was in a link has been put back in one
###########################################

foreach $seq (keys %currSource) {
    if (!$link{$seq}) { 
      $log->write_to("$seq not put back into a link\n"); 
      $error = 1;
    }
}


$camdb->close;

close (ACE);


# Close log files and exit
$log->mail();
print "Finished.\n" if ($verbose);
exit($error);			# return the error status


############# end of file ################

__END__

=pod

=head2   NAME - makesuperlinks.pl

=head1 USAGE

=over 4

=item makesuperlinks.pl [-options]

=back

makesuperlinks queries an ACEDB database and generates the SUPERLINK
objects to an .acefile. This uses the Overlap_right tags within the
database and a hard-coded hash of starting clones within this
script

makesuperlinks mandatory arguments:

=over 4

=item B<-acefile>, file to output ACE to

=back

makesuperlinks OPTIONAL arguments:

=over 4

=item B<-db text>, database mode. Only dumps acefiles for the named database. The default is ~wormpub/DATABASES/camace

=item B<-debug>, send output to specified user only 

=item B<-verbose>, verbose report

=item B<-help>, this help page

=back

=head1 AUTHOR (& person to blame)

=over 4

=item Dan Lawson dl1@sanger.ac.uk

=back

=cut
