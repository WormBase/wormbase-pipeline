#!/usr/local/bin/perl5.8.0 -w
#
# make_agp_file.pl
#
# by Dan Lawson (dl1@sanger.ac.uk)
#
# Last edited by: $Author: gw3 $
# Last edited on: $Date: 2005-12-20 13:52:23 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use File::Copy;
use File::Path;
use Storable;


#################################
# Command-line options          #
#################################

my $help;
my $debug;
my $test;      # uses test environment
my $quicktest; # same as $test but only runs one chromosome
my $database;
my $verbose;
my $store;

GetOptions ("help"         => \$help,
	    "debug:s"      => \$debug,
	    "test"         => \$test,
	    "verbose"      => \$verbose,
            "quicktest"    => \$quicktest,
	    "database:s"   => \$database,
	    "store:s"      => \$store
	   );

# check that -test and -quicktest haven't both been set.  Also if -quicktest is specified, 
# need to make -test true, so that test mode runs for those steps where -quicktest is meaningless
if($test && $quicktest){
  die "both -test and -quicktest specified, only one of these is needed\n";
}
($test = 1) if ($quicktest);

my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
  ($test = 1) if ($quicktest);
  $wormbase->set_test($test);   # set test in the wormbase object

}
else {
  ($test = 1) if ($quicktest);
  $wormbase = Wormbase->new( -debug   => $debug,
			     -test    => $test,
			   );
}

# Display help if required
&usage("Help") if ($help);


my $log = Log_files->make_build_log($wormbase);

##############################
# Script variables (run)     #
##############################

# Set up top level base directory which is different if in test mode
# Make all other directories relative to this
my $basedir  = $wormbase->basedir;
my $outdir    = $wormbase->autoace."/yellow_brick_road";
mkpath $outdir unless (-e $outdir);
my $gff_dir = $wormbase->gff_splits;


my @gff_files = ('I','II','III','IV','V','X');
@gff_files = ('III') if ($quicktest);


##############################
# getz query to build hashes #
 ##############################

# ID   AF016448   standard; DNA; INV; 44427 BP.
# SV   AF016448.1

our %seqver = ();
our %seqlen = ();

open (SEQUENCES, "getz -f \'id seqversion\' \"([emblnew-org:Caenorhabditis elegans] \& [emblnew-div:INV]) | ([embl-org:Caenorhabditis elegans] \& [embl-div:INV]) | ([emblrelease-org:Caenorhabditis elegans] \& [emblrelease-div:INV])\" |");
while (<SEQUENCES>) {
    if (/^ID\s+(\S+)\s+standard\; DNA\; INV\; (\d+)/) {
	$seqlen{$1} = $2;
    }
    if (/^SV\s+(\S+)\.(\d+)/) {
      $seqver{$1} = $2;
    }
}
close(SEQUENCES);

print "Finished assigning to hash\n";

#################################################################################
# Main Loop                                                                     #
#################################################################################

my $seq_len;
my $sv_acc;
my $sv_ver;

foreach my $chromosome (@gff_files) {

  our %stop=();
  our %start=();
  our %span=();
  our %acc=();
  our %clone=();
  our %ver=();
  my @report="";

  my $i = 1;
  my ($start,$stop,$clone,$acc,$gap_span,$offset,$span,$gpseq,$gspan,$limit,$span2get,$unique_stop) = "";
  my ($last_stop,$last_start);
  $last_start =0;
  $last_stop = 2;

  my $file = "$outdir/CHROMOSOME_$chromosome.agp";

  &error(1,$chromosome) unless ("-e $gff_dir/CHROMOSOME_${chromosome}.clone_acc.gff");
  
  # read data from gff file
  open (GFF, "<$gff_dir/CHROMOSOME_$chromosome.clone_acc.gff") or $log->log_and_die("cant open $gff_dir/CHROMOSOME_$chromosome.clone_acc.gff");
  while (<GFF>) {
    
    $seq_len = "";
    $sv_acc  = "";
    $sv_ver  = "";
    
    # parse from each line
    ($start, $stop) = (/(\d+)\s+(\d+)/); 
    ($clone, $acc) = (/Sequence \"(\S+)\" acc\=(\S+)/);
    
    # catch gaps
    if ($last_stop < $start) {
      $gap_span = $start - $stop{$i-1};
      
      if ($gap_span == 1) {
	push (@report,"Putative butt-end join\n");
      }
      else {
	print "[ " . ($stop{$i-1}+1) . " : " . ($start-1) ."] so insert padding characters over the gap\n" if ($debug);
	push (@report, "$chromosome\t$stop{$i-1}\t" . ($start - 1) . "\t$i\tN\t$gap_span\n");
	$start{$i}  = $stop{$i-1} + 1;
	$stop{$i}   = $start-1;
	$span{$i}   = $gap_span -1;
	$acc{$i}    = "gap";
	$clone{$i}  = "gap";
	$last_start = $start;
	$last_stop  = $stop;
	$i++;
      }
    }
    
    # subsumed under previous clone
    if (($start > $last_start) && ($stop < $last_stop)) {
      push (@report,"$clone is redundant\n");
      next;
    }
    
    $clone{$i} = $clone;
    $start{$i} = $start;
    $stop{$i}  = $stop;
    $acc{$i}   = $acc;
    $span{$i}  = $stop - $start + 1;
    
#	printf "%8s [%8d => %8d : %6d] $acc{$i}.$sv_ver\n", $clone{$i}, $start{$i}, $stop{$i}, $span{$i};

    # fudge to try to get sequence version using pfetch if getz fails
    if($seqver{$acc} eq ""){
      my $getz_acc = `/usr/local/pubseq/bin/pfetch $acc | grep ">"`;
      $getz_acc =~ s/.*\.(\d+)\s+.*/$1/;
      $ver{$1} = $getz_acc;
    }
    else{
      $ver{$i}    = $seqver{$acc};
    }

    $last_stop  = $stop;
    $last_start = $start;
    $i++;
    
  }
    
  close GFF;
    
  $start{$i} = $stop{$i-1} + 1;
  $limit = $i;
  
  print "\n";
  
  # write report lines (redundant clones, butt-ends, gaps) to the log file
  foreach (@report) {
    $log->write_to("$_");
  }
  @report = "";      # clean the report up
  
  open (OUT, ">$file");
  for ($i=1; $i<$limit;$i++) {
    $span2get = $start{$i+1} - $start{$i};
    
    if ($clone{$i} eq "gap") {
      my $msg = sprintf "%3d %8s\t[%8d => %8d] [%8d] : Pad with %6d bp of '-'s}\n", $i, $clone{$i}, $start{$i}, $stop{$i}, $start{$i+1},$span2get;
      $log->write_to($msg);

      $unique_stop = $start{$i} + $span2get -1;
      print OUT "$chromosome\t$start{$i}\t$unique_stop\t$i\tN\t$span2get\n";
    } 
    else {
      my $msg = sprintf "%3d %8s\t[%8d => %8d] [%8d] : Extract %6d bp from accession $acc{$i}, version $ver{$i}\n", $i, $clone{$i}, $start{$i}, $stop{$i}, $start{$i+1},$span2get;
      $log->write_to($msg);
      $unique_stop = $start{$i} + $span2get -1;
      print OUT "$chromosome\t$start{$i}\t$unique_stop\t$i\tF\t$acc{$i}.$ver{$i}\t1\t$span2get\t+\n";
    }
  }
  close OUT;

  # copy agp file to correct directory
  # copy command returns 0 for failed
  my $status = copy("$file", $wormbase->chromosomes."/CHROMOSOME_$chromosome.agp"); 
  $log->write_to("ERROR: Couldn't copy file: $!\n") if ($status == 0);
}


$log->mail();
print "Finished.\n" if ($verbose);
exit(0);
   
 ##############################
 # hasta luego                #
 ##############################

exit(0);


#################################################################################
### Subroutines                                                               ###
#################################################################################



##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################

sub error {
  my $error = shift;
  my $chromosome = shift;
  if ($error == 1){ 
    # No gff file to work from
    print "The gff file '$outdir/CHROMOSOME_${chromosome}.clone_acc.gff' doesn't exist.\n";
    exit(0);
  }
}


__END__

=pod

=head2   NAME - make_agp_file.pl

=head1 USAGE

=over 4

=item make_agp_file.pl [-options]

=back

make_agp_file produces agp format lists for each chromosome. These files
delineate the yellow brick road (aka golden path) through a tiling set
of EMBL/GenBank entries (with sequence versions).

make_agp_file requires gff files with accessions (this is made after running
GFFsplitter late in the build procedure).

make_agp_file.pl mandatory arguments:

=over 4

=item none,

=back

make_agp_file.pl  OPTIONAL arguments:

=over 4

=item -test

Uses test environment in ~wormpub/TEST_BUILD/

=back

=over 4

=item -quicktest

Will only run against one chromosome (for speed) which is CHROMOSOME_III

=back


=back

=cut
