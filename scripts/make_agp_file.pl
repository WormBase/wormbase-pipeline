#!/usr/local/bin/perl5.8.0 -w
#
# make_agp_file.pl
#
# by Dan Lawson (dl1@sanger.ac.uk)
#
# Last edited by: $Author: klh $
# Last edited on: $Date: 2012-09-10 10:55:44 $

use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use File::Copy;
use File::Path;
use Storable;

my $mfetch = "mfetch";

#################################
# Command-line options          #
#################################

my $help;
my $debug;
my $test;      # uses test environment
my $quicktest; # same as $test but only runs one chromosome
my $verbose;
my $store;

GetOptions ("help"         => \$help,
	    "debug:s"      => \$debug,
	    "test"         => \$test,
	    "verbose"      => \$verbose,
            "quicktest"    => \$quicktest,
	    "store:s"      => \$store
	   );

($test = 1) if ($quicktest);
my $wormbase;
if( $store ) {
  $wormbase = retrieve( $store ) or croak("cant restore wormbase from $store\n");
}
else {
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


#################################################################################
# Main Loop                                                                     #
#################################################################################

foreach my $chromosome (@gff_files) {

  our %stop=();
  our %start=();
  our %span=();
  our %acc=();
  our %clone=();
  our %ver=();
  my @report="";

  my $i = 1;
  my ($start,$stop,$clone,$acc,$ver,$gap_span,$offset,$span,$gpseq,$gspan,$limit,$span2get,$unique_stop) = "";
  my ($last_stop,$last_start);
  $last_start =0;
  $last_stop = 2;

  my $file = "$outdir/CHROMOSOME_$chromosome.agp";

  if (not "-e $gff_dir/CHROMOSOME_${chromosome}_clone_acc.gff") {
    $log->log_and_die("The gff file $outdir/CHROMOSOME_${chromosome}_clone_acc.gff doesn't exist.\n");
  }
  
  # read data from gff file
  open (GFF, "<$gff_dir/CHROMOSOME_${chromosome}_clone_acc.gff") 
      or $log->log_and_die("cant open $gff_dir/CHROMOSOME_${chromosome}_clone_acc.gff");

  while (<GFF>) {   
    # parse from each line
    ($start, $stop) = (/(\d+)\s+(\d+)/); 
    ($clone, $acc, $ver) = (/Sequence \"(\S+)\" acc\=(\S+)\s+ver=(\d+)/);
    
    # catch gaps
    if ($last_stop < $start) {
      $gap_span = $start - $stop{$i-1};
      
      if ($gap_span == 1) {
	push (@report,"Putative butt-end join\n");
      }
      else {
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
    $ver{$i}   = $ver;
    
    $last_stop  = $stop;
    $last_start = $start;
    $i++;
  }
    
  close GFF;

  $start{$i} = $stop{$i-1} + 1;
  $limit = $i;
  
  # write report lines (redundant clones, butt-ends, gaps) to the log file
  foreach (@report) {
    $log->write_to("$_");
  }
  @report = "";      # clean the report up
  
  open (OUT, ">$file");
  for ($i=1; $i<$limit;$i++) {
    $span2get = $start{$i+1} - $start{$i};
    
    if ($clone{$i} eq "gap") {
      my $msg = sprintf("%3d %8s\t[%8d => %8d] [%8d] : Pad with %6d bp of '-'s}\n", 
                        $i, $clone{$i}, $start{$i}, $stop{$i}, $start{$i+1},$span2get);
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

##################
# Check the files
##################

foreach my $chromosome (@gff_files) {

  $wormbase->check_file("$outdir/CHROMOSOME_$chromosome.agp", $log,
                        minsize => 18000,
			maxsize => 33000,
                        lines => ["^${chromosome}\\s+\\d+\\s+\\d+\\s+\\d+\\s+F\\s+\\S+\\s+1\\s+\\d+\\s+\\+"],
                        );
}

$log->mail();
exit(0);


#################################################################################
### Subroutines                                                               ###
#################################################################################


sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
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
