#!/usr/local/bin/perl5.6.1 -w
#
# make_agp_file.pl
# v0.3
#
# dl
#
# Usage : make_agp_file.pl
#
# Last edited by: $Author: krb $
# Last edited on: $Date: 2002-09-03 15:50:41 $


#################################################################################
# Initialise variables                                                          #
#################################################################################


$|=1;
use strict;
use vars qw ($debug $seq_len $sv_acc $sv_ver $opt_d $opt_h);
use Getopt::Std;
use lib '/wormsrv2/scripts';
use Wormbase;

 ##############################
 # Script variables (run)     #
 ##############################

my $maintainers = "All";
my $rundate     = `date +%y%m%d`;   chomp $rundate;
my $runtime     = `date +%H:%M:%S`; chomp $runtime;

my @gff_files = ('I','II','III','IV','V','X');

my $outdir    = "/wormsrv2/autoace/yellow_brick_road";

#my $datadir   = "/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS";
my $datadir    = "/wormsrv2/autoace/yellow_brick_road";

 ##############################
 # command-line options       #
 ##############################

getopts ('hd');
&error(0) if ($opt_h);

$debug = 1;
($maintainers = "dl1\@sanger.ac.uk") if ($debug);


 ##############################
 # getz query to build hashes #
 ##############################

# ID   AF016448   standard; DNA; INV; 44427 BP.
# SV   AF016448.1

our %seqver = ();
our %seqlen = ();

open (SEQUENCES, "getz -f \'id seqversion\' \"([emblnew-org:Caenorhabditis elegans] \& [emblnew-div:INV]) | ([embl-org:Caenorhabditis elegans] \& [embl-div:INV])\" |");
while (<SEQUENCES>) {
#    print "$_";
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

#my $last_start = 0;

foreach my $chromosome (@gff_files) {

  our %stop=();
  our %start=();
  our %span=();
  our %acc=();
  our %clone=();
  our %ver=();
  my @report="";


  my $i = 1;
  my ($start,$stop,$clone,$acc,$gap_span,$offset,$span,$gpseq,$gspan,$last_stop,$last_start,$limit,$span2get,$unique_stop) = "";

  my $file = "$outdir/CHROMOSOME_$chromosome.agp";
  $last_stop = 2;

  &error(1,$chromosome) unless ("-e $datadir/CHROMOSOME_${chromosome}.clone_acc.gff");
  
  # read data from gff file
  open (GFF, "<$datadir/CHROMOSOME_$chromosome.clone_acc.gff");
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

    if($seqver{$acc} eq ""){
      my $getz_acc = `/usr/local/pubseq/bin/pfetch $getz_id | grep ">"`;
      $getz_acc =~ s/.*\.//;
      $getz_acc =~ s/ Caenorh*//;
      $ver{$1} = $getz_acc;
    }

    
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

  open (LOG, ">${file}.log");
  
  # write report lines (redundant clones, butt-ends, gaps) to the log file
  foreach (@report) {
    print LOG $_;
  }
  @report = "";      # clean the report up
  
  open (OUT, ">$file");
  for ($i=1; $i<$limit;$i++) {
    $span2get = $start{$i+1} - $start{$i};
    
    if ($clone{$i} eq "gap") {
      printf LOG "%3d %8s\t[%8d => %8d] [%8d] : Pad with %6d bp of '-'s}\n", $i, $clone{$i}, $start{$i}, $stop{$i}, $start{$i+1},$span2get;
      $unique_stop = $start{$i} + $span2get -1;
      print OUT "$chromosome\t$start{$i}\t$unique_stop\t$i\tN\t$span2get\n";
    } 
    else {
      printf LOG "%3d %8s\t[%8d => %8d] [%8d] : Extract %6d bp from accession $acc{$i}, version $ver{$i}\n", $i, $clone{$i}, $start{$i}, $stop{$i}, $start{$i+1},$span2get;
      $unique_stop = $start{$i} + $span2get -1;
      print OUT "$chromosome\t$start{$i}\t$unique_stop\t$i\tF\t$acc{$i}.$ver{$i}\t1\t$span2get\t+\n";
    }
  }
  close LOG;
  close OUT;
  
  # copy agp file to correct directory
  system ("cp $file /wormsrv2/autoace/CHROMOSOMES/");
    
}

 ##############################
 # hasta luego                #
 ##############################

exit(0);


#################################################################################
### Subroutines                                                               ###
#################################################################################


sub error {
    my $error = shift;
    my $chromosome = shift;
    if ($error == 1){ 
	# No gff file to work from
	print "The gff file '$datadir/CHROMOSOME_${chromosome}.clone_acc.gff' doesn't exist.\n";
	exit(0);
    }
    elsif ($error == 0) {
        # Normal help menu
	exec ('perldoc',$0);
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

=back

=cut
