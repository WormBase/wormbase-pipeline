#!/usr/bin/env perl
#
# find_intergenic.pl
#
# Return intergenic region sequences, where "intergenic" is defined as maximal genomic regions containing no genes
#
# Note: due to gene clusters, and genes in the introns in other genes, the resulting number of regions will usually
# be significantly lower than the total gene count
#
# by Gary Williams
#
# Last updated by: $Author: klh $                      
# Last updated on: $Date: 2012-06-20 08:43:54 $        

use strict;                                      
use Getopt::Long;
use Carp;
use Storable;
use IO::Handle;
use Cwd;
use Ace;

use Bio::SeqIO;
use Bio::PrimarySeq;

use lib $ENV{'CVS_DIR'};
use Wormbase;
use Log_files;
use Sequence_extract;


######################################
# variables and command-line options #
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($output, $proximity, $side, $operons, $dbdir, $species);
#$output = "";       		# file to write output to
$proximity = 0;              # region around gene to restrict output to (<= 0 write complete intergenic region)
$side = "both";	        # either "both", "5", "3" side of gene to output
$operons = "include";	# either "include", "only", "no" intergenic regions inside operons

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "test"        => \$test,
	    "verbose"     => \$verbose,
	    "store:s"     => \$store,
	    "output=s"    => \$output,
	    "proximity=i" => \$proximity,
	    "side=s"      => \$side,
	    "operons=s"   => \$operons,
	    "database=s"  => \$dbdir,
	    "species:s"   => \$species,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species,
			     );
}

# Display help if required
&usage("Help") if ($help);

$dbdir          = $wormbase->test ? $wormbase->database("current") : $wormbase->autoace unless $dbdir; # Database path
my $gffdir      = $wormbase->gff_splits;        # GFF directory
my @chromosomes = $wormbase->get_chromosome_names; # chromosomes


# establish log file.
my $log = Log_files->make_build_log($wormbase);


# Check options are valid
if ($side ne "both" && $side ne "5" && $side ne "3") {
    die "Error: option -side=$side is invalid (requires 'both', '5' or '3')\n";
}

if ($operons ne "include" && $operons ne "only" && $operons ne "no") {
    die "Error: option -operons=$operons is invalid (requires 'include', 'only' or 'no')\n";
}
if (!defined($output)){
  $output = $wormbase->sequences."/intergenic_sequences.dna" unless $output;
}


########################################
# get transcripts out of the gff files #
########################################
   
my @line;			# GFF input line
my $gene_name;			# name of gene or operon
my $no_sequences = 0;		# count number of sequences written out


my $seq_obj = Sequence_extract->invoke($dbdir, undef, $wormbase);

my $seqio = Bio::SeqIO->new(-format => 'fasta',
                            -file   => ">$output");

foreach my $chromosome (@chromosomes) {

  # start with a fresh set of genes and operons for each chromosome
  my %gene = ();
  my %gene_count = ();

  $log->write_to("Processing chromosome $chromosome\n") if ($verbose);

  # loop through the CHROMOSOME GFF transcript file  
  # lines are like:
  # CHROMOSOME_X    gene    gene    1316    1935    .       +       .       Gene "WBGene00008351"

  print "Loop through Gene GFF file\n" if ($verbose);
  my $GFF = $wormbase->open_GFF_file("${\$wormbase->chromosome_prefix}${chromosome}", 'gene', $log);

  while (<$GFF>) {
    chomp;
    s/^\#.*//;
    next unless /\S/;
    @line = split /\t/;
    
    # we want genes and we may want operons
    if ($line[1] eq "gene") {

      ($gene_name) = ($line[8] =~ /Gene\s+\"(\S+)\"/);
      $gene_count{$gene_name}++;
      if ($gene_count{$gene_name} > 1) {
	print "WARNING: Gene name $gene_name is duplicated\n";
      }
      # NB store gene name with the start position, end position, strand, and operon status = 0
      $gene{$gene_name} = [$line[3], $line[4], $line[6], 0];
      $log->write_to("Gene : '$gene_name'\n") if ($verbose);
    } elsif ($line[1] eq "operon" && $operons ne "include") {
      # if $operons eq "include" we don't want to take any notice of operons, so don't read them in
      ($gene_name) = ($line[8] =~ /Operon\s+\"(\S+)\"/);
      # NB store operon name with the start position, end position and strand, and operon status = 1
      $gene{$gene_name} = [$line[3], $line[4], $line[6], 1];
      $log->write_to("Operon : '$gene_name'\n") if ($verbose);
    }
  }
  close ($GFF);

  ####################################################
  # Sort by start position then operons before genes #
  ####################################################

  my @sorted_genes = sort { 
    $gene{$a}->[0] <=> $gene{$b}->[0]
	or
    $gene{$b}->[3] <=> $gene{$a}->[3]
  } keys %gene;


  ########################
  # produce output files #
  ########################

  my $last_end = 0;	        # the position of the end of the previous gene
  my $operon_end = -1;		# end position of last operon found
  my $operon_start = -1;	# start position of last operon found
  # name of previous gene
  my $last_name = "startofsequence";
  my $last_strand = "";		# strand of the previous gene
  my $in_operon = 0;		# flag: true if in operon
  my $sequence;			# sequence of intergenic region
  my $width;			# amount of sequence to output
  my $seq_start;		# position to start writing sequence from
  my $print_start;		# human-readable start coordinate

  $log->write_to("Produce output file\n") if ($verbose);

  # debug count
  #my $count=0;

  foreach $gene_name (@sorted_genes) {

    # debug - only do first few
    #if ($count++ > 5) {last;}      
       
    my $start = $gene{$gene_name}->[0];
    my $end = $gene{$gene_name}->[1];
    my $strand = $gene{$gene_name}->[2];
    my $operon_status = $gene{$gene_name}->[3];

    # are we in an operon?
    if ($operon_status == 1) {
      $operon_start = $start;
      $operon_end = $end;
      $in_operon = 1;
    } elsif ($in_operon && $start > $operon_end) {
      $in_operon = 0;
    }

    $log->write_to("name: $gene_name chromosome: $chromosome start: $start end: $end strand: $strand\n") if ($verbose);

    # $operons eq "include" means we have no operon data (operons are totally ignored) and just output non-overlapping genes
    # $operons eq "no" means we simply treat operons like big genes and just output non-overlapping genes and operons
    # $operons eq "only" means we check to see if this is a gene within an operon before outputting
    if ((($operons eq "no" || $operons eq "include") && $start-1 > $last_end+1) ||
	($operons eq  "only" && $in_operon && $operon_status == 0 && $start > $operon_start && $start-1 > $last_end+1)) {
      $log->write_to("intergenic region: $last_end+1 - $start-1 in_operon: $in_operon\n") if ($verbose);
      if ($proximity <= 0) {
	# just print the intergenic region
	# get width of intergenic distance
	$width = $start-$last_end-1;
	$seq_start = $last_end;
	if ($seq_start eq "0") {$seq_start = "1";}
	if ($seq_start eq "1") {$print_start = $seq_start;}
	$print_start = $seq_start+1 unless ($seq_start eq "1"); # human-readable start coordinate

	$sequence = $seq_obj->Sub_sequence("${\$wormbase->chromosome_prefix}$chromosome", "$seq_start", "$width");

        my $seq = Bio::PrimarySeq->new(-id   => "${last_name}_${gene_name}",
                                       -desc => "$chromosome $print_start, len:$width", 
                                       -seq  => uc($sequence));
        $seqio->write_seq($seq);
	$no_sequences++;
      } else {

	#
	# get the regions around the start and end of the intergenic region separately
	# we look first at the last gene and then at this gene
	#

	# print the region near the last (first) gene
	if ($side eq "both" || ($last_strand eq "+" && $side eq "3") || ($last_strand eq "-" && $side eq "5")) {
	  # get width of intergenic distance
	  $width = $start-$last_end-1;
	  # use the smaller of $width or $proximity
	  $width = $width < $proximity ? $width : $proximity;
	  $seq_start = $last_end;
	  $print_start = $seq_start+1; # human-readable start coordinate
	  my $prime =  ($last_strand eq "+") ? "3" : "5";

	  $sequence = $seq_obj->Sub_sequence("${\$wormbase->chromosome_prefix}$chromosome", "$seq_start", "$width");
	  if ($last_strand eq "-") {
	    $sequence = $seq_obj->DNA_revcomp($sequence);
	  }

          my $seq = Bio::PrimarySeq->new(-id   => "$last_name.${prime}prime",
                                         -desc => "$chromosome $print_start, len:$width",
                                         -seq  => uc($sequence));
          
          $seqio->write_seq($seq);
	  $no_sequences++;

	}

	# print the region near this (second) gene
	if ($side eq "both" || ($strand eq "+" && $side eq "5") || ($strand eq "-" && $side eq "3")) {
	  # get width of intergenic distance
	  $width = $start-$last_end-1;
	  # use the smaller of $width or $proximity
	  $width = $width < $proximity ? $width : $proximity;
	  $seq_start = $start-$width-1;
	  $print_start = $seq_start+1; # human-readable start coordinate
	  my $prime =  ($strand eq "+") ? "5" : "3";
	  # output sequence from $width before the gene to the start
	  $sequence = $seq_obj->Sub_sequence("${\$wormbase->chromosome_prefix}$chromosome", "$seq_start", "$width");
	  if ($strand eq "-") {
	    $sequence = $seq_obj->DNA_revcomp($sequence);
	  }

          my $seq = Bio::PrimarySeq->new(-id   => "$gene_name.${prime}prime",
                                         -desc => "$chromosome $print_start, len:$width",
                                         -seq  => uc($sequence));
          
          $seqio->write_seq($seq);
	  $no_sequences++;

	}

      }
    }

    # get the next end to work from
    if ($operons eq "no" || $operons eq "include") {
      # don't want intergenic regions within operons, so treat operons like an overlapping gene
      if ($end > $last_end) {
	$last_end = $end;
	$last_name = $gene_name;
      }
    } elsif ($operons eq "only") {
      # only want to store end positions of genes
      if ($operon_status == 0 && $end > $last_end) {
	$last_end = $end;
	$last_name = $gene_name;
      }
    }
    # update the strand of the previous gene (not the previous operon - ignore these)
    if ($operon_status == 0) {
      $last_strand = $strand;
    }

  }

  # and don't forget the intergenic region after the last gene to the end of the chromosome
  if ($operons ne "only") {
    if ($proximity <= 0) {
      # just print the intergenic region
      # get width of intergenic distance
      $width = 50000;
      $seq_start = $last_end;
      $print_start = $seq_start + 1;

      # next if it is exactly at the end of the chromosome.
      next if $seq_start >= $seq_obj->{LENGTH}->{"${\$wormbase->chromosome_prefix}$chromosome"};

      $sequence = $seq_obj->Sub_sequence("${\$wormbase->chromosome_prefix}$chromosome", "$seq_start", "$width");

      if ($sequence) {
         my $seq = Bio::PrimarySeq->new(-id   => ">${last_name}_endofsequence",
                                         -desc => "$chromosome $print_start len:$width",
                                         -seq  => uc($sequence));
	$seqio->write_seq($seq);
	$no_sequences++;
      }
    } else {
      # print the region near the end of the last gene
      if ($side eq "both" || ($last_strand eq "+" && $side eq "3") || ($last_strand eq "-" && $side eq "5")) {
	# get width of intergenic distance
	$width = $proximity;
	$seq_start = $last_end;
	$print_start = $seq_start+1; # human-readable start coordinate
	my $prime =  ($last_strand eq "+") ? "3" : "5";

	# output sequence from end of gene to $width past the end
	$sequence = $seq_obj->Sub_sequence("${\$wormbase->chromosome_prefix}$chromosome", "$seq_start", "$width");
	if ($sequence) { # if the gene ends at the end of the chromosome there will not be a sequence after the gene
	  if ($last_strand eq "-") {
	    # get reverse-comp of sequence
	    $sequence = $seq_obj->DNA_revcomp($sequence);
	  }

          my $seq = Bio::PrimarySeq->new(-id   => ">$last_name.${prime}prime",
                                         -desc => "$chromosome $print_start",
                                         -seq  => uc($sequence));
          $seqio->write_seq($seq);
	  $no_sequences++;
	}
      }
    }
  }
}

$log->write_to("Finished Chromosome loop\n") if ($verbose);


####################################
# print some statistics to the log #
####################################

$log->write_to("\n\nFinished writing output file\n");
$log->write_to("Wrote $no_sequences sequences\n");
$log->write_to("gzipping $output\n");
$wormbase->run_command("gzip -f $output", $log);

$wormbase->check_file("$output.gz", $log,
		      minsize => 6182319,
);

$log->write_to("Finished.");
$log->mail();

exit(0);

##############################################################
#
# Subroutines
#
##############################################################

###############################
# Prints help and disappears  #
###############################

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

=head2 NAME - find_intergenic.pl

=head1 USAGE

=over 4

=item find_intergenic.pl [-options]

=back

find_intergenic.pl writes out intergenic sequences in Fasta format. By
default it will write out all sequences between gene transcripts that
do not overlap other genes on either strand. It can optionally be
restricted to write out only up to a specified length of region around
genes. It can be restricted to only write the region around the 5' or
around the 3' end of genes. If only the proximal region around the
genes are being printed out, then the sequence will be output in the
same sense as that gene. By default, if the whole intergenic region is
being printed, only the forward sense sequence is output. (What should
the sense of a sequence spanning the region between a forward and a
reverse sense gene be?).

find_intergenic.pl mandatory arguments:


=over 4

=item -output, Specifies the output file name

=back

find_intergenic.pl optional arguments:

=over 4

=item -proximity integer, Restricts the output to proximal regions.

=over 4

If the intergenic region is larger than the integer parameter, then
the output will be restricted to the length of that integer, if the
intergenic region is smaller than the integer parameter, then the
smaller value will be used. The default value is 0 which turns the
proximity restriction off. If the proximal intergenic regions of
adjacent genes overlap, then that overlapping section of the genome
will be output twice, once for each gene.

=back

=item -side "both" | "5" | "3", Output one side or other of the gene.

=over 4

Restricts the output to either just
the 5' side of the gene, or the 3' side, or (the default) writes both
sides of the gene. This only has effect if the '-proximity' option is
set to a non-zero size.

=back

=item -operons "include" | "only" | "no", Include operons or not.

=over 4

If the parameter is set to "include" (the default), then operons have
no effect on the output - they are ignored. If it is set to "only",
then only those genes that are within an operon will be considered, so
only intergenic regions of genes within operons will be output. If it
is set to "no", then intergenic regions within operons are ignored and
not output.

=back

=item -debug, Verbose/Debug mode

=item -test, Test mode, generate the acefile but do not upload them 

=item -help, Help pages

=back

=cut
