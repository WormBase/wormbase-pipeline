#!/usr/bin/env perl
#
# RNASeq_determine_strandedness.pl
# 
# by Gary Williams                        
#
# Small script to investigate the strandedness of RNASeq libraries by
# looking at the region of the ama-1 gene, or homologs
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2015-06-04 11:25:17 $      

use strict;
use lib $ENV{'CVS_DIR'};
use Carp;
use Modules::RNASeq;
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species, $specified_experiment);
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
            "verbose"    => \$verbose,
            "store:s"    => \$store,
            "species:s"  => \$species, # the default is elegans
	    "srx:s"      => \$specified_experiment,
	   );


$species = 'elegans' unless $species;

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
			   );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


######################################
# variables and command-line options # 
######################################

my $count_no_bam=0;
my $status;
my $data;


my $RNASeq = RNASeq->new($wormbase, $log, 0, 0);

$log->write_to("Get experiments from config\n");
my $data = $RNASeq->get_transcribed_long_experiments();
my $data_count = scalar keys %{$data};
print "Have $data_count transcribed long experiments in total\n";

# only want to look at the experiments that we have not looked at before or a specified experiment
my $filtered;
if (defined $specified_experiment) {
  if (exists $data->{$specified_experiment}) {
    $filtered->{$specified_experiment} = $data->{$specified_experiment};
  } else {
    $log->log_and_die("$specified_experiment is an unknown experiment accession.\n");
  }
} else {
  foreach my $expt (keys %{$data}) {
    if (!exists $data->{$expt}{'strandedness'}) {
      $filtered->{$expt} = $data->{$expt};
    }
  }
}
$data = $filtered;

$data_count = scalar keys %{$data};
print "Have $data_count transcribed long experiments with no strandedness information\n";

my %regions = ( # chromosome, start, end and sense of highly expressed region of orthologs of ama-1
	       'brenneri'  => ['43', 207605, 224629, '-'                ], # CBN32550
	       'briggsae'  => ['chrIV', 15412678, 15431288, '-'         ], # WBGene00027817 Cbr-rbp-1
	       'brugia'    => ['Bm_v4_ChrX_scaffold_001', 19536809, 19549237, '-'], # Bma-ama-1
	       'elegans'   => ['IV', 4248070, 4258173, '+'              ], # WBGene00000123 ama-1
	       'japonica'  => ['17301', 25817, 31400, '-'               ], # Cjp-ama-1
	       'ovolvulus' => ['OVOC_OM1b', 20162608, 20175080, '-'     ], # WBGene00239543 OVOC2734
	       'pacificus' => ['59', 541834, 551440, '-'                ], # Ppa-ama-1
	       'remanei'   => ['125', 5842, 20968, '+'                  ], # Cre-ama-1
	       'sratti'    => ['SRAE_scaffold2', 418838, 423508, '-'    ], # sratti ama-1
	      );

foreach my $experiment_accession (keys %{$data}) {
  
  my $experiment = $data->{$experiment_accession};
  my $strandedness;
  my $config = $experiment->{strandedness};
  my $experiment_accession = $experiment->{experiment_accession};
  my $RNASeqSRADir = $RNASeq->{RNASeqSRADir};
  my $RNASeqGenomeDir = $RNASeq->{RNASeqGenomeDir};
  my $Software = $RNASeq->{Software};
  my $alignmentDir = $RNASeq->{'alignmentDir'};
  my $status;
  my $library_type=''; # for paired reads FR, FF, RR etc.
  
  if (!-e "$RNASeqSRADir/$experiment_accession/$alignmentDir/accepted_hits.bam") {$count_no_bam++; next}
  
  my $accession=$experiment->{'experiment_accession'};
  my $library_strategy=$experiment->{'library_strategy'};
  my $library_selection=$experiment->{'library_selection'};
  my $library_layout=$experiment->{'library_layout'};
  my $library_source=$experiment->{'library_source'};
  my $study_accession=$experiment->{'secondary_study_accession'};
  my $instrument_platform=$experiment->{'instrument_platform'};
  
  if (!defined $config) {
    $config = '';
  }
  
  chdir "$RNASeqSRADir/$experiment_accession/$alignmentDir";
  
  my $chrom         = $regions{$species}[0];
  my $start         = $regions{$species}[1];
  my $end           = $regions{$species}[2];
  my $region_sense  = $regions{$species}[3];
  
  # now parse the BAM file looking at our region of interest where the ama-1 gene is
  my $forward_count=0;
  my $reverse_count=0;
  my $forward_first=0;
  my $forward_second=0;
  my $reverse_first=0;
  my $reverse_second=0;
  
  # Flag values
  my $sequence_paired=0x0001;
  my $sequence_mapped=0x0004;
  my $mate_unmapped=0x0008;
  my $strand_of_query=0x0010; # 1 for reverse
  my $strand_of_mate=0x002;
  my $first_read_in_pair=0x0040;
  my $second_read_in_pair=0x0080;
  
  my $percentage=100;

  my $samtools_view = "$Software/samtools/samtools view accepted_hits.bam '${chrom}:${start}-${end}'"; # look at the alignments in the selected region
  open (HITS, "$samtools_view |") || $log->log_and_die("can't run $samtools_view in determine_strandedness()\n");
  while (my $line = <HITS>) {
    
    # SRR068555.4396580       0       IV      4256715 255     1S35M   *       0       0       CCCAACATCTCCACGCGGATTCTCGTCGCCACAGTA    BCB?ACCBCBCBABBBBA8ABBBBBABBB@;;9?9>    NH:i:1  HI:i:1  AS:i:34 nM:i:0
    
    my ($bam_id, $bam_flags, $bam_chrom, $bam_pos) = ($line =~ /^(\S+)\s+(\d+)\s+(\S+)\s+(\d+)/);
    if ($bam_flags & $strand_of_query) { # find reverse alignments
      $reverse_count++;
      if ($bam_flags & $first_read_in_pair) {
	$reverse_first++;
      } elsif ($bam_flags & $second_read_in_pair) {
	$reverse_second++;
      }
    } else {
      $forward_count++;
      if ($bam_flags & $first_read_in_pair) {
	$forward_first++;
      } elsif ($bam_flags & $second_read_in_pair) {
	$forward_second++;
      }
    }
  }
  
  # if the ama-1 ortholog is in the reverse strand, then swap over the counts to the other strand
  if ($region_sense eq '-') {
    ($reverse_count, $forward_count) = ($forward_count, $reverse_count);
    ($reverse_first, $forward_first) = ($forward_first, $reverse_first);
    ($reverse_second, $forward_second) = ($forward_second, $reverse_second);
  }
  
  # calculate the strandedness and store it in the config file
  my $study_accession = $experiment->{secondary_study_accession};
  my $experiment_ini = $RNASeq->read_experiments_from_study_config($study_accession);
  
  # single ended
  if ($library_layout eq 'SINGLE') {
    if ($forward_count > 10 * $reverse_count && $forward_count > 50) {
      $strandedness = 'stranded';
      $percentage = $reverse_count * 100 / $forward_count;

    } elsif ($reverse_count > 10 * $forward_count && $reverse_count > 50) {
      $strandedness = 'reverse_stranded'; # reads are stranded, but in the opposite orientation to the genes
      $percentage = $forward_count * 100 / $reverse_count;
      
    } elsif ($forward_count > 50 && $reverse_count > 50) {
      $strandedness = 'unstranded';
      
    } else {
      $strandedness = 'unknown';
      
    }
    
  } elsif ($library_layout eq 'PAIRED') {
    if ($forward_first > 10 * $reverse_first && $forward_first > 50) {
      $strandedness = 'stranded';
      $percentage = $reverse_first * 100 / $forward_first;

      if ($reverse_second > 10 * $forward_second && $reverse_second > 50) {
	$library_type = 'fr'; # first read mate is 'f'orward, second read mate is 'r'eversed
      } elsif ($forward_second > 10 * $reverse_second && $forward_second > 50) {
	$library_type = 'ff';
      } else {
	$library_type = 'unknown';
      }
      
    } elsif ($reverse_first > 10 * $forward_first && $reverse_first > 50) {
      $strandedness = 'reverse_stranded'; # reads are stranded, but in the opposite orientation to the genes
      $percentage = $forward_first * 100 / $reverse_first;

      if ($reverse_second > 10 * $forward_second && $reverse_second > 50) {
	$library_type = 'rr'; # this is rarely used, I think
      } elsif ($forward_second > 10 * $reverse_second && $forward_second > 50) {
	$library_type = 'rf'; # not sure this is a real thing, but look for it anyway
      } else {
	$library_type = 'unknown';
      }
      
    } elsif ($forward_first > 50 && $reverse_first > 50) {
      $strandedness = 'unstranded';
      $library_type = 'unknown';
      
    } else {
      $strandedness = 'unknown';
      $library_type = 'unknown';
    }
    
    $experiment_ini->newval($experiment_accession, 'library_type', $library_type);
    
  } else { # library_layout not specified
    $strandedness = 'unknown';
  }
  
  
  $log->write_to("Experiment=$experiment_accession $study_accession $library_source $library_strategy $library_selection $library_layout $instrument_platform strand=$strandedness library_type=$library_type\nForward: $forward_count Reverse: $reverse_count\nForward first: $forward_first Forward second: $forward_second Reverse first: $reverse_first Reverse second: $reverse_second Percentage: $percentage\n");
  
  $experiment_ini->newval($experiment_accession, 'strandedness', $strandedness);
  $experiment_ini->newval($experiment_accession, 'strandedness_percent', $percentage);
  $experiment_ini->RewriteConfig;
  
  
}

print "Have $count_no_bam experiments with no BAM files\n";



$log->mail();
print "Finished.\n" if ($verbose);
exit(0);
