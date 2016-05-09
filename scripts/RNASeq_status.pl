#!/usr/bin/env perl
#
# RNASeq_status.pl
#
# Gives some handy stats on the current types of short-read libraries we know about.
# 
# by Gary Williams                        
#
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2014-11-12 15:36:23 $      

use strict;
use lib $ENV{'CVS_DIR'};
use Carp;
use Modules::RNASeq;

use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species, $all);
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
            "verbose"    => \$verbose,
            "store:s"    => \$store,
	    "species:s"  => \$species,
	    "all"        => \$all, # report not only the data used by the alignemnt pipeline, but also GENOMIC and everything else
);



$species = 'elegans' if !defined $species;

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

my %accessors = ($wormbase->species_accessors);
$accessors{$wormbase->species} = $wormbase;

foreach my $spec (sort { $accessors{$a}->full_name cmp $accessors{$b}->full_name } keys %accessors) {
  #$log->write_to("\nGenerating $spec RNASeq overview stats...\n");
  stats($spec, $accessors{$spec});
}   
  
$log->mail();
print "Finished.\n" if ($verbose);
exit(0);

######################################
# output the stats

sub stats {
  my ($spec, $wb) = @_;

  my $RNASeq = RNASeq->new($wb, $log);
  my $aligned = $RNASeq->get_transcribed_long_experiments();
  my @aligned_expts = keys %{$aligned};
  my $expts = $RNASeq->get_all_experiments();
  my %total_studies;
  my %categories;
  my %runs;
  my %studies;
  my %centers;
  my %instruments;
  my %layout;
  my $run_count=0;
  
  foreach my $experiment (keys %{$expts}) {
    if (! $all && !grep /$experiment/, @aligned_expts) {next} 
    my $study_accession = $expts->{$experiment}{study_accession};
    
    my $library_source = $expts->{$experiment}{library_source};
    my $library_strategy = $expts->{$experiment}{library_strategy};
    my $library_selection = $expts->{$experiment}{library_selection};
    
    my $library_layout = $expts->{$experiment}{library_layout};
    my $center_name  = $expts->{$experiment}{center_name};
    my $instrument_platform = $expts->{$experiment}{instrument_platform};
    
    $run_count += $expts->{$experiment}{run_count};
    
    
    $total_studies{$study_accession}=1;
    $categories{$library_source}{$library_strategy}{$library_selection}++; # experiment count for each category
    $runs{$library_source}{$library_strategy}{$library_selection} +=  $expts->{$experiment}{run_count}; # run counts for each category
    $studies{$library_source}{$library_strategy}{$library_selection}{$study_accession} = 1; # unique hash for each study in each category
    $centers{$center_name}++;
    $instruments{$instrument_platform}++;
    $layout{$library_layout}++;
    
  }
  
  $log->write_to("\n\n$spec Total RNASeq Experiment stats:\n\n");
  my $total_no_studies = keys %total_studies;
  $log->write_to("Studies:\t$total_no_studies\n");
  my $no_experiments = keys %{$expts};
  $log->write_to("Experiments:\t$no_experiments\n");
  $log->write_to("Runs:\t\t$run_count\n");
  my $paired_end = $layout{PAIRED};
  $log->write_to("Paired-end:\t$paired_end\n");
  my $single_end = $layout{SINGLE};
  $log->write_to("Single-end:\t$single_end\n");
  my @centers = (keys %centers);
  my $centers = join ', ', @centers;
  $log->write_to("Centres:\t$centers\n");
  
  my %trans_studies = ();
  my $trans_expts = 0;
  my $trans_runs = 0;
  $log->write_to("\nTRANSCRIPTOMIC Experiment stats:\n\n");
  foreach my $strategy (keys  %{$studies{'TRANSCRIPTOMIC'}}) {
    foreach my $selection (keys  %{$studies{'TRANSCRIPTOMIC'}{$strategy}}) {
      foreach my $study (keys  %{$studies{'TRANSCRIPTOMIC'}{$strategy}{$selection}}) {
	$trans_studies{$study} = 1;
      }
    }
  }
  foreach my $strategy (keys  %{$categories{'TRANSCRIPTOMIC'}}) {
    foreach my $selection (keys  %{$categories{'TRANSCRIPTOMIC'}{$strategy}}) {
      $trans_expts += $categories{'TRANSCRIPTOMIC'}{$strategy}{$selection};
      $trans_runs += $runs{'TRANSCRIPTOMIC'}{$strategy}{$selection};
    }
  }
  my $trans_no_studies = keys %trans_studies;
  $log->write_to("TRANSCRIPTOMIC Studies:\t$trans_no_studies\n");
  $log->write_to("TRANSCRIPTOMIC Expts:\t$trans_expts\n");
  $log->write_to("TRANSCRIPTOMIC Runs:\t$trans_runs\n");
  

  $log->write_to("\n\n");
  $log->write_to("Expt categories (some Studies will be counted more than once if they have Expts in more than one category)\n");
  my $string = sprintf("%-15s %-9s %-19s %-7s %-7s %-7s", "Source", "Strategy", "Selection", "Studies", "Expts", "Runs");
  $log->write_to("$string\n");
  $string = sprintf("%-15s %-9s %-19s %-7s %-7s %-7s", "------", "--------", "---------", "-------", "-----", "----");
  $log->write_to("$string\n");
  
  foreach my $source (keys %categories) {
    foreach my $strategy (keys %{$categories{$source}}) {
      foreach my $selection (keys %{$categories{$source}{$strategy}}) {
	
	my $studies = keys %{$studies{$source}{$strategy}{$selection}};
	my $no_expts = $categories{$source}{$strategy}{$selection};
	my $runs = $runs{$source}{$strategy}{$selection};
	$string = sprintf("%-15s %-9s %-18s %6s %6s %6s", $source, $strategy, $selection, $studies, $no_expts, $runs);
	$log->write_to("$string\n");
      }
    }
  }
}

