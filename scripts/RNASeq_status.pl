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
# Last updated on: $Date: 2014-10-09 12:41:45 $      

use strict;
use lib $ENV{'CVS_DIR'};
use Carp;
use Modules::RNASeq;

use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species);
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
            "verbose"    => \$verbose,
            "store:s"    => \$store,
);



$species = 'elegans';

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

  my $expts = $RNASeq->get_all_experiments();
  my %studies;
  my %categories;
  my %centers;
  my %instruments;
  my %layout;
  my $run_count=0;
  
  foreach my $experiment (keys %{$expts}) {
    my $study_accession = $expts->{$experiment}{study_accession};
    
    my $library_source = $expts->{$experiment}{library_source};
    my $library_strategy = $expts->{$experiment}{library_strategy};
    my $library_selection = $expts->{$experiment}{library_selection};
    
    my $library_layout = $expts->{$experiment}{library_layout};
    my $center_name  = $expts->{$experiment}{center_name};
    my $instrument_platform = $expts->{$experiment}{instrument_platform};
    
    $run_count += $expts->{$experiment}{run_count};
    
    
    $studies{$study_accession}=1;
    $categories{$library_source}{$library_strategy}{$library_selection}++;
    $centers{$center_name}++;
    $instruments{$instrument_platform}++;
    $layout{$library_layout}++;
    
  }
  
  $log->write_to("\n\n$spec RNASeq Experiment stats:\n\n");
  my $no_studies = keys %studies;
  $log->write_to("Studies:\t$no_studies\n");
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
  
  $log->write_to("\n\n");
  
  my $string = sprintf("%-15s %-9s %-19s %s", "Source", "Strategy", "Selection", "Count");
  $log->write_to("$string\n");
  $string = sprintf("%-15s %-9s %-19s %-5s", "------", "--------", "---------", "-----");
  $log->write_to("$string\n");
  
  foreach my $source (keys %categories) {
    foreach my $strategy (keys %{$categories{$source}}) {
      foreach my $selection (keys %{$categories{$source}{$strategy}}) {
	
	my $count = $categories{$source}{$strategy}{$selection};
	$string = sprintf("%-15s %-9s %-19s %-d", $source, $strategy, $selection, $count);
	$log->write_to("$string\n");
      }
    }
  }
}

