#!/usr/bin/perl

# Calculate the overall %identity of ortholog pairs
# and identity distribution across the genome
#  Todd Harris (harris@cshl.org)
#  DATE : 05 May 2003


use strict;
use Statistics::Descriptive;
use Statistics::TTest;
use Statistics::PointEstimation;
use lib '/Users/todd/projects/briggsae/lib';
use Parse;

my $arms = {
	    I =>   [5000000, 10000000, 15080475],
	    II =>  [5300000, 12000000, 15174006],
	    III => [4000000, 10500000, 13855084],
	    IV  => [4500000, 12500000, 17493985],
	    V   => [6000000, 15000000, 20916337],
	    X   => [7000000, 11000000, 17749735]
	   };

my $CB_GENES = '19528';
my $CE_GENES = '18808';

my $file = shift;
chomp $file;
my $parse   = Parse->new();
my $elegans = $parse->elegans('chrom');
my $orthologs = $parse->orthologs($file);

# Fetch all elegans genes and store them by their
# chrom and position in the chromosome

my $PRINT_STATS = shift;
my $PRINT_TTEST;
$PRINT_TTEST = ($PRINT_STATS) ? undef : 1;


printf ("%-12s %-12s %-12s %-12s %-12s\n",'right','left','arms','cluster','whole') if ($PRINT_STATS);
print "orthologs_seen:windows_scored:mean (stdev/SE)\n" if ($PRINT_STATS);

my $totals = get_all_elegans_genes();
for my $chrom (qw/I II III IV V X/) {
#  check_for_orthologs($chrom,undef);
}
check_for_orthologs(undef,'all');
#check_for_orthologs(undef,'autosomes');


# =========
#   SUBS
# =========
sub get_all_elegans_genes {
  my $totals = {};
  # Dump out a file for easy R analysis
  open OUT,"elegans_genes.dist";
  print OUT join("\t",'chrom','gene','region'),"\n";
  foreach my $chrom (qw/I II III IV V X/) {

    # Get all of the genes in the left and right arms
    foreach my $ref (@{$elegans->{$chrom}}) {
      my $gene = $ref->{protein};
      my $start = $ref->{chrom_start};
      my $switch = get_switch($start,$chrom);
      push (@{$totals->{$chrom}->{$switch}->{all_genes}},$gene);
      print OUT JOIN("\t",$chrom,$gene,$switch),"\n";
    }
  }
  close OUT;
  return $totals;
}



sub check_for_orthologs {
  my ($current_chrom,$flag) = @_;
  my @chroms;
  if ($current_chrom) {
    push (@chroms,$current_chrom);
  } else {
    push (@chroms,qw/I II III IV V/);
  }
  push (@chroms,'X') if ($flag eq 'all');
  
  my $stats = { whole_chrom => Statistics::Descriptive::Full->new(),
		arms        => Statistics::Descriptive::Full->new(),
		cluster     => Statistics::Descriptive::Full->new(),
		right       => Statistics::Descriptive::Full->new(),
		left        => Statistics::Descriptive::Full->new()};

  # Iterate over chromosomes, regions, and window spans
  # storing the number of orthologs
  my $compiled = {};
  open OUT,">ortholog_hist.stats";
  foreach my $chrom (@chroms) {
    $compiled = stuff($totals->{$chrom},$compiled);
  }
  
  $stats->{whole_chrom}->add_data(@{$compiled->{whole_chrom}});
  $stats->{cluster}->add_data(@{$compiled->{cluster}});
  $stats->{right}->add_data(@{$compiled->{right}});
  $stats->{left}->add_data(@{$compiled->{left}});
  $stats->{arms}->add_data(@{$compiled->{arms}});

  my $label = ($current_chrom) ? $current_chrom : $flag;
  # Print out t-tests for this group, cluster vs arms...
  print_ttest(\@{$compiled->{arms}},\@{$compiled->{cluster}},$label);
  
  if ($PRINT_STATS) {
    print_stats($stats,$label);
  }
  
  print "\n" if ($PRINT_STATS);
}



# May have already started loading the totals hash
sub stuff {
  my ($data,$compiled) = @_;
  $compiled = {} if !$compiled;

  # Iterate over all of the genes in each region of this
  # chromosome in 100 gene windows...
  for my $region (qw/left right arms cluster whole_chrom/) {
    my @genes;
    if ($region eq 'arms') {
      push (@genes,@{$data->{left}->{all_genes}},@{$data->{right}->{all_genes}});
    } elsif ($region eq 'whole_chrom') {
      push (@genes,@{$data->{left}->{all_genes}},@{$data->{right}->{all_genes}},
	    @{$data->{cluster}->{all_genes}});
    } else {
      push @genes,@{$data->{$region}->{all_genes}};
    }


    for (my $i=0;$i<=scalar @genes/100;$i++) {
      # How many of the 100 genes in this window are also orthologs?
      my $total_orthologs;
      my @window = @genes[$i*100..($i*100)+99];
      foreach (@window) {
	$total_orthologs++ if ($orthologs->{$_});
      }

      push (@{$compiled->{$region}},$total_orthologs);
      # print OUT $region,"\t",$total_orthologs,"\n";
    }
  }
  return $compiled;
}


sub get_switch {
  my ($test,$chrom) = @_;
  my ($left,$right,$total) = @{$arms->{$chrom}};
  my $switch;
  $switch = 'left'  if ($test < $left);
  $switch = 'right' if ($test > $right);
  $switch = 'cluster' if ($test >= $left and $test <= $right);
  return $switch;
}


sub print_stats {
  my ($stats,$flag) = @_;
  printf ("%-10s",$flag);

  for my $key (qw/right left arms cluster whole_chrom/) {
    printf ("%-.0f:%-.0f:%-.2f (%-.2f/%-.2f) ",
	    $stats->{$key}->sum(),
	    $stats->{$key}->count(),
	    $stats->{$key}->mean(),
	    $stats->{$key}->standard_deviation(),
	    $stats->{$key}->variance());
  }
  print "\n";
}


sub print_ttest {
  my ($list1,$list2,$label) = @_;
  return unless ($PRINT_TTEST);
  my $ttest = Statistics::TTest->new();
  $ttest->set_significance(95);
  $ttest->load_data(\@$list1,\@$list2);
  print "\nT-Test output: $label\n";
  $ttest->output_t_test();
}
