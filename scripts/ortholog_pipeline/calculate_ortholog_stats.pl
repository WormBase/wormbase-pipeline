#!/usr/bin/perl

# Calculate a variety of statistics on ortholog pairs
#  Todd Harris (harris@cshl.org)
#  DATE : 05 May 2003

use strict;
use Statistics::Descriptive;
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

my $total_genes;

my $dist      = {};
my $by_method = {};
my %methods;

my $ortho_density = calc_density_per_window();

foreach my $chrom_iterator (qw/I II III IV V X/) {
  # Fetch the boundaries for this chrom
  my $this_chrom = {};
  foreach my $ce (keys %$orthologs) {
    # Fetch out some data for easier access and storage later.
    my $chrom = $orthologs->{$ce}->{chrom};
    my $start = $orthologs->{$ce}->{ce_start};
    my $eval  = $orthologs->{$ce}->{evalue};
    my $per_id= $orthologs->{$ce}->{per_id};
    my $meth  = $orthologs->{$ce}->{meth};
    my $conf  = $orthologs->{$ce}->{conf};

    # Whoops! This gene is not on this chrom.
    next if $chrom ne $chrom_iterator;

    # Calculate some overall stats for all the orthologs
    $methods{$meth}++;
    push (@{$by_method->{per_id}->{$meth}},$per_id);
    push (@{$by_method->{eval}->{$meth}},$eval);
    push (@{$by_method->{conf}->{$meth}},$conf);
    if ($meth =~ /seg/) {
      push (@{$by_method->{per_id}->{all_best_mutuals}},$per_id);
      push (@{$by_method->{eval}->{all_best_mutuals}},$eval);
      push (@{$by_method->{conf}->{all_best_mutuals}},$conf);
    } elsif ($meth =~ /kb/) {
      push (@{$by_method->{per_id}->{all_kb_synteny}},$per_id);
      push (@{$by_method->{eval}->{all_kb_synteny}},$eval);
      push (@{$by_method->{conf}->{all_kb_synteny}},$conf);
    } else {
      push (@{$by_method->{per_id}->{all_synteny}},$per_id);
      push (@{$by_method->{eval}->{all_synteny}},$eval);
      push (@{$by_method->{conf}->{all_synteny}},$conf);
    }
    
    push (@{$by_method->{per_id}->{all}},$per_id);
    push (@{$by_method->{eval}->{all}},$eval);
    push (@{$by_method->{conf}->{all}},$conf);
    
    # Fetch out the start and stop positions so I cna place it
    my $switch = get_switch($start,$chrom);
    

    # This gives me the eval and per_id by region...
    push (@{$dist->{$chrom}->{$switch}->{per_id}},$per_id);
    push (@{$dist->{$chrom}->{$switch}->{eval}},$eval);
  }
}


# Print out overall statistics for percent identity by method
print "Ortholog percent identity by assignment method\n";
printf ("%-30s %-15s %-10s %-10s %-10s %-10s %-10s %10s %-10s %-10s %-10s %-10s %-10s\n",
	'method','total_genes','avg_per_id','SD','SE',
	#'eval','SD','SE',
	'conf','SD','SE',
	'%CB','%CE',
       );
$methods{all_best_mutuals}++;
$methods{all_synteny}++;
$methods{all_kb_synteny}++;
$methods{all}++;

my $genes_seen;
foreach my $method (sort keys %methods) {

  my $total = scalar @{$by_method->{per_id}->{$method}};
  $genes_seen += $total;
  my $id_stats = Statistics::Descriptive::Full->new();
  $id_stats->add_data(@{$by_method->{per_id}->{$method}});
  
  my $eval_stats = Statistics::Descriptive::Full->new();
  $eval_stats->add_data(@{$by_method->{eval}->{$method}});
  
  my $conf_stats = Statistics::Descriptive::Full->new();
  $conf_stats->add_data(@{$by_method->{conf}->{$method}});
		  
#  printf ("%-30s %-15s %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f %-8.2f\n",
  printf ("%-30s %-15s %-8.3f %-8.3f %-8.3f %-8.2f %-8.2f %-8.2f %-.2f %-.2f\n",
	  $method,$total,
	  $id_stats->mean(),$id_stats->standard_deviation(),$id_stats->variance(),
#	  $eval_stats->mean(),$eval_stats->standard_deviation(),$eval_stats->variance(),
	  $conf_stats->mean(),$conf_stats->standard_deviation(),$conf_stats->variance(),
	  ($total/$CB_GENES) * 100,($total/$CE_GENES) * 100
	 );
}
print "Total genes seen: $genes_seen\n";



sub calc_density_per_window {
  my $avg_density = {};
  
  # Get out 100 genes at a time for each segment of the chromosome
  # Do all chromosomes in one fell swoop
  open OUT,">dumped_stats/ortholog_distribution.stats";
  foreach my $chrom (keys %$elegans) {
    my ($left,$right,$total) = @{$arms->{$chrom}};
    my $switch;
    for (my $i=0;$i<=scalar @{$elegans->{$chrom}};$i+=99) {
      my @genes = @{$elegans->{$chrom}}[$i..$i+99];
      # Get the last gene in the list.
      # Where does it fall?
      my $end = $genes[-1]->{chrom_start};
      $switch = get_switch($end,$chrom);
      
      my $total_orthos;
      my $total_genes;
      my $count;
      foreach my $ref (@genes) {
	my $gene = $ref->{protein};
	# my $gene = $ref->{gene};
	my $start = $ref->{chrom_start};
	my $new_switch = get_switch($start,$chrom);
	
	# Have I moved to a new region?
	# If so, reset $i to start counting from there.
	if ($new_switch ne $switch) {
	  $i=$i+$count;
	  print OUT join("\t",$chrom,$switch,$total_orthos),"\n";
	  next;
	}
	
	# How many of the 100 genes in this window are also orthologs?
	$total_orthos++ if ($orthologs->{$gene} ne '');
	
	# Keep track of the total orthologs and total genes for this window
	$total_genes++;
	$count++;
      }
      print OUT join("\t",$chrom,$switch,$total_orthos),"\n";
    }   # Been through all the windows in these region...
  }
}





sub calc_stats {
  my $all_chroms = {};
  my $autosomes  = {};
  for my $chrom (qw/I II III IV V X/) {
    
    # Calculate the stats for each segment on this chromosome
    for my $switch (qw/left right cluster arms/) {
      my $stats   = Statistics::Descriptive::Full->new();
      
      # Collect left and right together for the arms
      my @data;
      if ($switch eq 'arms') {
	push (@data,@{$dist->{$chrom}->{left}->{windows}},
	      @{$dist->{$chrom}->{right}->{windows}});
      } else {
	push (@data,@{$dist->{$chrom}->{$switch}->{windows}});
	$stats->add_data(@data);
      }
    }
  }
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




# Now print out the stats for this data grouping
sub print_data {
  my ($data,$label) = @_;
  printf ("%-5s %-12s %-12s %-12s",'chrom','region','total_orthos','total_genes');
  printf ("%-25s %-25s\n",'avg_orthos/100_genes (SD/SE)','avg_per_id (SD/SE)');
  for my $switch (qw/left right cluster arms/) {
    # Fetch out the ortholog density for this region
    my $total_genes  = $ortho_density->{$label}->{$switch}->{total_genes};
    my $total_orthos = scalar @{$data->{$switch}->{per_id}};
    my $density_mean = sprintf ("%.2f",$ortho_density->{$label}->{$switch}->{mean});
    my $density_sd   = sprintf ("%.2f",$ortho_density->{$label}->{$switch}->{sd});
    my $density_var  = sprintf ("%.2f",$ortho_density->{$label}->{$switch}->{variance});
    
    # Fetch and calculate all the statistical values...
    #  for my $class (qw/per_id evalue/) {
    printf ("%-5s %-12s %-12s %-12s",$label,$switch,$total_orthos,$total_genes);
    printf ("%-25s",$density_mean . " ($density_sd/$density_var)");
    for my $class (qw/per_id/) {
      my $stats = Statistics::Descriptive::Full->new();
      $stats->add_data(@{$data->{$switch}->{$class}});
      my $mean  = sprintf ("%.2f",$stats->mean());
      my $stdev = sprintf ("%.2f",$stats->standard_deviation());
      my $var   = sprintf ("%.2f",$stats->variance());
      printf ("%-25s",$mean . " ($stdev/$var)");
      print "\t";
    }
    print "\n";
  }
}
