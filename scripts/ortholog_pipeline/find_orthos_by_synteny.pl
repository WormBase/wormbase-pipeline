#!/usr/bin/perl

# SHOULD REALLY save all the evalues and such when parsing the orthologs
# so that I can easily calculate ranges and distributions and such

=pod

=head1 NAME

  find_orthos_by_synteny

=head1 DESCRIPTION

  Assigns ortholog pairs to previously unassigned genes based on synteny.
  Ranges are bound by previously defined orthologs, or by user supplied
  ranges that will center on one specific ortholog
  
=head1 SYNOPSIS

   find_orthos_by_synteny
       --max_eval   The maximum e-value for significance of blast hits 
       --cutoff     Optional. The cutoff difference between high scoring pairs
       --method     A method tag for the output.  Purely descriptive.
       --range      The range in bp to examine.
       --orthologs  Optional. Path to predetermined file of orthologs. Genes in this
                    list will ee excluded from consideration during this analysis.
       --help       Display this documentation 

   The cutoff switch controls the threshold between two high-scoring
   blast hits.  The difference between the expectation value of the top
   hit and the second best hit must be greater than the --cutoff value
   in order to be assigned to an ortholog pair. The default value is 10^5.

   If range is provided, the script will center windows the size of
   --range on each ortholog pair, then examine those windows for
   potential orthologs.

=head1 AUTHOR

  Todd Harris (harris@cshl.org) 
  DATE : 05 May 2003 
  VERSION : $Id: find_orthos_by_synteny.pl,v 1.1 2004-05-21 10:33:38 ar2 Exp $

=head1 TODO

  # 1. check for hits in both directions (not checking reciprocals right now)
  # 2. check for possible paralogs outside of the interval

  Range based queries do not behave properly at the ends of supercontigs or chromosomes.
  This isn't a huge deal.  It really only affects the statistics.

  Should really check spans where there are one to many relationships.  Might be able
  to assign ortholog to these quite easily.

=head1 IMPORTANT NOTES

 Confidence values for synteny assigned ortholog pairs are not entirely
 accurate since they do not consider all blast hits, but only those within
 a specific region.  Furthermore, they only compare the ratio of evalues in
 a single direction.  USE THESE VALUES WITH SOME CAUTION.

=cut

#'

use Getopt::Long;
use Pod::Usage;
use Parse;
use Statistics::Descriptive;
use strict;

my ($max_eval,$cutoff,$method,$range,$orthos_file,$help);
GetOptions('max_eval=s' => \$max_eval,
	   'cutoff=s'   => \$cutoff,
	   'method=s'   => \$method,
	   'range=i'    => \$range,
	   'orthologs=s'=> \$orthos_file,
	   'help'       => \$help);
pod2usage(-verbose=>2) if $help;
pod2usage(0) unless ($method);

my $cutoff   ||= 100000;
my $max_eval ||= 1e-10;

# Additional globals
my $DEBUG   = 0;
my $ALERT   = 0;
my $SUCCESS = 1;
my $stats         = {};
my $new_orthologs = {};
my $parse = Parse->new();

open ERROR,">assignment_progress.$method";

# Fetch and order all elegans genes and relevant information
my $elegans = $parse->elegans('protein');
my @ordered_elegans = sort { $elegans->{$a}->{genome_start} <=> $elegans->{$b}->{genome_start} } keys %$elegans;

# Fetch and then order all orthologs according to their elegans gene
my $orthologs = $parse->orthologs_temp($orthos_file);
my $ordered_orthologs = order_orthologs(\@ordered_elegans);

# Fetch all of the briggsae genes on a given supercontig
# or alternatively a seperate hash keyed by gene
my $briggsae = $parse->briggsae('super');
my $briggsae_by_gene = $parse->briggsae('gene');

# Load blast lookup tables in both directions
my $elegans_blast  = $parse->blast('elegans');
my $briggsae_blast = $parse->blast('briggsae');

if ($range) {
  range_query();
} else {
  ortholog_span();
}

dump_new_orthologs();
print_stats();




# +++++++++++++++++++++
#   Begin subroutines
# +++++++++++++++++++++
sub ortholog_span {
  my $position = 0;
  foreach my $ce_left (@$ordered_orthologs) {
    my $ce_right = $ordered_orthologs->[$position+1];
    $position++;
    
    # Are we jumping to the next chromosome?
    next if ($elegans->{$ce_left}->{chrom} ne $elegans->{$ce_right}->{chrom});
    
    # 1. Fetch all of the elegans genes that occur between these two...
    #    Do this by fetching the individual indices of each ortholog pair
    my $ce_left_limit  = get_index($ce_left);
    my $ce_right_limit = get_index($ce_right);
    $stats->{spans_examined}++;
    
    print_alert("------> Evaluating new range: $ce_left $ce_right");
    my ($ce_range_start,$ce_range_stop,$ce_length) = get_range($elegans->{$ce_left}->{genome_start},
							       $elegans->{$ce_left}->{genome_stop},
							       $elegans->{$ce_right}->{genome_start},
							       $elegans->{$ce_right}->{genome_stop});
    
    # 2. Now, fetch out all elegans genes between these limits
    my $length = $ce_right_limit - $ce_left_limit;
    my @copy = @ordered_elegans;
    my @ce_candidates = splice(@copy,$ce_left_limit+1,$length-1);
    
    # 3. Look up the briggsae mates of the left limit and the right limit
    my $cb_left  = mate($ce_left);
    my $cb_right = mate($ce_right);
    
    # This is all just debugging code...
    print_debug("bounds: $ce_left $ce_right ($ce_left_limit-$ce_right_limit/$ce_range_start-$ce_range_stop)"
		. " $cb_left $cb_right");
    
    foreach (@ce_candidates) {
      print_debug("Ce genes in range: $_\t" . $elegans->{$_}->{genome_start} );
    }
    
    # 4. Fetch all briggsae genes on this supercontig that fall within this range.
    my ($cb_candidates,$cb_length,$flag) = get_cb_candidates($cb_left,$cb_right);
    if ($flag) { # I've seen a break in collinearity
      $stats->{collinearity_breaks}++;
      next;
    }
    
    my $total_ce = scalar @ce_candidates;
    my $total_cb = scalar @$cb_candidates;
    
    # Keep track of some general span statistics
    push (@{$stats->{spans}},[$cb_length,$ce_length,$total_cb,$total_ce]);
    $stats->{collinear_spans}++;
    
    # Are these two briggsae genes on the same supercontig?
    # It would be more efficient to place this above, but then I'd love the stats
    print_debug('Supercontig check: ' . $briggsae_by_gene->{$cb_left}->{supercontig} . ' ' . 
		$briggsae_by_gene->{$cb_right}->{supercontig});
    next unless ($briggsae_by_gene->{$cb_left}->{supercontig} eq
		 $briggsae_by_gene->{$cb_right}->{supercontig});
    
    # No briggsae or elegans genes within this span?
    # No sense in evaluating.
    next if (scalar @$cb_candidates == 0 || scalar @ce_candidates == 0);
    
    print_debug("GENES IN: ce/cb: $total_ce/$total_cb");
    
    # 4. Iterate through each of the briggsae genes
    #    Do any of them in the interval hit genes in the elegans interval?
    evaluate_interval($cb_candidates,\@ce_candidates);
  }
}


sub build_ce_range {
  my $anchor = shift;
  my $position = $elegans->{$anchor}->{genome_start};
  my $left_limit  = $position - ($range/2);
  my $right_limit = $position + ($range/2);
  
  my $chrom = $elegans->{$anchor}->{chrom};
  
  # Iterate over all of the genes in this chrom to see if they
  # fall in the span
  # This is extremely inefficient...
  my @ce_candidates;
  foreach (@ordered_elegans) {
    if ($elegans->{$_}->{genome_start} > $left_limit && $elegans->{$_}->{genome_stop} < $right_limit) {
      # Skip this gene if it's already assigned or newly assigned
      next if ($new_orthologs->{$_});
      next if ($orthologs->{$_});
      # Crap this is inefficient - skip to the next gene if
      # we are not on the appropriate chromosome
      next if $elegans->{$_}->{chrom} ne $chrom;
      push @ce_candidates,$_;
    }
  }
  return \@ce_candidates;
}


sub build_cb_range {
  my $anchor = shift;
  my $super  = $briggsae_by_gene->{$anchor}->{supercontig};
  my $position = $briggsae_by_gene->{$anchor}->{start};
  
  my $left_limit  = $position - ($range/2);
  my $right_limit = $position + ($range/2);
  $left_limit = ($left_limit < 0 ) ? $left_limit = 0 : $left_limit;
  
  my @cb_genes = @{$briggsae->{$super}};
  my @cb_candidates;
  foreach my $candidate (@cb_genes) {
    # Already assigned to an ortholog pair? Ignore it...
    my $gene = $candidate->{gene};
    # Skip this gene if its a new or pre-existing ortholog
    next if ($orthologs->{$gene});
    next if ($new_orthologs->{$gene});
    
    if ($candidate->{start} > $left_limit && $candidate->{stop} < $right_limit) {
      print_debug("\t\t$gene " . $candidate->{stop} . ' ' . $candidate->{start});
      push (@cb_candidates,$candidate);
    }
  }
  
  # This length isn't quite right... Doens't handle the ends of chromosomes at all
  return (\@cb_candidates,$right_limit-$left_limit);
}

sub range_query {
  my $position = 0;
  foreach my $ce_left (@$ordered_orthologs) {
    print_alert("------> Evaluating new range centered on $ce_left");
    # Collect all of the genes that are within range of this anchor
    my $ce_candidates = build_ce_range($ce_left);
    $stats->{spans_examined}++;
    
    my $ce_range_start = $elegans->{$ce_left}->{chrom_start} - ($range / 2);
    my $ce_range_stop  = $elegans->{$ce_left}->{chrom_start} + ($range / 2);
    my $ce_range_start = ( $ce_range_start < 0 ) ? 0 : $ce_range_start;
    my $ce_length = $ce_range_stop - $ce_range_stop;
    
    foreach (@$ce_candidates) {
      print_debug("Ce genes in range: $_\t" . $elegans->{$_}->{genome_start} );
    }
    
    # 3. Look up the briggsae mates of the left limit and the right limit
    my $cb_left  = mate($ce_left);
    
    # 4. Fetch all briggsae genes on this supercontig that fall within this range.
    my ($cb_candidates,$cb_length) = build_cb_range($cb_left);
    my $total_ce = scalar @$ce_candidates;
    my $total_cb = scalar @$cb_candidates;

    # Keep track of some span statistics like how many times each gene is sampled
    push (@{$stats->{spans}},[$cb_length,$ce_length,$total_cb,$total_ce]);
    
    # No briggsae or elegans genes within this span? No sense in evaluating...
    next if (scalar @$cb_candidates == 0 || scalar @$ce_candidates == 0);
    
    print_debug("GENES IN: ce/cb: $total_ce/$total_cb");
    
    # 4. Iterate through each of the briggsae genes
    #    Do any of them in the interval hit genes in the elegans interval?
    evaluate_interval($cb_candidates,$ce_candidates);
  }
}



sub dump_new_orthologs {
  # The new orthologs has is bi-directional.
  # Need to compress it into uniques

  # CHECK THIS!
  # This is a bug.  The list is actually getting multiple genes assigned
  # to the same mate during the analysis.
  # Here, they just blow each other away, irrespective of whichever is the best
  
  # The problem, I think, is in the hash lookup tables.
  my $flat = {};
  foreach ( keys %$new_orthologs) {
    my ($cb,$ce,$evalue,$confidence) = @{$new_orthologs->{$_}};
    $flat->{$ce} = [$cb,$ce,$evalue,$confidence];
  }
  #  my $flat = map { $new_orthologs->{$_}->[1] => [ $new_orthologs->{$_} ]  } keys %$new_orthologs;
  foreach (keys %$flat) {
    my ($cb,$ce,$evalue,$confidence) = @{$flat->{$_}};
    print join("\t",$cb,$ce,$evalue,$confidence,$method),"\n"; # placeholder for confidence
  }
}


sub print_stats {
  open OUT,">output/$method.stats";
  print OUT "Orthologs by synteny: $method\n";
  print OUT "Spans examined      : ",$stats->{spans_examined},"\n";
  print OUT "Collinear spans     : ",$stats->{collineary_spans},"\n";
  print OUT "Collinearity breaks : ",$stats->{collinearity_breaks},"\n";

  # Calculate the average span lengths
  my $ce_length = Statistics::Descriptive::Full->new();
  my $cb_length = Statistics::Descriptive::Full->new();
  my $ce_genes_in = Statistics::Descriptive::Full->new();
  my $cb_genes_in = Statistics::Descriptive::Full->new();
  foreach (@{$stats->{spans}}) {
    my ($cb,$ce,$total_cb,$total_ce) = @{$_};
    $ce_length->add_data($ce);
    $cb_length->add_data($cb);
    $ce_genes_in->add_data($total_ce);
    $cb_genes_in->add_data($total_cb);
  }

  print OUT " elegans bp in  : ",$ce_length->sum(),"\n";
  print OUT " briggsae bp in : ",$cb_length->sum(),"\n";
  
  print OUT " average span lengths\n";
  print OUT "\telegans  : " . $ce_length->mean() . ' (SD: ' . $ce_length->standard_deviation() . 
    "SE: " . $ce_length->variance() . ")\n";
  print OUT "\tbriggsae : " . $cb_length->mean() . ' (SD: ' . $cb_length->standard_deviation() . 
    "SE: " . $cb_length->variance() . ")\n";
  print OUT "\n";
  
  print OUT " average genes / span\n";
  print OUT "\telegans  : " . $ce_genes_in->mean() . ' (SD: ' . $ce_genes_in->standard_deviation() . 
    "SE: " . $ce_genes_in->variance() . ")\n";
  print OUT "\tbriggsae : " . $cb_genes_in->mean() . ' (SD: ' . $cb_genes_in->standard_deviation() . 
    "SE: " . $cb_genes_in->variance() . ")\n";
  print OUT "\n";
  
  print OUT "Extent of sampling :\n";
  my $ce_average = fetch_average_sampling('ce_sampled');
  my $cb_average = fetch_average_sampling('cb_sampled');  
  print OUT "\tC. elegans genes sampled  : ",scalar keys %{$stats->{ce_sampled}},"\n";
  print OUT "\t\teach gene sampled on average $ce_average\n";
  
  print OUT "\tC. briggsae genes sampled : ",scalar keys %{$stats->{cb_sampled}},"\n";
  print OUT "\t\teach gene sampled on average $cb_average\n";

  print OUT "\n\nSPANS: ce_length cb_length ce_genes_in cb_genes_in\n";
  foreach (@{$stats->{spans}}) {
    my ($cb_length,$ce_length,$total_cb,$total_ce) = @{$_};
    print OUT join("\t",$ce_length,$cb_length,$total_ce,$total_cb),"\n";
  }
  close OUT;
}

sub fetch_average_sampling {
  my $tag = shift;
  my $total;
  foreach (keys %{$stats->{$tag}}) {
    $total += $stats->{$tag}->{$_};
  }
  my $avg = $total / scalar keys %{$stats->{$tag}};
  return $avg;
}

sub mate      { 
  my $index = shift;
  my $mate = $orthologs->{$index}->{mate};
  return $mate;
}

sub get_index { 
  my $index = shift;
  my $val = $orthologs->{$index}->{index};
  return $val;
}

sub get_range {
  my ($start1,$stop1,$start2,$stop2) = @_;
  my @vals;
  push (@vals,$start1,$stop1,$start2,$stop2);
  my @sorted = sort { $a <=> $b } @vals;
  my $range_start  = $sorted[0];
  my $range_stop   = $sorted[-1];
  my $length = $range_stop - $range_start;
  return ($range_start,$range_stop,$length);
}

sub get_cb_candidates {
  my ($left,$right) = @_;
  my $cb_left = $briggsae_by_gene->{$left};
  my $cb_right = $briggsae_by_gene->{$right};
  my ($cb_range_start,$cb_range_stop,$cb_length) = get_range($cb_left->{start},
							     $cb_left->{stop},
							     $cb_right->{start},
							     $cb_right->{stop});

  my $supercontig = $briggsae_by_gene->{$left}->{supercontig};
  print_alert("\t--> Fetching cb candidates: $left $right $supercontig ($cb_range_start $cb_range_stop)");
  my @cb_genes = @{$briggsae->{$supercontig}};
  
  my @cb_candidates;
  foreach my $candidate (@cb_genes) {
    # Already assigned to an ortholog pair? Ignore it...
    my $gene = $candidate->{gene};
    next if ($orthologs->{$gene});
    next if ($new_orthologs->{$gene});
    if ($candidate->{start} > $cb_range_start && $candidate->{stop} < $cb_range_stop) {
      # Is there an intervening ortholog within this span?
      # If so, let's ignore it since this is a break in collinearity
      # I'll independently assay these spans in the range-based queries
      return (undef,undef,'break in collinearity')
	if ($orthologs->{$gene} || $new_orthologs->{$gene});

      print_debug("\t\t$gene " . $candidate->{stop} . ' ' . $candidate->{start});
      push (@cb_candidates,$candidate);
    }
  }
  return (\@cb_candidates,$cb_length);
}


sub evaluate_interval {
  print_alert("\t--> Evaluating interval...");
  my ($cb_interval,$ce_interval) = @_;
  
  # How many of the elegans genes in this interval does
  # this cb gene hit?

  # First, create a hash lookup for the elegans candidates
  my %ce_lookup = map { $_ => 1 } @$ce_interval;

  foreach (@$ce_interval) {
    $stats->{ce_sampled}->{$_}++;
  }

  # Keep track of all the genes in this interval and
  # all the genes that they hit in elegans
  my $possible_mates = {};
  
  foreach my $cb (@$cb_interval) {
    my $cb_gene = $cb->{gene};

    # Keep track of how many times each gene is sampled.
    $stats->{cb_sampled}->{$cb_gene}++;

    # 1. Fetch all of its blast hits
    print_alert("\t--> Checking for blast hits. query: $cb_gene");
    
    # Skip this gene if we've already assigned it to a new ortholog
    next if ($new_orthologs->{$cb_gene});
    
    my @hits = eval { @{$briggsae_blast->{$cb_gene}} };
    foreach my $hit (@hits) {
      my ($subject,$evalue) = @{$hit};
      
      # Skip this gene if we've already assigned its subject to an ortholog
      # This doesn't quite work out as expected...
      next if ($new_orthologs->{$subject});
      
      # 2. Does this subject fall in the interval?
      if (defined $ce_lookup{$subject}) {
	
	# 1. Does this hit fall below the maximum evalue?
	next unless ($evalue < $max_eval);
	
	# Save this query/subject pair as possible ortholog pair
	# keyed by the briggsae name
	push (@{$possible_mates->{$cb_gene}},[$cb_gene,$subject,$evalue]);
      }
    }
  }
  
  # How many briggsae genes hit elegans genes?
  # 1. Easiest case: One briggsae gene has some hits in elegans
  if (scalar keys %{$possible_mates} == 1) {
    foreach my $cb (keys %{$possible_mates}) {
      # Pass in the list of all the hits for this gene
      evaluate_one_to_ones(@{$possible_mates->{$cb}});
    }
  } else {
    # Multiple briggsae genes have hits in this interval
    # How many elegans genes are hit?
    my %elegans_count;
    foreach my $cb (keys %{$possible_mates}) {
      foreach my $hit (@{$possible_mates->{$cb}}) {
	my ($cb,$ce,$eval,@rest) = @{$hit};
	$elegans_count{$ce}++;
      }
    }
    # 2. Do they all hit the same gene (many-to-ones)?
    if (scalar keys %elegans_count == 1) {
      # This is a true many to one relationship
      # Simply compare the difference of the top hits to see if they are significant
      evaluate_many_to_ones(\%{$possible_mates});
      next;
    } elsif (scalar keys %elegans_count == scalar keys %$possible_mates) {
      # 3. Multiple one-to-ones within the interval
      #    The number of briggsae genes with hits in interval is the 
      #    same as the number of elegans genes seen.
      #    We have a bunch of one-to-one relationships within this interval
      #    Cool.  Let's just evaluate each in turn.
      
      # THIS IS NOT WORKING!  Many-to-ones are slipping in here.
      # The problem here is that I don't discriminate between
      # the quality of the hits, sometimes ending up with multiple
      # briggsae assigned to the same elegans gene
      foreach my $cb (keys %{$possible_mates}) {
	# evaluate_one_to_ones(@{$possible_mates->{$cb}});
      }
    } else {
      # The final case: we have multiple briggsae genes that hit multiple
      # elegans genes all in the same interval perhaps mixed in with true
      # one to ones.
      # These are pretty hairy to deal with making unambiguous assignment
      # of orthology difficult.
      next;
    }
  }
}


sub new_ortholog_pair {
  my ($temp,$confidence) = @_;
  my ($cb,$ce,$evalue) = @$temp;
  print_success("----- FOUND A NEW ORTHOLOG PAIR $cb $ce $evalue ------","$cb\t$ce");
  $new_orthologs->{$cb} = [$cb,$ce,$evalue,$confidence];
  $new_orthologs->{$ce} = [$cb,$ce,$evalue,$confidence];
}


sub compare_top_hits {
  my $hits = shift;
  print_alert("\t\t--> Comparing top hits...");

  # Sort all of the hits based on their evalue
  # Q: IS THIS THE RIGHT SORT ORDER?
  foreach (@$hits) {
    print_debug($_ . "\n");
  }
  my @sorted = sort { $a->[2] <=> $b->[2] } @$hits;
  my $best = $sorted[0]->[2];
  my $next = $sorted[1]->[2];
  compare_evals($best,$next,$sorted[0]);
}

# Best hit should be aray reference
# of cb,ce,eval
sub compare_evals {
  my ($best,$next,$best_hit) = @_;
  my $factor;
  my $confidence;
  if ($best == 0 && $next != 0) {
    # Code from avril -
    # What is happening here?  Are we just setting an arbitrarily large number??
    $factor     = 10000000000000000000000000000000000000000000000;
    $confidence = 0;
  } elsif ($best == 0 && $next == 0) {
    $factor     = 1;
  } else {
    $factor = $next/$best;
    $confidence = log10($next/$best);
  }
  
  if ($factor < $cutoff) {
    return ("? (several equally good hits)");
  } else {
    new_ortholog_pair($best_hit,$confidence);
  }
}



# Compare cases where I have found multiple briggsae genes that
# hit a single elegans gene
sub evaluate_many_to_ones {
  my $candidates = shift;
  print_alert("\t\t--> Evaluating many to ones...");
  # Fetch out all hits for all genes, then pick the two best hits.
  # scope the key for access later in sub
  my $key;
  my @all_hits;
  foreach my $cb (keys %$candidates) {
    foreach my $hit (@{$candidates->{$cb}}) {
      push (@all_hits,$hit);
      my $key  = $hit->[1];
    }
  }
  
  # Now sort the hits and fetch the top hits.
  my @sorted = sort { $a->[2] <=> $b->[2] } @all_hits;
  my $best = $sorted[0]->[2];
  my $next = $sorted[1]->[2];
  compare_evals($best,$next,$sorted[0]);
  return;
}


# Pass in an array of hits to evaluate
sub evaluate_one_to_ones {
  my @hits = @_;
  print_alert("\t\t--> Evaluating possible one-to-one relationship..");
  #  Does it only hit a single gene? (one-to-ones within the interval)
  # We've found a very clear cut new ortholog pair
  if (scalar @hits == 1) {
    my $hit  = @hits[0];
    # Arbitrarily set the confidence to 200 - not sure if this is appropriate
    my $confidence = 200;
    new_ortholog_pair($hit,$confidence);
  } else {
    # This single gene hits multiple genes in the interval below max_eval
    # Compare the top hits and see if I can unambiguosly assign an ortholog pair
    compare_top_hits(\@hits);
  }
  return;
}

sub print_alert {
  my $msg = shift;
  print STDERR "$msg\n" if ($ALERT);
  #  print STDERR -t STDOUT && !$ENV{EMACS} ? "\r" : "\n";
  return;
}

sub print_debug {
  my $msg = shift;
  print STDERR "\t$msg\n" if ($DEBUG);
}

sub print_success {
  my ($msg,$debug) = @_;
  print STDERR "\t$msg\n"  if ($SUCCESS);
  print ERROR "$debug\n" if ($SUCCESS);
}

# I am differentially fetching start and stop.
# For elegans, I am fetching chrom_start
# for briggsae, start....
sub get_positional_keys {
  my $source = shift;
  my ($start_key,$stop_key);
  if ($source eq 'elegans') {
    $start_key = 'chrom_start';
    $stop_key  = 'chrom_stop';
  } else {
    $start_key = 'start';
    $stop_key  = 'stop';
  }
  return ($start_key,$stop_key);
}



sub order_orthologs {
  print_alert("----> Ordering orthologs...");
  my $elegans_ordered = shift;
  my @ordered_orthos;
  my $position = 0;
  foreach (@$elegans_ordered) {
    if (defined $orthologs->{$_}) {
      push (@ordered_orthos,$_);
      # save the ordered elegans position of this ortholog
      $orthologs->{$_}->{index} = $position;
    }
    $position++;
  }
  return (\@ordered_orthos);
}





# Calculates the confidence threshold for ortholog pairs
# The arbitray best score is 200.
# Expects two values, fwd and bwd, which refer to each member of an ortho pair
# sub calc_confidence {
#   my ($fwd,$bwd) = @_;
#   my $top = '200';
  
#   # Get out the evalues for this top hit
#   # 'Natch, since these are already best mutuals,
#   # These two values should be the same (or clase to the same).
  
#   my $fwd_best_eval = $hits->{evalues}->{$fwd . '_' . $bwd};
#   my $bwd_best_eval = $hits->{evalues}->{$bwd . '_' . $fwd};
  
#   # Get out the second best_hits for both directions...
#   my ($fwd_next,$fwd_next_eval,$bwd_next,$bwd_next_eval) = get_second_bests($fwd,$bwd);
  
#   # Just return the top confidence value if there were no second hits
#   return ($top) if (!$fwd_next_eval || !$bwd_next_eval);
  
#   # If the expect value is 0.0, return arbitrary high score
#   return $top if ($fwd_best_eval eq '0.0' || $fwd_best_eval == 0);
#   return $top if ($bwd_best_eval eq '0.0' || $bwd_best_eval == 0);
  
#   # calculate the confidence value
#   # It should be the log of the second best hit over the next best hit
#   my $fwd_conf = log10( $fwd_next_eval / $fwd_best_eval );
#   my $bwd_conf = log10( $bwd_next_eval / $bwd_best_eval );
#   my $conf = ($fwd_conf < $bwd_conf) ? $fwd_conf : $bwd_conf;
#   return $conf;
# }

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}
