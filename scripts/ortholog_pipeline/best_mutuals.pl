#!/usr/bin/perl
# Assign best mutual match blast hits between two genomes
# Todd Harris (harris@cshl.org)
# 5.2003

=pod

=head1 NAME

  Calculate best reciprocal match ortholog pairs from raw blast output

=head1 DESCRIPTION

  This script determines clear best mutual match ortholog pairs
  between two genomes.  It expects standard tab-delimited blast output
  generated using the -m 8 command line flag.

=head1 SYNOPSIS

    best_mutuals
       --blast1     tab delimited blast output, A->B
       --blast2     tab delimited blast output, B->A
       --max_eval   The maximum e-value for significance of blast hits
       --method     Purely descriptive. Typically seg-on or seg-off
       --orphans    Optional boolean. Search for orphans
       --fasta1     Optional. Fasta protein file for searching for orphans
       --fasta2     Optional. Fasta protein file for searching for orphans
       --orthologs  Optional. Path to predetermined orthologs to exclude
       --cutoff     Optional. The cutoff difference between high scoring pairs
       --help       Display this documentation

   The file provided in blast1 will be treated as the reference file
   in the output.  If the --orphans flag is provided, orphans, or genes
   that have no hit below the max_eval, will be dumped out.  This requires
   seperate lists of the genes in each genome. Currently the script expects
   the fasta file that were used for the blast.

   A file of orthologs previously assigned can be passed to exclude
   these genes from consideration in additional ortholog assignments.
   Typically, this is used to step through various blast analyses of
   decreasing sensitivity.  For example, blast seg-on analyses may be
   parsed first without passing the orthologs flag, followed by
   parsing of seg-off blast analyses, where the orphan assignments
   from the seg-on analysis are provided. The format of the file
   should minimally be tab-delimited:
                    gene1    gene2 ...
   Both genes in the pair will be excluded from new ortholog assignments.

   The cutoff switch controls the threshold between two high-scoring
   blast hits.  The difference between the expectation value of the top
   hit and the second best hit must be greater than the --cutoff value
   in order to be assigned to an ortholog pair. The default value is 10^5.

=head1 AUTHOR

  Todd Harris (harris@cshl.org) 
  DATE : 05 May 2003 
  VERSION : $Id: best_mutuals.pl,v 1.2 2004-07-06 14:26:39 ar2 Exp $

=head1 TODO

=cut

#'

use Getopt::Long;
use Pod::Usage;
use Bio::Seq;
use Bio::SeqIO;
use IO::File;
use strict;

my $DEBUG = 0;

my ($blast1,$blast2,$max_eval,$orphans,$fasta1,$fasta2,$method,$orthologs,$cutoff,$help);
GetOptions('blast1=s'   => \$blast1,
	   'blast2=s'   => \$blast2,
	   'max_eval=s' => \$max_eval,
	   'orphans=s'  => \$orphans,
	   'fasta1=s'   => \$fasta1,
	   'fasta2=s'   => \$fasta2,
	   'method=s'   => \$method,
	   'orthologs=s'=> \$orthologs,
	   'cutoff=s'   => \$cutoff,
	   'help'       => \$help);
pod2usage(-verbose=>2) if $help;
pod2usage(0) unless ($blast1 && $blast2 && $max_eval && $method);

my %OLD_ORTHOLOGS = ();
my $hits = {};   # A hash reference to store best and second best hits...
my $cutoff ||= 100000;


# =====================================
#  MAIN
# =====================================
parse_orthologs() if ($orthologs);

# Sequentially parse the reciprocal blasts, tracking top hits for each gene
parse_blast($blast1,'elegans');
parse_blast($blast2,'briggsae');

# Reiterate through all of the blast hits in elegans
# attempting to find mutual best match orthologs

# This is kinda silly.  Should just read the keys of %BEST
# except that it contains values from both genomes.  Oh well.
my %SEEN = ();
open(IN,$blast1) || die "ERROR: cannot open $blast1.\n";
while (my $line = <IN>)  {
  my ($query,$subject,$evalue) = parse_line($line);
  
  # Have I come across this potential ortholog companion before?
  # If so, I've already looked for its ortholog mate.
  # Don't try again.
  next if ($SEEN{$subject});
  
  my ($ortholog,$confidence) = find_ortholog($subject);
  if ($ortholog =~ /\?/) {
    print STDERR "none $subject $ortholog\n" if ($DEBUG);
  } else {
    # Fetch the evalue for this ortholog pair
    my $evalue = get_evalue('best',$ortholog);
    print join("\t",$subject,$ortholog,$evalue,$confidence,$method),"\n";
  }
  $SEEN{$subject}++;
}
close(IN);


if ($orphans) {
  dump_orphans($fasta1,'elegans');
  dump_orphans($fasta2,'briggsae');
}



# =====================================
#  Subroutines
# =====================================
# Build a lookup table of previously found orthologs that should be excluded here.
sub parse_orthologs {
  open IN,$orthologs;
  while (my $line = <IN>) {
    chomp $line;
    my ($gene1,$gene2,$eval,$confidence,$method) = split("\t",$line);
    $OLD_ORTHOLOGS{$gene1}++;
    $OLD_ORTHOLOGS{$gene2}++;
  }
}

# Pass a flag to indicate the directionality of the blast.
# Expects either elegans or briggsae (whatever was used as the query)
sub parse_blast {
  my ($file,$flag) = @_;
  
  open(IN,"$file") || die "Couldn't open $file: $!\n";
  while (my $line = <IN>) {
    my ($query,$subject,$evalue) = parse_line($line);
    $query =~ tr/[a-z]/[A-Z]/;
    
    # Parse out the true elegans name from alternative splices
    # Need to flip flop based on directionality of the blast
    
    # Keep the splice variant identifier with the hit
    if (0) {
      if ($flag eq 'elegans') {
	$query = ($query =~ /(.*\.\d+).*/) ? $1 : $query;
      } else {
	$subject = ($subject =~ /(.*\.\d+).*/) ? $1 : $subject;
      }
    }
    store_best_hit($query,$subject,$evalue);
  }
  close(IN);
}


sub parse_line {
#83268048        CE00003 14      144     86      214     B0041.7 2       91      3.214670        19      14,86,2:121,191,0
  my $line = shift;
  chomp $line;
  my @temp = split(/\s+/,$line);
  
  my $evalue  = 10 * exp(-$temp[10]);
  if (substr($evalue,0,1) eq 'e') { $evalue = "1".$evalue;}
  my $query   = $temp[6];
  my $subject = $temp[1];
  return ($query,$subject,$evalue);
}


sub store_best_hit {
  my ($query,$subject,$evalue) = @_;
  # If we do not already have a top hit for this query,
  # enter it into the evalues hash (keyed with
  # a concatenation of the query and subject IDs)
  
  # Is the evalue for this hit above the signficance cutoff?
  return if ($evalue > $max_eval);
  
  # Are we parsing seg off analysis?
  # If so, has the current gene been assigned to an ortholog pair?
  # Let's just return and not consider this gene...
  if ($orthologs) {
    return if (exists $OLD_ORTHOLOGS{$query});
  }
  
  # Have we seen a best hit for this query yet?
  # If not store the best hit and the evalue for this pair.
  my $best_hit = fetch('best',$query);
  if (!($best_hit)) {
    stuff('best',$query,$subject,$evalue);
  } else {
    # We already have a top hit for this query
    my $best_evalue = get_evalue('best',$query);
    
    # Is this evalue lower than that already stored for the top hit?
    # If so, replace the top hit
    if ($evalue < $best_evalue) {
      stuff('best',$query,$subject,$evalue);
    } else {
      # The evalue of the current pair is higher than that of the best.
      # If we I haven't found a second best hit yet, then store this value
      # (but only if this isn't just another hit)
      my $second_best_hit = fetch('second',$query);
      if (!($second_best_hit) && $subject ne $best_hit) {
	stuff('second',$query,$subject,$evalue);
      } elsif ($second_best_hit) {
	# We've seen a second_best_hit
	# First, is it a duplicate/contained hit?

	if ($subject ne $best_hit && $subject ne $second_best_hit) {
	  # Is the current evalue lower than that of the one stored for
	  # the second best hit?  If so, replace the second best_hit evalue
	  my $second_best_evalue = get_evalue('second',$query);
	  if ($evalue < $second_best_evalue) {
	    stuff('second',$query,$subject,$evalue);
	  }
	}
      }
    }
  }
}


# Find the ortholog pair
sub find_ortholog {
  my $test_gene = shift;

  # Is there a best hit assigned for this gene?
  my $best_hit = fetch('best',$test_gene);
  if ($best_hit) {
    # Testing for reciprocity
    # Was the best hit of best_hit the test_gene?
    
    my $reciprocal_best = fetch('best',$best_hit);
    if (!$reciprocal_best) {
      # There is not hit for the best hit...wierd...
      return ("? (hits a gene with no hits)");
      
      # What is the best_hit for best_hit?  Ugh, the names...
    } elsif ($reciprocal_best eq $test_gene) {

      # THE BEST HIT OF $gene IS $best_hit, AND THE BEST HIT OF $best_hit IS $gene.
      # Check second best hits and see if their e-vals are less than $cutoff
      my $flag = compare_second_best_hits($test_gene,$best_hit);
      return ($flag) if ($flag);
      
      # Now look in the opposite direction.
      # Does the evalue of the best hit of $best_hit at least $cutoff_factor?
      my $flag = compare_second_best_hits($best_hit,$test_gene);
      return ($flag) if ($flag);
      
      # IF WE HAVE GOT THIS FAR, WE TAKE $gene AND $best_hit TO BE ORTHOLOGS:
      my $ortholog = $best_hit;
      
      # Calculate the confidence for this orthologous pair
      my $confidence = calc_confidence($best_hit,$test_gene);
      return ($ortholog,$confidence);
    } elsif ($reciprocal_best ne $test_gene) {
      return ("? ($test_gene hits a gene $best_hit that has a better hit $reciprocal_best");
    }
  } else {
    # Genes that have no hits below MAX_EVAL will end up here
    # since they will not be found in the BEST_HITS hash
    return ("? No hits < $max_eval");
  }
  
  # Should never really end up here
  return;
}

sub compare_second_best_hits {
  my ($gene,$best_hit) = @_;
  my $best_evalue = $hits->{evalues}->{$gene . '_' . $best_hit};
  my $second_best_hit = fetch('second',$gene);
  if ($second_best_hit) {
    my $second_best_evalue = $hits->{evalues}->{$gene . '_' . $second_best_hit};
    my $factor;
    if ($best_evalue == 0 && $second_best_evalue != 0) {
      # What is happening here?  Are we just setting an arbitrarily large number??
      $factor     = 10000000000000000000000000000000000000000000000;
    } elsif ($best_evalue == 0 && $second_best_evalue == 0) {
      $factor     = 1;
    } else {
      $factor     = $second_best_evalue/$best_evalue;
    }
    if ($factor < $cutoff) {
      return ("? (several equally good hits)");
    }
  }
  return;
}



sub dump_orphans {
  my ($file,$tag) = @_;
  
  my $seqin = Bio::SeqIO->new(-format=>'Fasta',-file=>"$file");
  open OUT,">$tag-orphans.out";
  my $date = `date`;
  print OUT "// $tag orphans\n";
  print OUT "// Generated: $date";
  print OUT "// Threshold: $max_eval\n";
  while (my $seqobj = $seqin->next_seq()) {

    # Create a BioSeqIO out object
    my $id = $seqobj->id;
    #  my $clean_id = ($id =~ /(.*\.\d+).*/) ? $1 : $id;

    # Does this ID exist in the BEST_HITS hash?  If not,
    # it's an orphan, baby
    if (!$hits->{$id}) {
      print OUT $id,"\n";
    }
  }
}


# Calculates the confidence threshold for ortholog pairs
# The arbitray best score is 200.
# Expects two values, fwd and bwd, which refer to each member of an ortho pair
sub calc_confidence {
  my ($fwd,$bwd) = @_;
  my $top = '200';
  
  # Get out the evalues for this top hit
  # 'Natch, since these are already best mutuals,
  # These two values should be the same (or clase to the same).
  
  my $fwd_best_eval = $hits->{evalues}->{$fwd . '_' . $bwd};
  my $bwd_best_eval = $hits->{evalues}->{$bwd . '_' . $fwd};
  
  # Get out the second best_hits for both directions...
  my ($fwd_next,$fwd_next_eval,$bwd_next,$bwd_next_eval) = get_second_bests($fwd,$bwd);
  
  # Just return the top confidence value if there were no second hits
  return ($top) if (!$fwd_next_eval || !$bwd_next_eval);
  
  # If the expect value is 0.0, return arbitrary high score
  return $top if ($fwd_best_eval eq '0.0' || $fwd_best_eval == 0);
  return $top if ($bwd_best_eval eq '0.0' || $bwd_best_eval == 0);
  
  # calculate the confidence value
  # It should be the log of the second best hit over the next best hit
  my $fwd_conf = log10($fwd_next_eval / $fwd_best_eval );
  my $bwd_conf = log10($bwd_next_eval / $bwd_best_eval );
  my $conf = ($fwd_conf < $bwd_conf) ? $fwd_conf : $bwd_conf;
  #  if ($conf > 200) {
  #    print_debug("$fwd,$bwd,$fwd_next_eval,$fwd_best_eval,$bwd_next_eval,$bwd_best_eval,$conf");
  #  }
  
  my $format = sprintf("%.3f",$conf);
  return $format;
}

sub log10 {
  my $n = shift;
  return log($n)/log(10);
}


# Return the second best subjects and evalues
sub get_second_bests {
  my ($fwd,$bwd) = @_;

  # There may not be seconds forwards OR bwds
  my $next_fwd_subject = $hits->{$fwd}->{second};
  my $next_fwd_eval    = $hits->{evalues}->{$fwd . '_' . $next_fwd_subject};
  
  my ($next_bwd_eval,$next_bwd_subject);
  if ($bwd) {
    $next_bwd_subject = $hits->{$bwd}->{second};
    $next_bwd_eval    = $hits->{evalues}->{$bwd . '_' . $next_bwd_subject};
  }
  return ($next_fwd_subject,$next_fwd_eval,$next_bwd_subject,$next_bwd_eval);
}




# Data storage and access methods
sub stuff {
  my ($tag,$query,$subject,$evalue) = @_;
  $hits->{$query}->{$tag} = $subject;
  $hits->{evalues}->{$query . '_' . $subject} = $evalue;
}

sub fetch {
  my ($tag,$query) = @_;
  my $value = eval { $hits->{$query}->{$tag} };
  return $value;
}

sub get_evalue {
  my ($tag,$query) = @_;
  # Fetch the subject since evalues are stored as query_subject concatenations
  my $subject = $hits->{$query}->{$tag};
  my $evalue  = $hits->{evalues}->{$query . '_' . $subject};
 return $evalue;
}


#package Hits;
#
#
#sub new {
#  my $self = shift;
#  my $this = bless {},$self;
#  return $this;
#}
#
#
#sub best_hit {
#  my ($self,$query) = @_;
#  return eval { $blast_hits->{$query}->{best_hit}->{subject} };
#}
#
#1;

