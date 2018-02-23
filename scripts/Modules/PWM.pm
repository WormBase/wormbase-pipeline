#!/software/bin/perl -w
#
# PWM.pm                       
# 
# by Gary Williams
#
# Do position scoring matrix evaluation (PWM, PSM, PSSM) of a sequence.
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2008-06-27 15:39:08 $      

=pod

=head1 NAME

 PWM

=head1 SYNOPSIS

 my $pwm = PWM->new;

 my $result = splice3($sequence, $position);

 my $result = splice5($sequence, $position);

 my $result = atg($sequence, $position);

=head1 DESCRIPTION

    This uses position weight matrices (PWM, PSM, PSSM) to evaluate a
    score for a base being the site of a 3' splice site, a 5' splice
    site or a START codon.

    It uses the matrix data from ACEDB in the /wgf directory.

    The correction of the logodds matrix with background probabilities
    is done using a simple constant ratio for each of the four bases.

    The position in the sequence uses normal perl coordinates starting
    at zero.

    Splice cut sites occur just after the base being looked at (in
    the reverse sense it is just before the base being looked at.)

    These routines will fail it anything other than A,C,G,T,N,- is in
    the sequence.

=head1 CONTACT

Gary gw3@sanger.ac.uk

=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are preceded with a _

=cut




package PWM;

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Carp;

use POSIX qw(log10);


=head2 

    Title   :   new
    Usage   :   my $pwm = PWM->new;
    Function:   initialises the matrices for the PWM object
    Returns :   PWM object;
    Args    :   

=cut

sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;

  # initialise the matrices
  $self->_init_background_probs;
  $self->_init_splice3;
  $self->_init_splice5;
  $self->_init_atg;

  return $self;
}

=head2 

    Title   :   _init_background_probs
    Usage   :   $self->_init_background_probs;
    Function:   initialises the background fractions hash for the four bases
    Returns :   
    Args    :   

=cut

sub _init_background_probs {
  my $self = shift;

  # relative frequencies of background sequence:
  # A 26290
  # C 14383
  # G 14380
  # T 26285
  my $total = 26290+14383+14380+26285;
  my %background_prob;
  $background_prob{'a'} = 26290/$total;
  $background_prob{'c'} = 14383/$total;
  $background_prob{'g'} = 14380/$total;
  $background_prob{'t'} = 26285/$total;
  $background_prob{'n'} = 0;
  $background_prob{'-'} = 0;
    
  $self->{background_prob} = \%background_prob;
}

=head2 

    Title   :   _init_splice3
    Usage   :   $self->_init_splice3;
    Function:   initialises the counts matrix for the 3-prime splice sites
    Returns :   
    Args    :   

=cut

sub _init_splice3 {
  my $self = shift;

  #siteType: intron3  refSeqs: genes  freqType: within  classDef: unique
  #startOff: -25  endOff: 5  numSymbs: 1  maxSymb: 5  numForced: 2

  my %matrix3;
  $matrix3{'startoffset'} = -24;	# offset from splice site to start of matrix - shifted so that splice is after the high-scoring base
  $matrix3{'endoffset'} = 6;	# offset from splice site to end of matrix - shifted so that splice is after the high-scoring base
  $matrix3{'a'} = [2999, 3044, 3078, 3150, 3056, 3181, 3420, 3747, 4225, 4182, 3637, 3108, 2752, 2799, 2682, 2960, 3276, 3516, 2313, 476, 67, 757, 240, 8192, 0, 3359, 2401, 2514, 2240, 2412, 2279];
  $matrix3{'c'} = [1328, 1289, 1154, 1135, 1244, 1169, 1175, 1046, 920, 1026, 999, 1113, 1089, 1076, 1241, 1223, 970, 648, 664, 236, 129, 1109, 6830, 0, 0, 1277, 1533, 1847, 2101, 2101, 2250];
  $matrix3{'g'} = [804, 791, 850, 798, 765, 685, 639, 570, 535, 461, 576, 654, 715, 775, 626, 606, 593, 575, 516, 144, 39, 595, 12, 0, 8192, 2539, 1301, 1567, 1861, 1517, 1568];
  $matrix3{'t'} = [3061, 3068, 3110, 3109, 3127, 3157, 2958, 2829, 2512, 2523, 2980, 3317, 3636, 3542, 3643, 3403, 3353, 3453, 4699, 7336, 7957, 5731, 1110, 0, 0, 1017, 2957, 2264, 1990, 2162, 2095];

  $self->{splice3} = $self->_get_logodds(%matrix3);
  $self->{rev_splice3} = $self->_get_revcomp($self->{splice3});
}


=head2 

    Title   :   _init_splice5
    Usage   :   $self->_init_splice5;
    Function:   initialises the counts matrix for the 5-prime splice sites
    Returns :   
    Args    :   

=cut

sub _init_splice5 {
  my $self = shift;

  #siteType: intron5, refSeqs: genes, freqType: within  classDef: unique
  #startOff: -5  endOff: 25  numSymbs: 1  maxSymb: 5  numForced: 2

  my %matrix5;
  $matrix5{'startoffset'} = -5;  # offset from splice site to start of matrix
  $matrix5{'endoffset'} = 25;	# offset from splice site to end of matrix
  $matrix5{'a'} = [2347, 2997, 2938, 3404, 4644, 1518, 0, 0, 4836, 5486, 837, 1632, 2189, 2278, 2355, 2493, 2675, 2695, 2718, 2702, 2832, 2727, 2797, 2845, 2796, 2884, 2857, 2799, 2866, 2906, 2884];
  $matrix5{'c'} = [1589, 1419, 1582, 1850, 1224, 583, 0, 14, 118, 588, 237, 801, 771, 889, 986, 1018, 1063, 1028, 1024, 1104, 1068, 1031, 1004, 992, 1000, 985, 1034, 1080, 1177, 1095, 1142];
  $matrix5{'g'} = [1937, 1299, 1372, 1562, 912, 4891, 8192, 0, 1890, 672, 6164, 589, 962, 1056, 827, 998, 1134, 1167, 1092, 1077, 1049, 1124, 1112, 1055, 1078, 1160, 1098, 1095, 1054, 1072, 1098];
  $matrix5{'t'} = [2319, 2477, 2300, 1376, 1412, 1200, 0, 8178, 1348, 1446, 954, 5170, 4270, 3969, 4024, 3683, 3320, 3302, 3358, 3309, 3243, 3310, 3279, 3300, 3318, 3163, 3203, 3218, 3095, 3119, 3068];
  
  $self->{splice5} = $self->_get_logodds(%matrix5);
  $self->{rev_splice5} = $self->_get_revcomp($self->{splice5});
}


=head2 

    Title   :   _init_atg
    Usage   :   $self->_init_atg;
    Function:   initialises the counts matrix for the START codons
    Returns :   
    Args    :   

=cut

sub _init_atg {
  my $self = shift;

  #siteType: atg  refSeqs: genes  freqType: within  classDef: unique
  #startOff: -9  endOff: 11  numSymbs: 1  maxSymb: 5  numForced: 3  forcedPos:   0 1 2 jump: 1

  my %matrix_atg;
  $matrix_atg{'startoffset'} = -9;  # offset from splice site to start of matrix
  $matrix_atg{'endoffset'} = 11;	# offset from splice site to end of matrix
  $matrix_atg{'a'} = [19, 13, 21, 13, 12, 29, 31, 20, 22, 48,  0,  0, 13, 16,  9, 17, 15,  9, 13, 22, 13];
  $matrix_atg{'c'} = [11, 16,  7, 10, 16, 11,  3, 10, 12,  0,  0,  0,  8, 15, 13, 12, 10, 21,  7, 10, 21];
  $matrix_atg{'g'} = [5,   8,  3, 11,  5,  5, 12,  8,  5,  0,  0, 48, 17, 13, 14,  9,  5, 11, 17,  3,  6];
  $matrix_atg{'t'} = [13, 11, 17, 14, 15,  3,  2, 10,  9,  0, 48,  0, 10,  4, 12, 10, 18,  7, 11, 13,  8];

  $self->{atg} = $self->_get_logodds(%matrix_atg);
  $self->{rev_atg} = $self->_get_revcomp($self->{atg});

}



=head2 

    Title   :   _get_logodds
    Usage   :   $log_odds_matrix_href = $self->_get_logodds(%matrix5);
    Function:   gets the log-odds scores matrix from the matrix of counts
    Returns :   the hash-ref of the log-odds matrix
    Args    :   the counts matrix to be converted

=cut

sub _get_logodds {
  my ($self, %count_matrix) = @_;

  # get the log-odds scores matrix from the matrix of counts
  # log (prob_ij / background_prob_i)
  # where:
  # prob_ij is the probability of observing symbol i at position j of the motif
  # background_prob_i is the probability of observing the symbol i in a background model
  my %logodds;
  $logodds{'startoffset'} = $count_matrix{'startoffset'};
  $logodds{'endoffset'} = $count_matrix{'endoffset'};
  my $len = @{$count_matrix{'a'}} - 1;
  foreach my $base ('a', 'c', 'g', 't') {
    foreach my $pos (0 .. $len) {
      my $total = $count_matrix{'a'}->[$pos] + $count_matrix{'c'}->[$pos] + $count_matrix{'g'}->[$pos] + $count_matrix{'t'}->[$pos];
      if ($count_matrix{$base}->[$pos] == 0) {$count_matrix{$base}->[$pos] = 1;} # can't take the log10 of zero, so we make this 1
      my $prob_ij = $count_matrix{$base}->[$pos] / $total;
      $logodds{$base}->[$pos] = POSIX::log10($prob_ij / $self->{background_prob}{$base});
    }
  }
  return \%logodds;
}


=head2 

    Title   :   _get_revcomp
    Usage   :   $self->{rev_splice5} = $self->_get_revcomp($self->{splice5});
    Function:   gets the log-odds matrix for the reverse sense
    Returns :   the hash-ref of the reverse complement log-odds matrix
    Args    :   the hash-ref of the log-odds matrix to be converted

=cut

sub _get_revcomp {
  my ($self, $matrix_href) = @_;
  my %matrix = %{$matrix_href};
  
  my %revcomp;
  my $len = @{$matrix{'a'}} - 1;
  $revcomp{'startoffset'} = - $matrix{'endoffset'};
  $revcomp{'endoffset'} = - $matrix{'startoffset'};
  $revcomp{'a'} = [reverse @{$matrix{'t'}}];
  $revcomp{'c'} = [reverse @{$matrix{'g'}}];
  $revcomp{'g'} = [reverse @{$matrix{'c'}}];
  $revcomp{'t'} = [reverse @{$matrix{'a'}}];

  return \%revcomp;
}



=head2 

    Title   :   _evaluate
    Usage   :   $result = $self->_evaluate($self->{splice3}, $pos, $seq);
    Function:   evaluates the log-odds matrix against the given sequence at the given position
    Returns :   the score of the evaluation
    Args    :   hash-ref of the log-odds matrix, 
                position in the sequence to look at,
                sequence to evaluate, 

=cut

# evaluate the matrix against the sequence at the given position
sub _evaluate {
  my ($self, $matrix_href, $pos) = @_; # don't copy $seq to a local variable, use $_[3] for speed

  my $result = 0;

  if ($pos + $matrix_href->{'startoffset'} < 0) {return -999;}
  if ($pos + $matrix_href->{'endoffset'} > length($_[3]) - 1) {return -999;}

  for (my $i = 0; $i < @{$matrix_href->{'a'}}; $i++) {
    my $next_position = $pos + $matrix_href->{'startoffset'} + $i;
    my $chr = lc substr($_[3], $next_position, 1);
    if ($chr ne 'a' && $chr ne 'c' && $chr ne 'g' && $chr ne 't') {return -999;}
    $result += $matrix_href->{$chr}->[$i];
  }

  return $result;
}

=head2 

    Title   :   splice3_by_ref
    Usage   :   $result = $pwm->splice3_by_ref($position, $sense, $sequence);
    Function:   returns the score for a 3-prime splice site after the given position in the sequence_ref
    Returns :   the splice site score
    Args    :   reference of sequence to evaluate, 
                position in the sequence to look at (starting from zero)
                [optional] sense '+' or '-' ('+' is the default)

=cut

sub splice3_by_ref {
  my ($self, $pos, $sense) = @_;  # don't copy $seq to a local variable, use $_[3] for speed
  if (!defined $sense) {$sense = '+';}

  if ($sense eq '+') {
    return $self->_evaluate($self->{splice3}, $pos, $_[3]);
  } else {
    return $self->_evaluate($self->{rev_splice3}, $pos, $_[3]);
  }
}

=head2 

    Title   :   splice5_by_ref
    Usage   :   $result = $pwm->splice5_by_ref($position, $sense, $sequence);
    Function:   returns the score for a 5-prime splice site after the given position in the sequence_ref
    Returns :   the splice site score
    Args    :   reference of the sequence to evaluate, 
                position in the sequence to look at (starting from zero)
                [optional] sense '+' or '-' ('+' is the default)

=cut

sub splice5_by_ref {
  my ($self, $pos, $sense) = @_;  # don't copy $seq to a local variable, use $_[3] for speed
  if (!defined $sense) {$sense = '+';}

  if ($sense eq '+') {
    return $self->_evaluate($self->{splice5}, $pos, $_[3]);
  } else {
    return $self->_evaluate($self->{rev_splice5}, $pos, $_[3]);
  }
}

=head2 

    Title   :   atg_by_ref
    Usage   :   $result = $pwm->atg_by_ref($position, $sense, $sequence);
    Function:   returns the score for a START codon starting at the given position in the sequence_ref
    Returns :   the START codon score
    Args    :   reference of the sequence to evaluate, 
                position in the sequence to look at (starting from zero)
                [optional] sense '+' or '-' ('+' is the default)

=cut

sub atg_by_ref {
  my ($self, $pos, $sense) = @_;  # don't copy $seq to a local variable, use $_[3] for speed
  if (!defined $sense) {$sense = '+';}

  if ($sense eq '+') {
    return $self->_evaluate($self->{atg}, $pos, $_[3]);
  } else {
    return $self->_evaluate($self->{rev_atg}, $pos, $_[3]);
  }
}

=head2 

    Title   :   splice3
    Usage   :   $result = $pwm->splice3($sequence, $position);
    Function:   returns the score for a 3-prime splice site after the given position in the sequence
    Returns :   the splice site score
    Args    :   sequence to evaluate, 
                position in the sequence to look at (starting from zero)
                [optional] sense '+' or '-' ('+' is the default)

=cut

sub splice3 {
  my ($self, $seq, $pos, $sense) = @_;
  if (!defined $sense) {$sense = '+';}

  if ($sense eq '+') {
    return $self->_evaluate($self->{splice3}, $pos, $seq);
  } else {
    return $self->_evaluate($self->{rev_splice3}, $pos, $seq);
  }
}

=head2 

    Title   :   splice5
    Usage   :   $result = $pwm->splice5($sequence, $position);
    Function:   returns the score for a 5-prime splice site after the given position in the sequence
    Returns :   the splice site score
    Args    :   sequence to evaluate, 
                position in the sequence to look at (starting from zero)
                [optional] sense '+' or '-' ('+' is the default)

=cut

sub splice5 {
  my ($self, $seq, $pos, $sense) = @_;
  if (!defined $sense) {$sense = '+';}

  if ($sense eq '+') {
    return $self->_evaluate($self->{splice5}, $pos, $seq);
  } else {
    return $self->_evaluate($self->{rev_splice5}, $pos, $seq);
  }
}

=head2 

    Title   :   atg
    Usage   :   $result = $pwm->atg($sequence, $position);
    Function:   returns the score for a START codon starting at the given position in the sequence
    Returns :   the START codon score
    Args    :   sequence to evaluate, 
                position in the sequence to look at (starting from zero)
                [optional] sense '+' or '-' ('+' is the default)

=cut

sub atg {
  my ($self, $seq, $pos, $sense) = @_;
  if (!defined $sense) {$sense = '+';}

  if ($sense eq '+') {
    return $self->_evaluate($self->{atg}, $pos, $seq);
  } else {
    return $self->_evaluate($self->{rev_atg}, $pos, $seq);
  }
}


=head2 

    Title   :   setup_matrix
    Usage   :   $pwm->setup_matrix("mymatrix", @alignment);
    Function:   sets up the PWM called $name based on the aligned (no spaces) set of sequences
    Returns :   
    Args    :   name of matrix
                arary of sequences of examples of the motif to use to create the array
=cut
sub setup_matrix {
  my ($self, $name, @seqs) = @_;

  # get the matrix of counts
  my %count_matrix;
  foreach my $seq (@seqs) {
    my $i = 0;
    foreach my $base (split //, $seq) {
      $count_matrix{lc $base}->[$i++]++;
    }
  }
  $count_matrix{'startoffset'} = 0;  # offset from splice site to start of matrix
  $count_matrix{'endoffset'} = @{$count_matrix{'a'}} - 1;	# offset from splice site to end of matrix
  
  # make the log-odds matrix and its reverse complement
  $self->{$name} = $self->_get_logodds(%count_matrix);
  $self->{"rev_$name"} = $self->_get_revcomp($self->{$name});
}

=head2 

    Title   :   get
    Usage   :   $result = $pwm->get("mymatrix", $seq, $pos);
    Function:   returns the score for a named matrix at a specified position
    Returns :   the evaluated matrix score
    Args    :   name of the matrix,
                sequence to evaluate, 
                position in the sequence to look at (starting from zero)
                [optional] sense '+' or '-' ('+' is the default)

=cut

sub get {
  my ($self, $name, $seq, $pos, $sense) = @_;
  if (!defined $sense) {$sense = '+';}

  if ($sense eq '+') {
    return $self->_evaluate($self->{$name}, $seq, $pos);
  } else {
    return $self->_evaluate($self->{"rev_$name"}, $seq, $pos);
  }
}

=head2 

    Title   :   get
    Usage   :   $result = $pwm->atg($sequence, $threshold);
    Function:   returns a list of lists of position and score where the score is >= $threshold
    Returns :   the lists of lists of position and score
    Args    :   name of the matrix,
                sequence to evaluate, 
                threshold of score to return
                [optional] sense '+' or '-' ('+' is the default)

=cut

sub search {
  my ($self, $name, $seq, $threshold, $sense) = @_;
  if (!defined $sense) {$sense = '+';}

  my @result;
  for (my $pos = 0; $pos < length $seq; $pos++) {
    my $result = $self->get($name, $seq, $pos, $sense);
    if ($result >= $threshold) {
      push @result, [$pos, $result];
    }
  }

  return @result;
}


1;
