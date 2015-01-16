=pod 

=head1 NAME

Feature_mapper

=head1 SYNOPSIS

 my $mapper      = Feature_mapper->new($database, $refresh);
 my @map_info    = $mapper->map_feature("AH6","actgtacgtagcgagcaccgatcaggacgag","tgactagcggacagcgagcagctagctgat");
 print "@map_info\n";

 will print; 
 >AH6 3456 3567
 >

 my @genes = $mapper->check_overlapping_CDS("AH6",3456,23567);
 print "@genes\n";

 >AH6.1 AH6.t1 AH6.2
 >

=head1 DESCRIPTION

  This object is used to map features to the genome based on matching flanking sequences.  
  A sequence object is given as a "seed" so that the script knows where to start looking.
  If both flanking sequences cant be found in that seq, the seq is extended by 1000 bp in either direction up to 6000 bp each way.

  Inherits from Sequence_extract

=head1 CONTACT

Anthony  ar2@sanger.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. 
Internal methods are usually preceded with a _

=cut

package Feature_mapper;

use lib $ENV{'CVS_DIR'};

use strict;
use Sequence_extract;
use Carp;
use String::Approx qw(aindex aslice adist);

our (@ISA);

@ISA = ('Sequence_extract');

##########################################

=head2 new

  Title   :   new
  Usage   :   my $mapper = Feature_mapper->new($database,1);
  Returns :   ref to Feature_mapper obj
  Args    :   Database  (optional) - which database to use. Default is current_DB
              refresh   default is NULL - connect to the database and update coordinates
  Requires:    Database must have CHROMOSOMES directory with dna of each chromosome in CHROMOSOME_*.dna (FASTA format)

=cut

sub new {
    my ($class,$database,$refresh,$wormbase) = @_;
    my $self = Sequence_extract->invoke($database, $refresh,$wormbase);

    bless $self, $class;

    $self->wormbase($wormbase);

    return $self;
}

##########################################

=head2 map_feature

  Title   :   map_feature
  Usage   :   my @map_info = $mapper->map_feature("AH6","actgtacgtagcgagcaccgatcaggacgag","tgactagcggacagcgagcagctagctgat");
  Returns :   smallest sequence object that contains the given flanking seq and coords
              of last bp of left flank and first bp of right flank. E.g. for a 4 bp feature:

                        flank_L                          flank_R
              actgtacgtagcgagcaccgatcaggacgagXXXXtgactagcggacagcgagcagctagctgat
                                            ^    ^                       
    
              and for a 0-bp feature:


              actgtacgtagcgagcaccgatcaggacgagTGACTACGTGCTATGCAGCGAGCAT
                                            ^^

              Note therefore that for 0-bp features (e.g. SL-sites), the returned coords define
              the conventional extent of the feature and are thus ready for use; for non-0-bp
              features however, the caller will need to so some +1/-1 adjustment to the coords
              to make them define the extent of the feature itself.
=cut


sub map_feature {
  my ($self, $seq, $flank_L, $flank_R, $min_len, $max_len) = @_;
  
  if (not $flank_L and not $flank_R) {
    die "Cannot call map_feature with both flanks empty\n";
  } elsif (not $flank_L or not $flank_R) {
    return $self->map_feature_single_flank($seq, $flank_L, $flank_R);
  }
  
  my ($left_coord, $right_coord);
  
  # convert to chromosome coords
  my ($chr_id, $chr_st, $chr_en) = $self->LocateSpanUp($seq);
  # for reverse-oriented assembly components, coords come back in reverse order
  if ($chr_en < $chr_st) {
    ($chr_st, $chr_en) = ($chr_en, $chr_st);
  }

  my ($max_extension, $extend_by) = (5, 5000);
  
  for(my $cur_extension = 0; $cur_extension < $max_extension; $cur_extension++) {
    $chr_st -= ($cur_extension * $extend_by);
    $chr_en += ($cur_extension * $extend_by);

    $chr_st = 1 if $chr_st < 1;
    $chr_en = $self->Superlink_length($chr_id) if $chr_en > $self->Superlink_length($chr_id);

    my $reg_len = ($chr_en - $chr_st + 1);

    #print STDERR "Mapping to region of $chr_id, starting at $chr_st and extending for $reg_len bases\n";

    my $dna = $self->Sub_sequence($chr_id, $chr_st - 1, $reg_len);
    my @start_ends = $self->_check_for_match($dna, $flank_L, $flank_R, $min_len, $max_len);

    if (@start_ends) {
      my ($dna_left_coord, $dna_right_coord) = @{$start_ends[0]};
      my ($chr_left_coord, $chr_right_coord) = ($dna_left_coord + $chr_st - 1, $dna_right_coord + $chr_st - 1);

      #print " Feature maps to $chr_id $chr_left_coord, $chr_right_coord\n";
      my $orig_strand = ($chr_right_coord >= $chr_left_coord) ? "+" : "-";
      
      # need to map down to a region that contains the feature AND the flanking sequence
      my ($left_offset, $right_offset);
      
      if ($orig_strand eq '+') { 
        $left_offset = 1 - length($flank_L);
        $right_offset = length($flank_R) - 1;
      } else {
        $left_offset = length($flank_L) - 1;
        $right_offset = 1 - length($flank_R);
      }
      
      #printf " Adjusted to %d %d\n", $chr_left_coord + $left_offset, $chr_right_coord + $right_offset;
      
      my ($seq, $left_coord, $right_coord) = 
          $self->LocateSpan($chr_id, $chr_left_coord + $left_offset, $chr_right_coord + $right_offset);
      
      #printf " Got back %s %d %d\n", $seq, $left_coord, $right_coord;
      
      my $new_strand = ($left_coord <= $right_coord) ? "+" : "-";
      
      if ($orig_strand eq $new_strand) {
        $left_coord -= $left_offset;
        $right_coord -= $right_offset;
      } else {
        $left_coord -= $right_offset;
        $right_coord -= $left_offset;
      }

      #printf " Adjusted to %s %d %d\n", $seq, $left_coord, $right_coord;

      return ($seq, $left_coord, $right_coord);
    }
  }
  
  return 0;
}


=head2 map_feature_single_flank

  Title   :   map_feature
  Usage   :   my @map_info = $mapper->map_feature_to_end("AH6","actgtacgtagcgagcaccgatcaggacgag",undef)
              my @map_info = $mapper->map_feature_to_end("AH6",undef, "actgtacgtagcgagcaccgatcaggacgag")

  Description: this method is used to map features that are defined by a single flank only, where the
               location extend from that flank from the start (or end) of the given sequence

  Returns :   smallest sequence object that contains the given flanking seq and the extent of 
              the feature. The end coordinate will be given as 1-past the end, or 1 before the beginning
              as approprate. 

              The caller will need to so some +1/-1 adjustment to the coords
              to make them define the extent of the feature itself.
=cut

sub map_feature_single_flank {
  my ($self, $seq, $flank_L, $flank_R) = @_;
 
  # convert to chromosome coords
  my ($chr_id, $chr_st, $chr_en) = $self->LocateSpanUp($seq);

  # for reverse-oriented assembly components, coords come back in reverse order
  if ($chr_en < $chr_st) {
    ($chr_st, $chr_en) = ($chr_en, $chr_st);
  }

  my $reg_len = ($chr_en - $chr_st + 1);
  my $dna = lc($self->Sub_sequence($chr_id, $chr_st - 1, $reg_len));
  $flank_L = lc($flank_L);
  $flank_R = lc($flank_R);
  
  if ($flank_L and not $flank_R) {
    # feature extends from left flank to end of region (on forward strand) or start of region (on reverse strand)
    my ($left_coord) = $self->_check_for_match_left($dna, $flank_L, 1);
    if (defined $left_coord) {
      # extends from match position to end of region, on forward strand
      return ($seq, $left_coord, $chr_en + 1);
    } else {
      my $flank_R_rev = $self->DNA_revcomp($flank_L);
      # flank now behaves like a right flank
      ($left_coord) = $self->_check_for_match_right($dna, $flank_R_rev, -1);
      if (defined $left_coord) {
        # extends from match position to start of region, on reverse strand
        return ($seq, $left_coord, $chr_st - 1);
      }
    }
  } elsif ($flank_R and not $flank_L) {
    my ($right_coord) = $self->_check_for_match_right($dna, $flank_R, -1);
    if (defined $right_coord) {
      # extends from start of region to match position, on forward strand
      return ($seq, $chr_st - 1, $right_coord);
    } else {
      my $flank_L_rev = $self->DNA_revcomp($flank_R);
      ($right_coord) = $self->_check_for_match_left($dna, $flank_R, 1);
      if (defined $right_coord) {
        # extends from match position to end of region, on reverse strand
        return ($seq, $right_coord, $chr_en + 1, $right_coord);
      }
    }
  }

  return 0;
}


##########################################

=head2 _check_for_match

 NOTE: This function takes raw DNA seq as a string ie "aagtcaatcggatatgtgatggagctagctgatcgatcgtgc . . . . . etc ";

  Title   :   _check_for_match
  Usage   :   $self->_check_for_match($dna,"actgtacgtagcgagcaccgatcaggacgag","tgactagcggacagcgagcagctagctgat");
  Returns :   coords of last bp of left flank and first bp of right flank. E.g. for a 4 bp feature:

                        flank_L                          flank_R
              actgtacgtagcgagcaccgatcaggacgagXXXXtgactagcggacagcgagcagctagctgat
                                            ^    ^                       
    
              and for a 0-bp feature:


              actgtacgtagcgagcaccgatcaggacgagTGACTACGTGCTATGCAGCGAGCAT
                                            ^^

  Args    :   dna as string, 2 flanking seqs as strings

=cut

sub _check_for_match {
  my ($self, $dna, $flank_L, $flank_R, $min_len, $max_len) = @_;

  my (@good_pairs);
  
  # make sure all in same case !
  $dna = lc($dna);
  $flank_L = lc($flank_L);
  $flank_R = lc($flank_R);
  
  # map longest flank first, to narrow search space
  if (length($flank_L) >= length($flank_R)) {
    foreach my $left_c ($self->_check_for_match_left($dna, $flank_L, 1)) {
      #print STDERR "  found lmatch at coordinate $left_c\n";
      foreach my $right_c ($self->_check_for_match_right($dna, $flank_R, 1, $left_c)) {
        #print STDERR "  found rmatch at coordinate $right_c\n";
        my $flen = $right_c - $left_c - 1;
        
        if ($flen >= 0 and 
            (not defined $min_len or $flen >= $min_len) and
            (not defined $max_len or $flen <= $max_len)) {
          push @good_pairs, [$left_c, $right_c];
        }
      }
    }
  } else {
    foreach my $right_c ($self->_check_for_match_right($dna, $flank_R, -1)) {
      foreach my $left_c ($self->_check_for_match_left($dna, $flank_L, -1, $right_c - length($flank_L))) {
        my $flen = $right_c - $left_c - 1;
        
        if ($flen >= 0 and 
            (not defined $min_len or $flen >= $min_len) and
            (not defined $max_len or $flen <= $max_len)) {
          push @good_pairs, [$left_c, $right_c];
        }
      }
    }
  }

  if (not @good_pairs) {
    # try again on the reverse strand
    #print STDERR "Trying on reverse strand\n";
    my ($flank_L_rev, $flank_R_rev) = ($self->DNA_revcomp($flank_R), $self->DNA_revcomp($flank_L));

    if (length($flank_L_rev) >= length($flank_R_rev)) {
      foreach my $left_c ($self->_check_for_match_left($dna, $flank_L_rev, 1)) {
        #print STDERR "  found lmatch at coordinate $left_c\n";
        foreach my $right_c ($self->_check_for_match_right($dna, $flank_R_rev, 1, $left_c)) {
          #print STDERR "    found rmatch at coordinate $right_c\n";

          my $flen = $right_c - $left_c - 1;
          
          if ($flen >= 0 and 
              (not defined $min_len or $flen >= $min_len) and
              (not defined $max_len or $flen <= $max_len)) {
            push @good_pairs, [$right_c, $left_c];
          }
        }
      }
    } else {
      foreach my $right_c ($self->_check_for_match_right($dna, $flank_R_rev, -1)) {
        foreach my $left_c ($self->_check_for_match_left($dna, $flank_L_rev, -1, $right_c - length($flank_L_rev))) {
          my $flen = $right_c - $left_c - 1;
          
          if ($flen >= 0 and 
              (not defined $min_len or $flen >= $min_len) and
              (not defined $max_len or $flen <= $max_len)) {
            push @good_pairs, [$right_c, $left_c];
          }
        }
      }
    }
  }

  return @good_pairs;
}


sub _check_for_match_left {
  my ($self, $dna, $flank, $direction, $start_looking) = @_;

  my @matches;
  if ($direction < 0) {
    $start_looking = length($dna) if not defined $start_looking;

    while ((my $loc = rindex($dna, $flank, $start_looking)) >= 0) {
      push @matches, $loc + length($flank);
      $start_looking = $loc - 1;
    }
  } else {
    $start_looking = 0 if not defined $start_looking;

    while ((my $loc = index($dna, $flank, $start_looking)) >= 0) {
      push @matches, $loc + length($flank);
      $start_looking = $loc + 1;
    }
  }

  return sort { $a <=> $b } @matches;
}


sub _check_for_match_right {
  my ($self, $dna, $flank, $direction, $start_looking) = @_;

  my @matches;
  if ($direction < 0) {
    $start_looking = length($dna) if not defined $start_looking;

    while ((my $loc = rindex($dna, $flank, $start_looking)) >= 0) {
      push @matches, $loc + 1;
      $start_looking = $loc - 1;
    }
  } else {
    $start_looking = 0 if not defined $start_looking;

    while ((my $loc = index($dna, $flank, $start_looking)) >= 0) {
      push @matches, $loc + 1;
      $start_looking = $loc + 1;
    }
  }
  
  return sort { $a <=> $b } @matches;
}



##########################################

=head2 check_overlapping_CDS

  Title   :   check_overlapping_CDS
  Usage   :   my @genes = $mapper->check_overlapping_CDS("$seq", 1000, 2000);
  Returns :   array of CDSs or RNA genes falling in that region
  Args    :   any seq obj as string, coordinates relative to that seq obj

=cut

sub check_overlapping_CDS
  {
    my ($self, $seq, $x, $y) = @_;
    my $chromosome;

    # read in data from GFF files but ONLY if not already there ie do 1st time.
    unless ( defined $self->{'CHROM2GENE_POS'} ) {
      my $gff_dir = $self->{'DATABASE'}."/GFF_SPLITS";  # hard path not using Wormbase.pm
      croak "no GFF files in $gff_dir\n" unless (-e "$gff_dir");
      my @chromosomes = qw( I II III IV V X );
      my @methods = qw(curated miRNA snoRNA tRNA snRNA snlRNA rRNA non_coding_transcript scRNA asRNA piRNA lincRNA ncRNA miRNA_primary_transcript 7kncRNA);

      foreach my $chrom (@chromosomes) {
	foreach my $method (@methods) {	  
	  my $gff = "$gff_dir/CHROMOSOME_${chrom}_$method.gff";
	  open (GFF,"<$gff") or croak "cant open $gff\n";
	  my $chromosome = "CHROMOSOME_$chrom";

	  while (<GFF>) {
	    # CHROMOSOME_I curated CDS 222722  223159  . + . CDS "Y48G1BM.3" wp_acc=CE26120
	    # look for just protein or RNA genes by examining GFF_source and GFF_feature
	    next unless (/CDS/ or /primary_transcript/);
	    my @data = split;
	    if (($data[1] eq "curated"                  && $data[2] eq "CDS") ||
		($data[1] eq "Non_coding_transcript"    && $data[2] eq "nc_primary_transcript") ||
		($data[1] eq "miRNA_mature"             && $data[2] eq "miRNA") ||
		($data[1] eq "miRNA_precursor"          && $data[2] eq "pre_miRNA") ||
		($data[1] eq "snoRNA"                   && $data[2] eq "snoRNA") ||
		($data[1] eq "tRNA"                     && $data[2] eq "tRNA") ||
		($data[1] eq "snRNA"                    && $data[2] eq "snRNA") ||
		($data[1] eq "rRNA"                     && $data[2] eq "rRNA") ||
		($data[1] eq "ncRNA"                    && $data[2] eq "ncRNA") ||
		($data[1] eq "7kncRNA"                  && $data[2] eq "7kncRNA") ||
		($data[1] eq "asRNA"                    && $data[2] eq "asRNA") ||
 		($data[1] eq "lincRNA"                  && $data[2] eq "lincRNA") ||
 		($data[1] eq "piRNA"                    && $data[2] eq "piRNA") ||
		($data[1] eq "scRNA"                    && $data[2] eq "scRNA")) {
	      $data[9] =~ s/\"//g;
	      my $gene = $data[9];
	      my $end5 = $data[3];
	      my $end3 = $data[4];
	      $self->{'CHROM2GENE_POS'}->{"$chromosome"}->{"$gene"} = [$end5, $end3];
	    }
	  }
	  close GFF;
	}
      }
    }

    # determine if chromosome is being used as seq obj
    if( $self->{'single_chrom'}->{"$seq"} ) {
      $chromosome = $self->{'single_chrom'}->{"$seq"};
    }
    elsif( $seq =~ /CHROMOSOME_\d+/ ) {
      $chromosome = $seq;
    }
    else{
      ($chromosome,$x) = $self->Coords_2chrom_coords($seq,$x);
      ($chromosome,$y) = $self->Coords_2chrom_coords($seq,$y);
    }

    # if allele on -ve strand 3' coord will be bigger than 5'
    $self->swap(\$x, \$y) if( $y < $x ) ;

    my @gene_hits;
    foreach my $gene ( keys %{$self->{'CHROM2GENE_POS'}->{"$chromosome"}} ) {

      my $gene_x = $self->{'CHROM2GENE_POS'}->{"$chromosome"}->{$gene}->[0];
      my $gene_y = $self->{'CHROM2GENE_POS'}->{"$chromosome"}->{$gene}->[1];
      if (
	  ( $x > $gene_x and $x < $gene_y ) or # 5' end of allele is in gene
	  ( $y > $gene_x and $y < $gene_y ) or # 3' end of allele is in gene
	  ( $x < $gene_x and $y > $gene_y ) # whole gene removed
	 ) {
	
	push(@gene_hits,$gene);	
      }
    }
    return @gene_hits;
  }

##########################################

=head2 _get_flanking_sequence

  Title   :   _get_flanking_sequence
  Usage   :   my @flank_seq = $mapper->get_flanking_sequence("$seq", 1000, 2000);
  Returns :   array of two uppercase sequence strings of 30 bases or 
              more which uniquely define a region
    	      if a unique flanking sequence cannot be produced within the 
              bounds of this clone and the next larger sequence object 
              (superlink or chromosome) needs to be used, then 'undef' 
              is returned.
    	      This routine is the inverse of map_feature().
  Args    :   Sequence_object - any seq obj as string, 
              int, int - start/end coordinates of the region relative to that seq obj
	      int - (optional) the minimum length you want the flanking sequences to be, defaults to 30
              int - flag to not do a check for uniqueness
              int - flag to allow non-unique flanking sequences shorter than the minimum length at the ends of the sequence 

              (this last flag is intended for use in species with
              contig assemblies where Features are likely to be near
              the end of contigs and the contig sequence object
              supplied will be the top-level sequence object)

 IMPORTANT NOTE: The returned left flank will end at coord $pos1, and the returned right-flank will
                  start at $pos2 - in other words, the returned flanks will overlap the given extent
                  by 1bp on each side. Therefore, in order to get the flanking sequence for a feature with 
                  extent defined by $x - $y, this method should be called with $pos1 = $x-1 and $pos2 - $x+1. 
                  The reason for this hoopla is to be able to deal with 0-bp features (which can only be
                  defined with a 2bp extent)
=cut



sub _get_flanking_sequence {
  my ($self, $clone, $pos1, $pos2, $min_len, $no_unique_check, $short_end_flanks_allowed) = @_;

  # get sequence of clone
  my $seq = $self->Sub_sequence($clone);
  my $revseq = $self->DNA_revcomp($seq);
  my $len = length $seq;

  # convert to computer coords
  $pos1--;
  $pos2--;

  # are we in reverse sense? (i.e. reversed order of positions)
  if ($pos1 > $pos2) {
    ($seq, $revseq) = ($revseq, $seq);
    # get reverse coords
    $pos1 = $len-$pos1-1;
    $pos2 = $len-$pos2-1;
  }

  # set the initial flanking lengths
  if (! defined $min_len) {$min_len = 30}
  my $flank1 = $min_len;
  my $flank2 = $min_len;

  # loop to extend the sequence
  my $matches1 = 2;              # force at least one test of the flank by saying the last (imaginary) test found 2 matches
  my $matches2 = 2;
  my $flankseq1;
  my $flankseq2;
  while ($matches1 > 1 or $matches2 > 1) {

    # Check if the flanks are going off the ends of the sequence
    if ($pos1-$flank1 < 0) {
      if ($short_end_flanks_allowed) {
        $matches1=1; # declare that we have a unique match to allow exit from the while loop
        $flank1 = $pos1+1; # set the flank length to cover the sequence to the begining
	# the feature mapper finds the position of the longest flank first, so ensure that the other (unique) flank is longer than this one
	if ($flank2 < $flank1) {$flank2 = $flank1+1} 
      } else {
        return undef;
      }
    }
    if ($pos2+$flank2 >= $len) {
      if ($short_end_flanks_allowed) {
        $matches2=1; # declare that we have a unique match to allow exit from the while loop
        $flank2 = $len-$pos2; # set the flank length to cover the sequence to the end
	# the feature mapper finds the position of the longest flank first, so ensure that the other (unique) flank is longer than this one
	if ($flank1 < $flank2) {$flank1 = $flank2+1} 
      } else {
	return undef;
      }
    }

    # get flanking sequences
    $flankseq1 = substr($seq, $pos1-$flank1+1, $flank1);
    $flankseq2 = substr($seq, $pos2, $flank2);

    last if $no_unique_check;

    # find the number of matches. If we find a single match, good; but then also
    # need to check reverse strand because it may map there also
    if ($matches1 > 1) {
      $matches1 = $self->_matches($seq, $flankseq1);
      if ($matches1 == 1) {
	$matches1 += $self->_matches($revseq, $flankseq1);
      }
    }
    if ($matches2 > 1) {
      $matches2 = $self->_matches($seq, $flankseq2);
      if ($matches2 == 1) {
	$matches2 += $self->_matches($revseq, $flankseq2);
      }
    }

    # if there are more than one match, extend the length of the flank
    if ($matches1 > 1) {
      $flank1++;
    }
    if ($matches2 > 1) {
      $flank2++;
    }
  }

  # report the unique flanking sequences
  return (lc($flankseq1), lc($flankseq2));

}

=head2 get_flanking_sequence_for_feature

  Title   :   get_flanking_sequence_for_feature
  Usage   :   my @flank_seq = $mapper->get_flanking_sequence($seq, 100, 105);
  Returns :   array of two uppercase sequence strings of 30 bases or 
              more which uniquely define a region
    	      if a unique flanking sequence cannot be produced within the 
              bounds of this clone and the next larger sequence object 
              (superlink or chromosome) needs to be used, then 'undef' 
              is returned.
  Args    :   1. any seq obj as string, 
              2,3. start/end coordinates of the feature you wish to generate flanks for, relative to that seq obj
              4. whether or not given feature is 0-length
	      5. (optional) the minimum length you want the flanking sequences to be, defaults to 30
              6. (optional) flag to not do a check for uniqueness
              7. (optional) flag to allow non-unique flanking sequences shorter than the minimum length at the ends of the sequence 


  NOTE: when abs(end - start) == 1, the code cannot tell whether the supplied extent is
        a 2-bp feature or a 0-bp feature (lying in between the two basws). This information
        should therefore be supplied as the third arg
=cut

sub get_flanking_sequence_for_feature {
  my ($self, $clone, $start, $end, $is_zero_length, $min_flank_length, $no_unique_check, $short_end_flanks_allowed) = @_;

  if (abs($end - $start) != 1 or not $is_zero_length) {
    if ($start <= $end) {
      $start--;
      $end++;
    } else {
      $start++;
      $end--;
    }
  }

  return $self->_get_flanking_sequence($clone, $start, $end, $min_flank_length, $no_unique_check, $short_end_flanks_allowed);
  
}

##########################################
=head2 suggest_fix

  Title   :   suggest_fix
              used when map_feature() has failed. It suggests a possible fix for the mapping
  Usage   :   my @flank_seq = $mapper->suggest_fix($feature_id, $length, $seq, $flank_L, $flank_R, $version, @mapping_data);
  Returns :   $seq - clone or superlink
              suggested new $flank_L, $flank_R
              explanation for the suggested change
              1 = fixed, 0 = not fixed
  Args    :   Feature ID
              is_zero_length - 1 if the given feature is a 0-length feature (e.g. SL1, SL2), 0/undef otherwise
              any seq obj as string, e.g 'AC3'
              left and right flanking sequences that were tried
              version of database
              @mapping_data for remapping between currentdb and the current release

=cut

sub suggest_fix {

  my ($self, $feature_id, $expected_length, $clone, $flank_L, $flank_R, $mapper, $cache) = @_;
  my @result;
  my $FIXED = 1;
  my $NOT_FIXED = 0;

  # We can try to use the GFF file positions of features in the previous release
  # look for feature in currentdb
  # remap to current coordinates

  my $current_db = $self->wormbase->database('current');
  my $dir = "$current_db/CHROMOSOMES";

  if (not exists $cache->{$feature_id} or $cache->{$feature_id} == 0) {
    foreach my $chromosome ($self->wormbase->get_chromosome_names(-prefix => 1)) {
      my $gff = "$dir/${chromosome}.gff";
      open (GFF, "< $gff") || die "Can't open GFF file $gff\n";
      
      while (my $line = <GFF>) {
        chomp $line;
        if ($line =~ /^\s*$/) {next;}
        if ($line =~ /^#/) {next;}
        my @f = split /\t/, $line;
        if (!$f[8] || ($f[8] !~ /Feature \"\S+\"/)) {next;}
        
        my ($fid) = ($f[8] =~ /^Feature \"(\S+)\"/);
        if ($fid eq $feature_id or exists $cache->{$fid}) {
          $cache->{$fid} = [$f[0], $f[3], $f[4], $f[6]];
        }
      }
    }
  }
  
  if (ref($cache->{$feature_id})) {
    my ($chromosome, $start, $end, $sense) = @{$cache->{$feature_id}};

    my ($indel, $change);
    ($start, $end, $sense, $indel, $change) = $mapper->remap_gff($chromosome, 
                                                                 $start, 
                                                                 $end, 
                                                                 $sense); 
    
    # when getting new flanks, the returned flanks will include 1bp on each side of the 
    # extent that we supply. We therefore need to "extend" the feature by 1bp in each direction.
    
    # 0bp features appear in the GFF files as 2bp features. The GFF extent for these 
    # already includes 1bp of sequence. These therefore do not need adjustment. But how can 
    # we tell these from real 2bp features? The only option is to check the given min_len and
    # max_len
    if (abs($end - $start) != 1 or not defined $expected_length or $expected_length != 0) {
      $start--;
      $end++;
    }
    
    # get clone or superlink in the 4 Kb region around the feature
    my $left_coord = $start - 2000;
    my $right_coord = $end + 2000;
    my $new_clone;
    
    ($new_clone, $left_coord, $right_coord) = $self->LocateSpan($chromosome, $left_coord, $right_coord);
    $left_coord += 2000;
    $right_coord -= 2000;
    
    if ($sense eq '+') {
      ($flank_L, $flank_R) = $self->_get_flanking_sequence($new_clone, $left_coord, $right_coord);
    } else {
      ($flank_L, $flank_R) = $self->_get_flanking_sequence($new_clone, $right_coord, $left_coord);
    }
    if (defined $flank_L) {
      return ($new_clone, $flank_L, $flank_R, "New flanking sequences remapped from position in previous WormBase release", $FIXED);
    }
  }


  # If that didn't work, try to repair the existing flanking sequences

  my $dna = lc($self->Sub_sequence($clone));
  $flank_L = lc($flank_L);
  $flank_R = lc($flank_R);

  # check for non-ACGT characters in the flanking sequences
  if ($flank_L =~ s/\s//g || $flank_R =~ s/\s//g) {
    return ($clone,  $flank_L, $flank_R, "Space characters were found in the flanking sequences and corrected", $FIXED);
  } elsif ($flank_L =~ s/[^ACGT]//g || $flank_R =~ s/[^ACGT]//g) {
    return ($clone,  $flank_L, $flank_R, "Non-ACGT characters were found in the flanking sequences and corrected", $FIXED);
  }
  
  # finally, if one of the flanking sequences maps uniquely, and a given feature length
  # is given, suggest another flank
  if (defined $expected_length) {
    my ($flank_L_rev) = $self->DNA_revcomp($flank_L);
    my ($flank_R_rev) = $self->DNA_revcomp($flank_R);

    my @left_hits_for = $self->_check_for_match_left($dna, $flank_L, 1);
    my @left_hits_rev = $self->_check_for_match_right($dna, $flank_L_rev, -1);

    my @right_hits_for = $self->_check_for_match_right($dna, $flank_R, -1);
    my @right_hits_rev = $self->_check_for_match_left($dna, $flank_R_rev, 1);

    my $total_hits += scalar(@left_hits_for) + scalar(@left_hits_rev) + scalar(@right_hits_for) + scalar(@right_hits_rev);

    if ($total_hits == 1) {
      my ($left_c, $right_c);
      if (@left_hits_for) {
        ($left_c) = @left_hits_for;
        $right_c = $left_c + $expected_length + 1;
      } elsif (@right_hits_for) {
        ($right_c) = @right_hits_for;
        $left_c = $right_c - $expected_length - 1;
      } elsif (@left_hits_rev) {
        ($left_c) = @left_hits_rev;
        $right_c = $left_c - $expected_length - 1;
      } elsif (@right_hits_rev) {
        ($right_c) = @right_hits_rev;
        $left_c = $right_c + $expected_length + 1;
      }

      ($flank_L, $flank_R) = $self->_get_flanking_sequence($clone, $left_c, $right_c);
      
      if (defined $flank_L and defined $flank_R) {
        return ($clone, $flank_L, $flank_R, "Found unique match for one flank, extracted other based on position", $FIXED);
      }
    }

  }


  # can't do anything else
  return ($clone,  $flank_L, $flank_R, "Can't suggest a fix for this", $NOT_FIXED);

}


=head2 remap_and_generate_new_flanks_for_feature

  Title   :   remap_and_generate_new_flanks_for_feature
  Usage   :   my $hash_ref = $fmapper->remap_and_generate_new_flanks_for_features($assembly_mapper,
                 'Variation',
                 1,
                 \@feature_list);
  Returns :   Hash ref, indexed by feature/variation id, containing new flanks and remapped chromosomal
              location of the feature
  Args    :   1. Remap_Sequence_Change oject
              2. Class of features ('Variation' or 'Feature')
              3. whether or not given list contains 0-length features (which need coordinate adjustment)
	      4. Feature arrayreg

  NOTE: when abs(end - start) == 1, the code cannot tell whether the supplied extent is
        a 2-bp feature or a 0-bp feature (lying in between the two basws). This information
        should therefore be supplied as the third arg

  NOTE 2: Since the Remap_Sequence_Change only works sensible for C.elegans, so does this method
=cut

sub remap_and_generate_new_flanks_for_features {  
  my ($self, $assembly_mapper, $class, $zero_length, $min_flank_len, $feats) = @_;

  my $current_db = $self->wormbase->database('current');
  my $dir = "$current_db/CHROMOSOMES";

  my %feats;
  map { $feats{$_} = {} } @$feats;

  foreach my $chromosome ($self->wormbase->get_chromosome_names(-prefix => 1)) {
    my $gff = "$dir/${chromosome}.gff";
    open (GFF, "< $gff") || die "Can't open GFF file $gff\n";

    while (my $line = <GFF>) {
      next if $line =~ /^\#/;

      my @f = split /\t/, $line;
      next if $f[8] !~ /$class \"\S+\"/;

      my ($fid) = ($f[8] =~ /$class \"(\S+)\"/);
      next if not exists $feats{$fid};

      my ($chr, $start, $end, $sense) = ($f[0], $f[3], $f[4], $f[6]);

      my ($indel, $change);
      ($start, $end, $sense, $indel, $change) = $assembly_mapper->remap_gff($chr, 
                                                                            $start, 
                                                                            $end, 
                                                                            $sense); 


      if ($start > $end) {
        ($start, $end) = ($end, $start);
      }

      # when getting new flanks, the returned flanks will include 1bp on each side of the 
      # extent that we supply. We therefore need to "extend" the feature by 1bp in each direction.

      # but 0bp features appear in the GFF files as 2bp features. The GFF extent for these 
      # already includes 1bp of sequence. These therefore do not need adjustment. But how can 
      # we tell these from real 2bp features? The only option is to check the given min_len and
      # max_len
      if ($end - $start != 1 or not $zero_length) {
        $start--;
        $end++;
      }

      # we want the flanks to be unique within the chromosome, so we get them from the chr
      # and then map down
      my ($new_flank_L, $new_flank_R);

      if ($sense eq '+') {
	($new_flank_L, $new_flank_R) = $self->_get_flanking_sequence($chr, $start, $end, $min_flank_len);
      } else {
	($new_flank_L, $new_flank_R) = $self->_get_flanking_sequence($chr, $end, $start, $min_flank_len);
      }

      if (defined $new_flank_L) {
        $feats{$fid}->{left_flank} = uc($new_flank_L);
        $feats{$fid}->{right_flank} = uc($new_flank_R);
        $feats{$fid}->{chr}   = $chr;
        $feats{$fid}->{start} = $start;
        $feats{$fid}->{end}   = $end;
        $feats{$fid}->{strand}   = $sense;
      }
    }
    close (GFF);
  }

  foreach my $fid (@$feats) {
    if (not exists $feats{$fid}->{left_flank}) {
      delete $feats{$fid};
    }
  }

  return \%feats;
}
                          
##########################################

=head2 _approx_matches

  Title   :   _approx_matches
              finds the approximate matches of a string in the clone sequence
  Usage   :   my (@matches) = _approx_matches($sequence, $flanking_sequence)
  Returns :   array of (match start position, match end position, edit distance)
  Args    :   clone sequence string, flanking sequence string

=cut

sub _approx_matches () {
  my ($self, $seq, $flank) = @_;

  my @matches;

  # [i] - ignore case modifier
  # ($index, $size)   = String::Approx::aslice("pattern", [ modifiers ])

  $_ = $seq;
  my $offset=0;
  my ($pos, $size);
  while (($pos, $size) = aslice($flank, ["i"])) { # aslice() searches for $flank in $_
    $offset += $pos+$size;
    push @matches, [$pos, $offset-1];
    $_ = substr($seq, $offset);
  }

  return @matches;
}




##########################################
                                                                                                                                                      
=head2 _matches

  Title   :   _matches
  Usage   :   my $matches = _matches($sequence, $flanking_sequence)
  Returns :   list of positions of unique matches of the flanking_sequence in the sequence
  Args    :   clone sequence string, flanking sequence string

=cut

sub _matches () {
  my ($self, $seq, $flank) = @_;

  my @matches;

  my $pos = -1;
  while (($pos = index($seq, $flank, $pos)) > -1) {
    push @matches, $pos;
    $pos++;
  }
  return @matches;
}


sub wormbase {
  my ($self, $wb) = @_;

  if (defined $wb) {
    $self->{_wormbase} = $wb;
  }

  return $self->{_wormbase};
}

1;
