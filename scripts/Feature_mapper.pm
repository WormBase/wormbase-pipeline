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
use Modules::Remap_Sequence_Change;

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

sub new
  {
    my ($class,$database,$refresh,$wormbase) = @_;
    my $self = Sequence_extract->invoke($database, $refresh,$wormbase);

    bless $self, $class;
    return $self;
}

##########################################

=head2 map_feature

  Title   :   map_feature
  Usage   :   my @map_info = $mapper->map_feature("AH6","actgtacgtagcgagcaccgatcaggacgag","tgactagcggacagcgagcagctagctgat");
  Returns :   smallest sequence object that contains the given flanking seq and coords
  Args    :   sequence obj as string, 2 flanking seqs as strings

=cut


sub map_feature {
  my ($self, $seq, $flank_L, $flank_R) = @_;
  
  my ($left_coord, $right_coord);
  
  # convert to chromosome coords
  my ($chr_id, $chr_st, $chr_en) = $self->LocateSpanUp($seq);
  # for reverse-oriented assembly components, coords come back in reverse order
  if ($chr_en < $chr_st) {
    ($chr_st, $chr_en) = ($chr_en, $chr_st);
  }

  my ($max_extension, $extend_by) = (5, 1000);
  
  for(my $cur_extension = 0; $cur_extension < $max_extension; $cur_extension++) {
    $chr_st -= ($cur_extension * $extend_by);
    $chr_en += ($cur_extension * $extend_by);

    $chr_st = 1 if $chr_st < 1;
    $chr_en = $self->Superlink_length($chr_id) if $chr_en > $self->Superlink_length($chr_id);

    my $reg_len = ($chr_en - $chr_st + 1);

    my $dna = $self->Sub_sequence($chr_id, $chr_st - 1, $reg_len);
    my ($dna_left_coord, $dna_right_coord) = $self->_check_for_match($dna,$flank_L, $flank_R);
    
    if (defined $dna_left_coord and defined $dna_right_coord) {
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
  
  print "Cant map to $seq\n";
  return 0;
}

##########################################

=head2 _check_for_match

 NOTE: This function takes raw DNA seq as a string ie "aagtcaatcggatatgtgatggagctagctgatcgatcgtgc . . . . . etc ";

  Title   :   _check_for_match
  Usage   :   $self->_check_for_match($dna,"actgtacgtagcgagcaccgatcaggacgag","tgactagcggacagcgagcagctagctgat");
  Returns :   coords of feature being mapped ie 1 base after the left flank ends and 1 before right flank starts

                        flank_L                                                     flank_R
              actgtacgtagcgagcaccgatcaggacgag-----------------------------tgactagcggacagcgagcagctagctgat
                                             ^                           ^ 
  Args    :   dna as string, 2 flanking seqs as strings


                        flank_L             vv       flank_R
              actgtacgtagcgagcaccgatcaggacgagTGACTACGTGCTATGCAGCGAGCAT
                                            RL
=cut

sub _check_for_match 
  {
    my ($self, $dna, $flank_L, $flank_R) = @_;
    my ($rev_left,$rev_right,$offset);
    my ($match_left,$match_right,$span);

    # make sure all in same case !
    $dna = lc($dna);
    $flank_L = lc($flank_L);
    $flank_R = lc($flank_R);

    my $dna_length = length $dna;

    # check forward strand
    if ($dna =~ /$flank_L/i) {
      $match_left = length ($`);
    }
    if ($dna =~ /$flank_R/i) {
      $match_right = length ($`) + 1;
    }

    if( $match_left and $match_right ) {
      $match_left += (length $flank_L);
    }

    # check rev strand
    else  { 
      $rev_left     = $self->DNA_revcomp($flank_L);
      $rev_right    = $self->DNA_revcomp($flank_R);

      if ($dna =~ /$rev_left/) {
	$offset = length ($`);
	$match_left = $offset +1 ;
	
	#only try the right if the left is success
	if ($dna =~ /$rev_right/) {
	  $offset = length ($`);
	  $match_right = $offset + (length $flank_R);
	}
      }
    }
    
    if($match_left and $match_right ){ 
      return ($match_left,$match_right);
    }
    else {
      return undef;
    }
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
      my @methods = qw(curated miRNA snoRNA tRNAscan-SE-1.23 tRNAscan-SE-1.3 snRNA rRNA non_coding_transcript scRNA ncRNA miRNA_primary_transcript);

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
		($data[1] eq "miRNA_mature_transcript"  && $data[2] eq "miRNA") ||
		($data[1] eq "curated"                  && $data[2] eq "miRNA_primary_transcript") ||
		($data[1] eq "snoRNA_mature_transcript" && $data[2] eq "snoRNA") ||
		($data[1] eq "tRNAscan-SE-1.23"         && $data[2] eq "tRNA") ||
		($data[1] eq "snRNA_mature_transcript"  && $data[2] eq "snRNA") ||
		($data[1] eq "rRNA"                     && $data[2] eq "rRNA_primary_transcript") ||
		($data[1] eq "Non_coding_transcript"    && $data[2] eq "nc_primary_transcript") ||
		($data[1] eq "ncRNA"                    && $data[2] eq "ncRNA_primary_transcript") ||
		($data[1] eq "scRNA_mature_transcript"  && $data[2] eq "scRNA")) {
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

=head2 get_flanking_sequence

  Title   :   get_flanking_sequence
  Usage   :   my @flank_seq = $mapper->get_flanking_sequence("$seq", 1000, 2000);
  Returns :   array of two uppercase sequence strings of 30 bases or 
              more which uniquely define a region
    	      if a unique flanking sequence cannot be produced within the 
              bounds of this clone and the next larger sequence object 
              (superlink or chromosome) needs to be used, then 'undef' 
              is returned.
    	      This routine is the inverse of map_feature().
  Args    :   any seq obj as string, 
              start/end coordinates of the region relative to that seq obj
	      (optional) the minimum length you want the flanking sequences to be, defaults to 30

=cut

sub get_flanking_sequence
{

  my ($self, $clone, $pos1, $pos2, $min_len) = @_;


  # get sequence of clone
  my $seq = $self->Sub_sequence($clone);
  my $len = length $seq;
                                                                                                                                                      
  # convert to computer coords
  $pos1--;
  $pos2--;

  # are we in reverse sense? (i.e. reversed order of positions)
  if ($pos1 > $pos2) {
    $seq = $self->DNA_revcomp($seq);
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
  while ($matches1 > 1 || $matches2 > 1) {
    # Can't get unique first flank in sequence object $clone ending at $pos1
    if ($pos1-$flank1+1 < 0) {return undef;}
    # Can't get unique second flank in sequence object $clone starting at $pos2
    if ($pos2+$flank2 > $len) {return undef;}
    # get flanking sequences
    $flankseq1 = substr($seq, $pos1-$flank1+1, $flank1);
    $flankseq2 = substr($seq, $pos2, $flank2);
                                                                                                                                                      
    # find the number of matches
    $matches1 = $self->_matches($seq, $flankseq1);
    $matches2 = $self->_matches($seq, $flankseq2);
                                                                                                                                                      
    #print uc($flankseq1) . " ($matches1) " . uc($flankseq2) ." ($matches2)\n";
                                                                                                                                                      
    # if there are more than one match, extend the length of the flank
    if ($matches1 > 1) {
      $flank1++;
    }
    if ($matches2 > 1) {
      $flank2++;
    }
  }
                                                                                                                                                      
  # report the unique flanking sequences
  return (uc($flankseq1), uc($flankseq2));
                                                                                                                                                      
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
              expected length of the feature, or -1 if the type of feature does not have a single fixed length
              any seq obj as string, e.g 'AC3'
              left and right flanking sequences that were tried
              version of database
              @mapping_data for remapping between currentdb and the current release

=cut

sub suggest_fix
{

  my ($self, $feature_id, $length, $clone, $flank_L, $flank_R, $version, @mapping_data) = @_;
  my @result;
  my $FIXED = 1;
  my $NOT_FIXED = 0;

  # We can try to use the GFF file positions of features in the previous release
  # look for feature in currentdb
  # remap to current coordinates

  my $dir = "/nfs/wormpub/DATABASES/current_DB/CHROMOSOMES";
  foreach my $chromosome qw(I II III IV V X) {
    my $gff = "$dir/CHROMOSOME_${chromosome}.gff";
    open (GFF, "< $gff") || die "Can't open GFF file $gff\n";

    while (my $line = <GFF>) {
      chomp $line;
      if ($line =~ /^\s*$/) {next;}
      if ($line =~ /^#/) {next;}
      my @f = split /\t/, $line;
      if (!$f[8] || ($f[8] !~ /Feature \"$feature_id\"/)) {next;}
      my ($chromosome, $start, $end, $sense) = ($f[0], $f[3], $f[4], $f[6]);
      my ($indel, $change);
      ($start, $end, $sense, $indel, $change) = Remap_Sequence_Change::remap_gff($chromosome, $start, $end, $sense, $version - 1, $version, @mapping_data);
      # get clone or superlink in the 4 Kb region around the feature
      my $left_coord = $start - 2000;
      my $right_coord = $end + 2000;
      my $new_clone;
      # shift the coords in zero-based coordinates
      ($new_clone, $left_coord, $right_coord) = $self->LocateSpan($chromosome, $left_coord-1, $right_coord-1);
      $left_coord += 2000;
      $right_coord -= 2000;
      $left_coord++; # convert back to human coords
      $right_coord++;
      if ($sense eq '+') {
	($flank_L, $flank_R) = $self->get_flanking_sequence($new_clone, $left_coord, $right_coord);
      } else {
	($flank_L, $flank_R) = $self->get_flanking_sequence($new_clone, $right_coord, $left_coord);
      }
      if (defined $flank_L) {
	return ($new_clone, $flank_L, $flank_R, "New flanking sequences remapped from position in previous WormBase release", $FIXED);
      }
      last;
    }
    close (GFF);
  }


  # If that didn't work, try to repair the existing flanking sequences

  my $dna = uc $self->Sub_sequence($clone);
  $flank_L = uc $flank_L;
  $flank_R = uc $flank_R;


  # check for non-ACGT characters in the flanking sequences
  if ($flank_L =~ s/\s//g || $flank_R =~ s/\s//g) {
    @result = ($clone,  $flank_L, $flank_R, "Space characters were found in the flanking sequences and corrected", $FIXED);
    return;
  } elsif ($flank_L =~ s/[^ACGT]//g || $flank_R =~ s/[^ACGT]//g) {
    @result = ($clone,  $flank_L, $flank_R, "Non-ACGT characters were found in the flanking sequences and corrected", $FIXED);
    return;
  }
  

  # which flank has an exact match in which orientation?
  my @matchesR = $self->_matches($dna, $flank_R);
  my @matchesL = $self->_matches($dna, $flank_L);

  # check rev strand
  my $rev_flank_R = $self->DNA_revcomp($flank_R);
  my $rev_flank_L = $self->DNA_revcomp($flank_L);
  my @rev_matchesR = $self->_matches($dna, $rev_flank_R);
  my @rev_matchesL = $self->_matches($dna, $rev_flank_L);

  # two missing flanking sequences: nothing can be done
  if (!@matchesR && !@matchesL && !@rev_matchesR && !@rev_matchesL) {
    @result = ($clone,  $flank_L, $flank_R, "Neither of the flanking sequences map - nothing can be done with this", $NOT_FIXED);
    return @result;
  }


  # more than one match on both flanks
  if (scalar @matchesR > 1 && scalar @matchesL > 1 || scalar @rev_matchesR > 1 && scalar @rev_matchesL > 1 ) {
    @result = ($clone,  $flank_L, $flank_R, "Both flanks match more than once - nothing can be done with this", $NOT_FIXED);
    return @result;
  }


  # more than one match on one flank
  # this needs work to check each of the multiple hits to see which is closest to the flank with one hit
  if (scalar @matchesR == 1 && scalar @matchesL > 1 || scalar @rev_matchesR == 1 && scalar @rev_matchesL > 1 ) {
    @result = ($clone,  $flank_L, $flank_R, "Left flank matches more than once - nothing can be done with this", $NOT_FIXED);
    return @result;

  } elsif (scalar @matchesL == 1 && scalar @matchesR > 1 || scalar @rev_matchesL == 1 && scalar @rev_matchesR > 1) {
    @result = ($clone,  $flank_L, $flank_R, "Right flank matches more than once - nothing can be done with this", $NOT_FIXED);
    return @result;
  }

  # is the feature a fixed size?
  if ($length ne -1) {
    my ($start, $end);
    if (scalar @matchesR == 1) {
      $end = $matchesR[0];
      $start = $end - $length -1;
      ($flank_L, $flank_R) = $self->get_flanking_sequence($clone, $start+1, $end+1); # convert pos to human coords
      return ($clone,  $flank_L, $flank_R, "Left flank sequence has been updated based on the right flank position", $FIXED);

    } elsif (scalar @matchesL == 1) {
      $start = $matchesL[0] + length ($flank_L) - 1;
      $end = $start + $length +1;
      ($flank_L, $flank_R) = $self->get_flanking_sequence($clone, $start+1, $end+1); # convert pos to human coords
      return ($clone,  $flank_L, $flank_R, "Right flank sequence has been updated based on the left flank position", $FIXED);

      
    } elsif (scalar @rev_matchesR == 1) {
      $end = $rev_matchesR[0] + length ($flank_R) - 1;
      $start = $end + $length +1;
      ($flank_L, $flank_R) = $self->get_flanking_sequence($clone, $start+1, $end+1); # convert pos to human coords
      return ($clone,  $flank_L, $flank_R, "Left flank sequence has been updated based on the right flank position (reverse sense)", $FIXED);

    } elsif (scalar @rev_matchesL == 1) {
      $start = $rev_matchesL[0];
      $end = $start - $length -1;
      ($flank_L, $flank_R) = $self->get_flanking_sequence($clone, $start+1, $end+1); # convert pos to human coords
      return ($clone,  $flank_L, $flank_R, "Right flank sequence has been updated based on the left flank position (reverse sense)", $FIXED);

    }
  }

  # can't do anything else
  return ($clone,  $flank_L, $flank_R, "Can't suggest a fix for this", $NOT_FIXED);

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
    print "match = pos: $pos to $offset, size: $size\n";
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




1;
