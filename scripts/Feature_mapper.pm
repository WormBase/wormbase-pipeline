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

use Sequence_extract;
use Carp;

@ISA = ('Sequence_extract');

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


=head2 map_feature

  Title   :   map_feature
  Usage   :   my @map_info = $mapper->map_feature("AH6","actgtacgtagcgagcaccgatcaggacgag","tgactagcggacagcgagcagctagctgat");
  Returns :   smallest sequence object that contains the given flanking seq and coords
  Args    :   sequence obj as string, 2 flanking seqs as strings

=cut


sub map_feature
  {
    my ($self, $seq, $flank_L, $flank_R) = @_;
    my $dna = $self->Sub_sequence($seq);
    my ($left_coord, $right_coord) = $self->_check_for_match($dna,$flank_L, $flank_R);

  NOT_MAPPED:
    while( ! (defined $right_coord) ) {
      my $extension = 1000;
      while( $extension < 5001 ) {
	$dna = $self->Sub_sequence($seq,"-$extension","+$extension");
	($left_coord, $right_coord) = $self->_check_for_match($dna,$flank_L, $flank_R);

	
	if ( $left_coord and $right_coord) {
	  # make coordinate adjsutments for -+ extension
	  # at this stage they'll be something like AH6 (-+2000), 1500, 37500)

	  $left_coord -= $extension;
	  $right_coord -= $extension;
	  ($seq,$left_coord, $right_coord) = $self->LocateSpan($seq,$left_coord,$right_coord);
	  last NOT_MAPPED ;
	}
	$extension += 1000;
      }
      print "cant map to $seq +- $extension\n";
      return 0;
    }

    return ($seq,$left_coord, $right_coord);
  }


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
      my @methods = qw(curated miRNA snoRNA tRNAscan-SE-1.23 snRNA rRNA non_coding_transcript scRNA ncRNA);

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
	    if (($data[1] eq "curated"               && $data[2] eq "CDS") ||
		($data[1] eq "miRNA"                 && $data[2] eq "miRNA_primary_transcript") ||
		($data[1] eq "snoRNA"                && $data[2] eq "snoRNA_primary_transcript") ||
		($data[1] eq "tRNAscan-SE-1.23"      && $data[2] eq "tRNA_primary_transcript") ||
		($data[1] eq "snRNA"                 && $data[2] eq "snRNA_primary_transcript") ||
		($data[1] eq "rRNA"                  && $data[2] eq "rRNA_primary_transcript") ||
		($data[1] eq "Non_coding_transcript" && $data[2] eq "nc_primary_transcript") ||
		($data[1] eq "ncRNA"                 && $data[2] eq "ncRNA_primary_transcript") ||
		($data[1] eq "scRNA"                 && $data[2] eq "scRNA_primary_transcript")) {
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
              coordinates of the region relative to that seq obj

=cut

sub get_flanking_sequence
{

  my ($self, $clone, $pos1, $pos2) = @_;


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
  my $flank1 = 30;
  my $flank2 = 30;
                                                                                                                                                      
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
                                                                                                                                                      
=head2 _matches

  Title   :   _matches
  Usage   :   my $matches = _matches($sequence, $flanking_sequence)
  Returns :   the number of unique matches of the flanking_sequence in the sequence
  Args    :   any seq obj as string, coordinates relative to that seq obj

=cut

sub _matches () {
  my ($self, $seq, $flank) = @_;
                                                                                                                                                      
  my $matches = 0;
                                                                                                                                                      
  my $pos = -1;
  while (($pos = index($seq, $flank, $pos)) > -1) {
    $matches++;
    $pos++;
  }
                                                                                                                                                      
  return $matches;
}




1;
