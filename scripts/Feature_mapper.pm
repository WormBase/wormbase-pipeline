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

use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"
  : glob("~ar2/wormbase/scripts");

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
    my $class = shift;
    my $database = shift;
    my $refresh = shift;

    my $self = Sequence_extract->invoke($database, $refresh);

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

=cut

sub _check_for_match
  {
    my ($self, $dna, $flank_L, $flank_R) = @_;
    my ($rev_left,$rev_right,$offset);
    my ($match_left,$match_right,$span);

    my $dna_length = length $dna;

    # check forward strand
    if ($dna =~ /$flank_L/i) {
      $match_left = length ($`);
    }
    if ($dna =~ /$flank_R/i) {
      $match_right = length ($`);
    }

    if( $match_left and $match_right ) {
      $self->swap(\$match_left, \$match_right) if( $match_left > $match_right );

      $match_left += (length $flank_L) + 1;
    }

    # check rev strand
    else  { 
      $rev_left     = $self->DNA_revcomp($flank_L);
      $rev_right    = $self->DNA_revcomp($flank_R);

      if ($dna =~ /$rev_left/) {
	$offset = length ($`);
	$match_left = $offset;
	
	#only try the right if the left is success
	if ($dna =~ /$rev_right/) {
	  $offset = length ($`);
	  $match_right = $offset + (length $flank_R) + 1 ;
	
	  $self->swap(\$match_left, \$match_right);
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
      croak "no GFF files in $self->{DATABASE}/CHROMOSOMES\n" unless (-e "$self->{DATABASE}/CHROMOSOMES");
      my @chromosomes = qw( I II III IV V X );

      foreach (@chromosomes) {
	my $gff = "$self->{'DATABASE'}/CHROMOSOMES/CHROMOSOME_$_.gff";
	open (GFF,"<$gff") or croak "cant open $gff\n";
	my $chromosome = "CHROMOSOME_$_";

	while (<GFF>) {
	  # CHROMOSOME_I curated CDS 222722  223159  . + . CDS "Y48G1BM.3" wp_acc=CE26120
	  if ( ( /curated/ and /CDS/) or (/RNA/ and /Transcript/) ) {
	    my @data = split;
	    next unless( "$data[2]" eq "CDS" ); #using full gff so only need 'curated' 'CDS'?
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


    ($chromosome,$x) = $self->Coords_2chrom_coords($seq,$x);
    ($chromosome,$y) = $self->Coords_2chrom_coords($seq,$y);

    # if allele of > 1bp on - strand 3' coord will be bigger than 5'
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





1;
