package Feature_mapper;

use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts"
  : glob("~ar2/wormbase/scripts");

use Sequence_extract;

@ISA = ('Sequence_extract');

sub new
  {
    my $class = shift;
    my $database = shift;
    my $refresh = shift;

    my $self = Sequence_extract->invoke($database, $refresh);

    bless $self, $class;
    print "now doing $class (Feat_map) constructor\n\n\n";

    return $self;
}

sub map_feature
  {
    my ($self, $seq, $flank_L, $flank_R) = @_;
    my $dna = $self->Sub_sequence($seq);

    my ($rev_left,$rev_right,$offset);
    my ($match_left,$match_right,$span);

    my $dna_length = length $dna;

    print "length of $seq dna is $dna_length\n";
    # scan forward strand
    if ($dna =~ /$flank_L/) {
      $match_left = length ($`);
    }
    if ($dna =~ /$flank_R/) {
      $match_right = length ($`);
    }

    if( $match_left and $match_right ) {
      $self->swap(\$match_left, \$match_right) if( $match_left > $match_right );

      $match_left += (length $flank_L) + 1;
    }

    # match to rev strand
    else  { 
      $rev_left     = $self->DNA_revcomp($flank_L);
      $rev_right    = $self->DNA_revcomp($flank_R);

      if ($dna =~ /$rev_left/) {
	$offset = length ($`);
	$match_left = $dna_length - $offset + 1;
      }
      if ($dna =~ /$rev_right/) {
	$offset = length ($`);
	$match_right = $dna_length - $offset - length $flank_R;
      }
    }

    return ($seq,$match_left,$match_right);
}

1;
