package Strand_transformer;

use Carp;
use strict;

sub new 
  {
    my $class = shift;
    my $self = {};
    $self->{'start'} = shift;
    $self->{'end'} = shift;

    bless( $self, $class );
    return $self;
  }

# this transforms neg strand coords to be as if fwd
sub transform_neg_coord
  {
    my $self = shift;
    my $coord = shift;
    die "sub transform_strand_coord requires a coord\n" unless $coord;
    return ( $coord * -1 + $self->{'end'} + $self->{'start'} );
  }

# revert to original strandedness
sub revert_to_neg
  {
    my $self = shift;
    my $coord = shift;
    die "sub transform_strand_coord requires a coord\n" unless $coord;
    return ( (($coord - $self->{'start'}) - $self->{'end'})/-1) ;
  }

1;
