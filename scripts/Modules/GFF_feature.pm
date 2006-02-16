package GFF_feature;

use lib $ENV{'CVS_DIR'} ;
use Carp;

sub new
  {
    my( $class, $name, $start, $end, $score, $strand, $f_start, $f_end ) = @_;

    my $self = {};
    $self->{'name'}   = $name   if $name;
    $self->{'start'}  = $start  if $start;
    $self->{'end'}    = $end    if $end;
    $self->{'score'}  = $score  if $score;
    $self->{'strand'} = $strand if $strand;
    $self->{'f_start'}= $f_start if $f_start;
    $self->{'f_end'}  = $f_end  if $f_end;

    bless ( $self, $class );
    return $self;
  }

sub start 
  {
    my $self = shift;
    my $start = shift;
    $self->{'start'} = $start if $start;
    return $self->{'start'};
  }

sub end 
  {
    my $self = shift;
    my $end = shift;
    $self->{'end'} = $end if $end;
    return $self->{'end'};
  }

sub f_start 
  {
    my $self = shift;
    my $f_start = shift;
    $self->{'f_start'} = $f_start if $f_start;
    return $self->{'f_start'};
  }

sub f_end 
  {
    my $self = shift;
    my $f_end = shift;
    $self->{'f_end'} = $f_end if $f_end;
    return $self->{'f_end'};
  }

sub name
  {
    my $self = shift;
    my $name = shift;
    $self->{'name'} = $name if $name;
    return $self->{'name'};
  }

sub strand
  {
    my $self   = shift;
    my $strand = shift;
    $self->{'strand'} = $strand if $strand;
    return $self->{'strand'};
  }

sub score
  {
    my $self   = shift;
    my $score = shift;
    $self->{'score'} = $score if $score;
    return $self->{'score'};
  }

sub print 
  {
    my $self = shift;
    print STDERR $self->name," ",$self->start," ",$self->end,"\n";
  }

sub chromosome
  {
    my $self   = shift;
    my $chromosome = shift;
    $self->{'chromosome'} = $chromosome if $chromosome;
    return $self->{'chromosome'};
  }

sub method
  {
    my $self   = shift;
    my $method = shift;
    $self->{'method'} = $method if $method;
    return $self->{'method'};
  }

sub feature
  {
    my $self   = shift;
    my $feature = shift;
    $self->{'feature'} = $feature if $feature;
    return $self->{'feature'};
  }

1;
