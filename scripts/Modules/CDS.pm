package CDS;

use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'} ;
use Carp;
use Modules::SequenceObj;
use Modules::Transcript;

@ISA = qw( SequenceObj );

sub new
  {
    my $class = shift;
    my $name = shift;
    my $exon_data = shift;
    my $strand = shift;

    my $self = SequenceObj->new($name, $exon_data, $strand);
    my $transcript = Transcript->new( $name, $exons_data, $strand);

    push (@{$self->{'transcripts'}},$transcript);

    bless ( $self, $class );
    return $self;
  }

sub map_cDNA
  {
    my $self = shift;
    my $cdna = shift;

    # check strandedness
    return 0 if $self->strand ne $cdna->strand;

    #this must overlap - check exon matching
    foreach $transcript ( $self->transcripts ) {
      $transcripts->check_exon_match( $cdna );
    }
      return 1;
  }

sub transcripts
  {
    my $self = shift;
    my $transcript = shift;

    # if a new one is passed add
    if( $transcript ) {
      push (@{$self->{'transcripts'}},$transcript);
    }

    return @{$self->{'transcripts'};
  }

1;
