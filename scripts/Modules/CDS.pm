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
    my $chromosome = shift;

    my $self = SequenceObj->new($name, $exon_data, $strand);

    my $transcript = Transcript->new( $name, $exon_data, $strand);
    $transcript->chromosome( $chromosome ) if $chromosome;
    push (@{$self->{'transcripts'}},$transcript);
    
    bless ( $self, $class );
    $self->chromosome( $chromosome ) if $chromosome;

    return $self;
  }

# $CDS->map_cDNA results in calls to Transcript->map_cDNA for each transcript. One of which will be derived from the initial CDS structure.

sub map_cDNA
  {
    my $self = shift;
    my $cdna = shift;

    # check strandedness
    return 0 if $self->strand ne $cdna->strand;
    return 0 if $self->strand ne "+";

    #this must overlap - check exon matching

    my $matches_me = 0;
    foreach $transcript ( $self->transcripts ) {
      $matches_me = 1 if ($transcript->map_cDNA( $cdna ) == 1);
    }

    if( $matches_me == 0 ) {
      # check against just CDS structure  
      if( $self->start > $cdna->end ) {
	return 0;
      }
      elsif( $cdna->start > $self->end ) {
	return 0;
      }
      else {
	#this must overlap - check exon matching
	if( $self->check_exon_match( $cdna ) ) {
	  # if this cdna matches the CDS but not the existing transcripts create a new one
	  my $transcript = Transcript->new( $self->name, $self->exon_data, $self->strand);
	  $transcript->add_matching_cDNA($cdna);
	  
	  $self->transcripts($transcript);
	  $matches_me = 1;
	}
      }
    }
    $self->add_matching_cDNA($cdna) if $matches_me == 1;
    return $matches_me;
  }

sub transcripts
  {
    my $self = shift;
    my $transcript = shift;

    # if a new one is passed add
    if( $transcript ) {
      push (@{$self->{'transcripts'}},$transcript);
    }
    return @{$self->{'transcripts'}};
  }

sub add_matching_cDNA
  {
    my $self = shift;
    my $cdna = shift;
    print STDERR $cdna->name," matches ",$self->name,"\n";
    push( @{$self->{'matching_cdna'}},$cdna->name);
  }

sub report
  {
    my $self = shift;
    my $fh = shift;
    my $coords = shift;

    print $fh "\nCDS : \"",$self->name,"\"\n";
    foreach (@{$self->{'matching_cdna'}}) {
      print $fh "Matching_cDNA \"$_\"\n";
    }

    foreach  ( $self->transcripts ) {
      print $fh "Corresponding_transcript \"",$_->name,"\"\n";
    }

    foreach ( $self->transcripts ) {
      $_->report($fh, $coords);
    }
  }

1;
