=pod 

=head1 NAME

 CDS

=head1 SYNOPSIS

 my $cds = CDS->new($name,\%exons,"+","I",$transformer);
 my $match = $cds->map_cDNA( $cdna );
 $cds->add_matching_cDNA ($cdna ) if ( $match == 1 );
 $cds->add_3_UTR( $cdna2 );

 $cds->gene_start( 123405 );
 my $gene_span = $cds->gene_start - $cds->gene_end;

 $cds->transcripts( $new_transcript );
 my $transcripts = $cds->transcripts;

 $cds->report( *FH, $coords, $transformer );

=head1 DESCRIPTION

 This object represents a CDS in the transcript building process.  It controls the attaching of matching cDNAs and generation of new transcript objects when required.  Much of how it works in done through passing things on to its child transcript objects ( see Transcript.pm )

Inherits from SequenceObj ( SequenceObj.pm )

=head1 CONTACT

Anthony  ar2@sanger.ac.uk


=head1 METHODS

=cut


package CDS;

use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'} ;
use Carp;
use Modules::SequenceObj;
use Modules::Transcript;
use strict;

our @ISA = qw( SequenceObj );

=head2 new

    Title   :   new
    Usage   :   CDS->new($name,\%exons,"+","I",$transformer);
    Function:   Creates new CDS object
    Returns :   ref to self
    Args    :   name - string
                hash ref of exon structure
                strand as string
                chromosome as string
                ref to transformer obj

=cut

sub new
  {
    my $class = shift;
    my $name = shift;
    my $exon_data = shift;
    my $strand = shift;
    my $chromosome = shift;
    my $transformer = shift;

    my $self = SequenceObj->new($name, $exon_data, $strand);
    bless ( $self, $class );

    $self->transform_strand($transformer,"transform") if ( $self->strand eq "-" );

    $self->transformer( $transformer );
    my $transcript = Transcript->new( $name, $self);
    $self->transcripts( $transcript );

    $self->gene_start( $transcript->start );
    $self->gene_end( $transcript->end );

    if ($chromosome) {
      $transcript->chromosome( $chromosome ) if $chromosome;
      $self->chromosome( $chromosome ) if $chromosome;
    }

    return $self;
  }

# $CDS->map_cDNA results in calls to Transcript->map_cDNA for each transcript. One of which will be derived from the initial CDS structure.

=head2 map_cDNA

    Title   :   map_cDNA
    Usage   :   $cds->map_cDNA( $cdna )
    Function:   check if the passed sequence object matches the exon structure if itself.  If it matches the cds but not any existing transcript then a new transcript will be constructed
    Returns :   1 if match 0 otherwise
    Args    :   sequence_object


=cut

sub map_cDNA
  {
    my $self = shift;
    my $cdna = shift;

    # check strandedness
    if( $SequenceObj::debug ) { # class data
      ( print $self->name," has no strand\n" and return 0 ) unless $self->strand;
      ( print $cdna->name," has no strand\n" and return 0 ) unless $cdna->strand;
    }
    return 0 if $self->strand ne $cdna->strand;

    my $matches_me = 0;
    foreach my $transcript ( $self->transcripts ) {
      if ($transcript->map_cDNA( $cdna ) == 1) { 
	$matches_me = 1;
      }
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
	return 0 unless ($self->check_features($cdna) == 1);
	#this must overlap - check exon matching
	if( $self->check_exon_match( $cdna ) ) {
	  # check reciprocal CDS -> cdna
	  if( $cdna->check_exon_match( $self )) {
	    # if this cdna matches the CDS but not the existing transcripts create a new one
	    # append .x to indicate multiple transcripts for same CDS.
	    my $transcript_count = scalar($self->transcripts);
	    my $new_name;
	    if( $transcript_count == 1 ) {
	      # rename the original as .1
	      $new_name = $self->name . ".$transcript_count";
	      my @transcripts = $self->transcripts;
	      $transcripts[0]->name("$new_name");
	    }
	    $transcript_count++;
	    $new_name = $self->name . ".$transcript_count";
	    
	    my $transcript = Transcript->new( $new_name, $self);
	    $transcript->chromosome( $self->chromosome );

	    # this will add the new cDNA through the correct method, ensuring that exons are extended accordingly.
	    $transcript->map_cDNA($cdna);
	    # now recheck cDNAs already matched to CDS to new transcript
#	    my $matched_cDNAs = $self->matching_cDNAs;
#	    foreach my $match ( @{$matched_cDNAs} ) {
#	      $transcript->map_cDNA($match);
#	    }

	    # add new transcript to CDS obj
	    $self->transcripts($transcript);
	    $matches_me = 1;
	  }
	}
      }
    }
    $self->add_matching_cDNA($cdna) if $matches_me == 1;
    return $matches_me;
  }

=head2 transcripts

    Title   :   transcripts
    Usage   :   $cds->transcripts
    Function:   Add to / query existing transcripts for this cds object
    Returns :   array of refs to transcript obj
    Args    :   none or new transcript obj

=cut

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

=head2 add_matching_cDNA

    Title   :   add_matching_cDNA
    Usage   :   $cds->add_matching_cDNA( $cdna )
    Function:   add matching_cDNA to list
    Returns :   nothing
    Args    :   cdna object

=cut

sub add_matching_cDNA
  {
    my $self = shift;
    my $cdna = shift;
    #print STDERR $cdna->name," matches ",$self->name,"\n";
    push( @{$self->{'matching_cdna'}},$cdna);
  }
=head2 

    Title   :   add_3_UTR
    Usage   :   $cds->add_3_UTR( $cdna )
    Function:   Adds cdna identified by paired read matches ( passes them on to transcript objs )
    Returns :   nothing
    Args    :   cdna object

=cut

sub add_3_UTR
  {
    my $self = shift;
    my $cdna = shift;
    foreach my $transcript ( $self->transcripts ) {
      $transcript->add_3_UTR( $cdna );
    }
    $self->add_matching_cDNA( $cdna );
  }

=head2 report

    Title   :   report
    Usage   :   $cds->report
    Function:   print out relevant data to passed file hadle.  Calls report on transcript objs
    Returns :   nothing
    Args    :   filehandle to print to
                Coords_convert object ref
                Transformer object ref

=cut

sub report
  {
    my $self = shift;
    my $fh = shift;
    my $coords = shift;
    my $transformer = shift;

    #$fh = STDOUT unless defined $fh;

    print $fh "\nCDS : \"",$self->name,"\"\n";
    foreach (@{$self->{'matching_cdna'}}) {
      print $fh "Matching_cDNA \"",$_->name,"\"Inferred_Automatically \"transcript_builder.pl\"\n";
    }

    foreach  ( $self->transcripts ) {
      print $fh "Corresponding_transcript \"",$_->name,"\"\n";
    }

    foreach (@{$self->{'matching_cdna'}}) {
      print $fh "\nSequence :\"",$_->name,"\n";
      print $fh "Matching_CDS ",$self->name," Inferred_Automatically \"transcript_builder.pl\"\n";
    }

    foreach ( $self->transcripts ) {
      $_->report($fh, $coords);
    }
  }

=head2 gene_start

    Title   :   gene_start
    Usage   :   $cds->gene_start
    Function:   sets / returns start of the gene - ie 5 prime most coord of all transcripts
    Returns :   int
    Args    :   none or candidate new start coord

=cut

sub gene_start
  {
    my $self = shift;
    my $start = shift;

    if( $start ) {
      if( $self->{'gene_start'} ) {
	$self->{'gene_start'} = $start if ( $start < $self->{'gene_start'} ) ;
      }
      else {
	$self->{'gene_start'} = $start;
      }
    }
    return $self->{'gene_start'};
  }

=head2 gene_end

    Title   :   gene_end 
    Usage   :   $cds->gene_end
    Function:   sets / returns end of the gene - ie 3 prime most coord of all transcript
    Returns :   int
    Args    :   none or candidate new end coord

=cut

sub gene_end
  {
    my $self = shift;
    my $end = shift;

    if( $end )  {
      if( $self->{'gene_end'} ) { 
	$self->{'gene_end'} = $end if $end > $self->{'gene_end'};
      }
      else {
	$self->{'gene_end'} = $end;
      }
    }
    return $self->{'gene_end'};
  }

1;
