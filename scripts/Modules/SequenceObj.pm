=pod 

=head1 NAME

 SequenceObj

=head1 SYNOPSIS

 my $seq = SequenceObj->new($name,\%exons,"+");
 $seq->sort_exons
 $seq->exon_data
 $seq->transform_strand($transformer,'transform');
 $seq->check_exon_match( $cdna );

 my $exons = $seq->( exon_data ); 
 foreach (keys %{$exons} ) {
   print "$_ $exons->{$_}\n";
 }

 my $first_exon = $seq->first_exon; # array ref
 my $last_exon  = $seq->last_exon; # array ref

 $self->_cDNA_wholelyInExon( $cdna );
 $self->_exon_that_ends( $exon_end );

 my $name   = $seq->name;
 my $start  = $seq->start;
 my $end    = $seq->end;
 my $strand = $seq->strand;

 $seq->mapped( $cds ); sets the CDS to which the mRNA is attached

 # feature data is added / queried via specific name method that in turn all call 'feature' with the specific type set
 # no distinction is made between SL1 and SL2. Although there are methods for both they both end up as SL
 my $SL = $seq->( SL1 );
 $seq->polyA_site( [ 182772,  182773,  "WBsf01634" ] ); array of coords and name is passed


=head1 DESCRIPTION

 This object represents a SequenceObj for use in the transcript_builder.pl script.  It is a generic object that stores exon structure as a hash 

 'exons' => HASH(0x143160ad0)
            4900023 => 4900309
            4899865 => 4899943
            4900377 => 4900592


and a sorted array of arrays 

 'sorted_exons' => ARRAY(0x143160b80)
      0  ARRAY(0x14315ffe8)
         0  4899865
         1  4899943
      1  ARRAY(0x14315fe88)
         0  4900023
         1  4900309
      2  ARRAY(0x14315fec8)
         0  4900377
         1  4900592

Also stores coordinate info - start, end  and strand. and features associated with the sequence

Inherited by CDS.pm Transcript.pm 

=head1 CONTACT

Anthony  ar2@sanger.ac.uk


=head1 METHODS

=cut


package SequenceObj;

use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'} ;
use Carp;
use strict;

# new 
# expects name , %-ref of exon start/ends and strand ( +/- )

our $debug; # class variable - available to all instances of SequenceObj and classes that inherit from it.

=head2 new

    Title   :   new
    Usage   :   SequenceObj->new($name,\%exons,"+");
    Function:   Creates new SequenceObj object
    Returns :   ref to self
    Args    :   name - string
                hash ref of exon structure
                strand as string

=cut

sub new
  {
    my $class = shift;
    my $name = shift;
    my $exon_data = shift; # \%
    my $strand = shift;

    my $self = {};
    if ($name) {
      $self->{'name'}   = $name;

      my ($start, $end);
      my @tmp;
      foreach ( keys %{$exon_data} ) {
	if ( !(defined $start) or $start > $_ ) {
	  $start = $_;
	}
	if ( !(defined $end) or $end < $$exon_data{$_} ) {
	  $end = $$exon_data{$_};
	}
	$self->{'exons'}->{$_} = $$exon_data{$_};
	push(@tmp,[($_,$$exon_data{$_})]);
      }
      @{$self->{'sorted_exons'}} = sort { $a->[0] <=> $b->[0] } @tmp;
      $self->{'start'} = $start;
      $self->{'end'}   = $end;
      $self->{'strand'}= $strand ;
    }
    bless ( $self, $class );
    return $self;
  }

=head2 sort_exons

    Title   :   sort_exons
    Usage   :   $seq->sort_exons
    Function:   create sorted array of exon arrays
                 [0]->( 1200, 1250 )
                 [1]->( 1300, 1350 )
                 [2]->( 1400, 1500 ) etc
    Returns :   nothing
    Args    :   none
               

=cut

sub sort_exons
  {
    my $self = shift;
    my ($start, $end);
    my $exon_data = $self->exon_data;

    my @tmp;
    foreach ( keys %{$exon_data} ) {
      if( !(defined $start) or $start > $_ ) {
	$start = $_;
      }
      if( !(defined $end) or $end < $$exon_data{$_} ) {
	$end = $$exon_data{$_};
      }
      $self->{'exons'}->{$_} = $$exon_data{$_};
      push(@tmp,[($_,$$exon_data{$_})]);
    }
    @{$self->{'sorted_exons'}} = sort { $a->[0] <=> $b->[0] } @tmp;
    $self->start( $start );
    $self->end( $end );
  }

=head2 transform_strand

    Title   :   transform_strand
    Usage   :   $seq->transform_strand( $transformer,'transform')
                $seq->transform_strand( $transformer,'revert')
    Function:   convert ( and revert ) negative strand coords to pseudo fwd so that same code exon comparison will work
    Returns :   nothing
    Args    :   Strand_transformer
                direction ( 'transform' or 'revert' )

=cut

sub transform_strand
  {
    my $self = shift;
    my $transformer = shift;
    my $direction = shift;

    my %tmp_exons;
    foreach ( keys %{$self->exon_data} ){
      my ($key,$value);
      if( $direction eq "transform" ){
	# swap exon start / end too
	$value = $transformer->transform_neg_coord( $_);
	$key   = $transformer->transform_neg_coord( $self->exon_data->{$_});
      }
      elsif ( $direction eq "revert" ) {
	# swap exon start / end too
	$value = $transformer->revert_to_neg( $_);
	$key   = $transformer->revert_to_neg( $self->exon_data->{$_});
      }
      else { die "need a transformation direction\n";	   }
	
      $tmp_exons{$key} = $value;
    }
    $self->exon_data(\%tmp_exons);

    $self->sort_exons;

    # transform feature data ( SL1 etc ).
    foreach my $feature ( keys %{ $self->{'feature'}} ) {
      $self->{'feature'}->{"$feature"} = [( $transformer->transform_neg_coord( $self->{'feature'}->{"$feature"}->[1]),
					    $transformer->transform_neg_coord( $self->{'feature'}->{"$feature"}->[0]),
					    $self->{'feature'}->{"$feature"}->[2]
					  )];
    }
  }

=head2 check_exon_match

    Title   :   check_exon_match
    Usage   :   $seq->check_exon_match( $cdna )
    Function:   check that the passed SequenceObj derived object has valid matching exon structure to this.  This comparison is mainly for comparing cDNAs to CDSs so has criteria such as 'final exon of gene so allow extension past end (UTR)'.  Also used for the reverse check cDNA ->cds
             It writes a flag in the exon data of the passed object indicating the match type of each exon, for use in Transcript.pm
    Returns :   1 for match; 0 for fail
    Args    :   SequenceObj
               

=cut

sub check_exon_match 
  {
    my $self = shift;
    my $cdna = shift;

    #use match_code to record match type eg exact exon match to make using this cdna in later steps easy

    #check if cDNA exon fits with gene model
    foreach my $exon ( @{$cdna->sorted_exons}) {
      my $cExonStart = $exon->[0];
      my $gExonS;
      
      # do cDNA and gene share exon start position
      if ( $self->{'exons'}->{"$cExonStart"} ) {
	if ($self->{'exons'}->{"$cExonStart"} == $cdna->{'exons'}->{"$cExonStart"} ) {
	  #exact match
	  print "\tExact Match\n" if $debug;
	  $exon->[2] = 1;
	}
	#is this final gene exon
	elsif ( $cExonStart == $self->last_exon->[0] ) {
	  if( $cdna->exon_data->{"$cExonStart"} > $self->last_exon->[1] ) {
	    print "\tMatch - last SeqObj exon\n" if $debug ;
	    $exon->[2] = 2;
	  }
	  elsif ( $cdna->exon_data->{"$cExonStart"} == $cdna->end ) {
	    print "Match - cDNA final exon ends in final gene exon\n"  if $debug;
	    $exon->[2] = 7;
	  }
	  else {
	    print STDERR "MISS : cDNA splices in last SeqObj exon\n" if $debug ;
	    return 0;
	  }
	}
	# or final cDNA exon?
	elsif ( $cExonStart == $cdna->last_exon->[0] ) {
	  # . . must terminate within gene exon
	  if ( $cdna->{'exons'}->{"$cExonStart"} > $self->{'exons'}->{"$cExonStart"} ) {
	    print STDERR "\tMISS - ",$cdna->name," $cExonStart => ",$cdna->{'exons'}->{$cExonStart}," extends over gene exon boundary\n" if $debug ;
	    return 0;
	  } else {
	    print "\tMatch - last cDNA exon\n" if $debug ;
	    $exon->[2] = 3;
	  }
	}
	else {
	  print STDERR "\tMISS -  ",$cdna->name," $cExonStart => ",$cdna->{'exons'}->{$cExonStart}," extends over gene exon boundary\n"  if $debug;
	  return 0;
	}
      }
      # do cDNA and gene share exon end position
      elsif ( ( $gExonS = $self->_exon_that_ends( $cdna->{'exons'}->{"$cExonStart"} ) and ($gExonS != 0) ) ) {
	#	# shared exon end
	
	if ( $gExonS == $self->first_exon->[0] ) { #is this the 1st gene exon 
	  if ( $cExonStart == $cdna->first_exon->[0] ) { # also cDNA start so always match
	    print "\tMatch - 1st exons end in same place\n" if $debug ;
	    $exon->[2] = 4;
	  }
	  elsif ( $cExonStart < $self->first_exon->[0] ) { # cDNA exon overlap 1st gene exon
	    print "\tMatch - cDNA exon covers 1st gene exon\n" if $debug ;
	    $exon->[2] = 5;
	    # extends 5'
	  }
	  else {
	    print STDERR "\tMISS - cDNA exon splices in gene exon\n" if $debug ;
	    print STDERR "\t\t",$cdna->name," $cExonStart => ",$cdna->exon_data->{$cExonStart},"\n" if $debug ;
	    print STDERR "\t\t",$self->name," $gExonS => ",$self->exon_data->{$gExonS},"\n"  if $debug;
	    return 0;
	  }
	}
	# exon matched is not 1st of SeqObj
	elsif ( ($cExonStart == $cdna->first_exon->[0] ) and # start of cDNA
		($cExonStart >$gExonS ) ) { # . . . is in SeqObj exon
	  print"\tMatch - 1st exon of cDNA starts in exon of SeqObj\n" if $debug ;
	  $exon->[2] = 6;
	} 
	else {
	  print STDERR "MISS - exon ",$cdna->name," : $cExonStart => ",$cdna->exon_data->{$cExonStart}," overlaps start of gene exon : $gExonS => ",$self->exon_data->{$gExonS},"\n"  if $debug;
	  return 0;
	}
      }# cDNA_wholelyInExon
      elsif ( $self->_cDNA_wholelyInExon($cdna) ) {
	print "Match cDNA contained in exon\n"  if $debug;
	$exon->[2] = 7;
      }
      # single exon gene contained in cDNA
      elsif( ( $cExonStart < $self->first_exon->[0] ) and
	     ( $cdna->exon_data->{$cExonStart} > $self->last_exon->[1] ) and
	     ( scalar keys %{$self->exon_data} == 1 ) # single exon gene
	   ) {
	print "Match single exon gene contained in cDNA\n"  if $debug;
	$exon->[2] = 13;
      }
      # cDNA exon overlaps gene 1st exon start and terminate therein
      elsif( ( $cExonStart == $cdna->last_exon->[0] ) and #  last exon of cDNA
	     ( $cExonStart < $self->first_exon->[0] ) and 
	     ( $cdna->last_exon->[1] > $self->first_exon->[0] and $self->first_exon->[0] <$self->first_exon->[1] )
	   ) {
	print "\tcDNA final exon overlaps first exon of gene and end therein\n" if $debug ;
	$exon->[2] = 8;
      }
      # cDNA exon starts in final gene exon and continues past end
      elsif( ($cdna->start > $self->last_exon->[0]) and 
	     ($cdna->start < $self->last_exon->[1]) and 
	     ($cdna->first_exon->[1] > $self->last_exon->[1] )
	   ) {
	print "MATCH : final cDNA exon starts in final gene exon and continues past end\n"  if $debug;
	$exon->[2] = 9;
      }
      # exon lies outside of CDS ( but other parts of cDNA overlap it )
      elsif( $exon->[1] < $self->start ) {
	print "MATCH : 5\'UTR exon\n" if $debug ;
	$exon->[2] = 10;
      }
      elsif( $exon->[0] > $self->end ) {
	print "MATCH : 3\'UTR exon\n"  if $debug;
	$exon->[2] = 11;
      }
      else {
	# doesnt match
	print STDERR $cdna->name," doesnt match ",$self->name,"\n" if $debug;
	return 0;
      }
    }
    return 1;
  }

=head2 _cDNA_wholelyInExon

    Title   :   _cDNA_wholelyInExon
    Usage   :   $self_cDNA_wholelyInExon
    Function:   internal method to determine if a passed SequenceObj lies completely within an exon of this SequenceObj
    Returns :   1 if does lie within exon; 0 otherwise
    Args    :   SequenceObj
               

=cut

sub _cDNA_wholelyInExon
  {
    my $self = shift;
    my $cdna = shift;

    foreach ( keys %{$self->exon_data} ) {
      if ( $cdna->start > $_ and $cdna->end < $self->{'exons'}->{$_} ) {
	return 1;
      }
    }
    return 0;
  }


=head2 _exon_that_end

    Title   :  _exon_that_end
    Usage   :  $self->_exon_that_end ( 20000 )
    Function:  internal method to find it this SequenceObj has an exon that ends with passed coord 
    Returns :   exon start coord or 0
    Args    :   coordinate as int
               

=cut

# expects coord of exon end to compare. Returns start of exon if match , or else 0 
sub _exon_that_ends
  {
    my $self = shift;
    my $exon_end = shift;
    foreach (keys %{$self->exon_data} ) {
      return $_ if( $self->exon_data->{$_} == $exon_end );
    }
    return 0;
  }

=head2 exon_data

    Title   :  exon_data
    Usage   :  exon_data
    Function:  $seq->exon_data( %exons )  
    Returns :   hash ref of exons
    Args    :   hash or nothing
               

=cut

sub exon_data
  {
    my $self = shift;
    my $new_exons = shift;
    $self->{'exons'} = $new_exons if $new_exons;
    return $self->{'exons'};
  }

=head2 first_exon

    Title   :   first_exon
    Usage   :   $seq->first_exon
    Returns :   returns array of first exon coords ( 1250, 1300 )
    Args    :   none
               

=cut

sub first_exon 
  {
    my $self = shift;
    return $self->{'sorted_exons'}->[0];
  }

=head2 last_exon

    Title   :  last_exon
    Usage   :  $seq->last_exon
    Returns :  returns array of last exon coords ( 1450, 1600 )
    Args    :  none
               

=cut

sub last_exon 
  {
    my $self = shift;
    return $self->{'sorted_exons'}->[-1];
  }

=head2 start

    Title   :  start
    Usage   :  $seq->start( '1250' )
    Function:  get / set start of object
    Returns :  start coord as int
    Args    :  coord as int or none
               

=cut

sub start 
  {
    my $self = shift;
    my $start = shift;
    $self->{'start'} = $start if $start;
    return $self->{'start'};
  }

=head2 end

    Title   :  end
    Usage   :  $seq->end( '1600' )
    Function:  get / set end of object
    Returns :  end coord as int
    Args    :  coord as int or none
               

=cut

sub end 
  {
    my $self = shift;
    my $end = shift;
    $self->{'end'} = $end if $end;
    return $self->{'end'};
  }

=head2 name

    Title   :   name
    Usage   :  $seq->name( "F45G2.2.1" )
    Function:   get / set name
    Returns :   name as string
    Args    :   name as string or none
               

=cut

sub name
  {
    my $self = shift;
    my $name = shift;
    $self->{'name'} = $name if $name;
    return $self->{'name'};
  }

=head2 strand

    Function:  get / set strand assignment
    Returns / Args : "+" or "-" as string
               

=cut

sub strand
  {
    my $self = shift;
    return $self->{'strand'};
  }

=head2 sorted_exons

    Title   :  sorted_exons
    Usage   :  my $second_exon = $seq->sorted_exons->[1]
               my $second_exon_start = $second_exon->[0]
    Function:  get sorted exons 
    Returns :  sorted array of exon arrays ( sees synopsis )
    Args    :  none
               

=cut

sub sorted_exons
  {
    my $self = shift;
    return $self->{'sorted_exons'};
  }

=head2 mapped

    Function:  get / set function for assigning Transcript / CDS that this is assiged to ( cDNA method really  !)
    Returns / Args :  CDS object

=cut

#this is really a cDNA specific method but cant be bothered to create new class yet !
sub mapped
  {
    my $self = shift;
    my $CDS = shift;
    $self->{'CDS'} = $CDS if( $CDS );
    
    return $self->{'CDS'};
  }

=head2 chromosome

    Function:   get / set chromosome assignment
    Returns / Args : eg "I" as string

=cut

sub chromosome
  {
    my $self = shift;
    my $chromosome = shift;
    $self->{'chromosome'} = $chromosome if $chromosome;
    return $self->{'chromosome'};
  }

=head2 transformer

    Title   :   transformer ( see Strand_transformer.pm )
    Usage   :   $seq->transformer( $transformer )
    Function:   get / set Strand_transformer 
    Returns :   Strand_transformer ref
    Args    :   
               

=cut

sub transformer 
  {
    my $self = shift;
    my $transformer = shift;
    $self->{'transformer'} = $transformer if $transformer;

    return $self->{'transformer'};
  }

=head2 debug 

    Title   :   debug
    Usage   :   $seq->debug
    Function:   set debug class variable
    Returns :   $debug value
    Args    :   defined of none
               

=cut

sub debug
  {
    my $self = shift;
    my $set = shift;
    $debug = 1 if $set;
    return $debug;
  }

=head2 features

    Title   :  features
    Usage   :  $seq->feature
    Function:  get array of features associated with this SequenceObj 
    Returns :  array of features
    Args    :  none
               

=cut

# just returns the WBfeature names
sub features
  {
    my $self = shift;
    my @features;
    foreach my $feature ( keys %{$self->{feature}} )  {
      push( @features, $self->{feature}->{$feature}->[2] );
    }
    return @features; 
  }

=head2 feature

    Title   :  feature 
    Usage   :  $self->feature( 'SL',\@coords );
    Function:  add / query specific feature type
    Returns :  array ( coord coord WBsf_id )
    Args    :  type ('SL', polyA_site etc ) and array ( coord coord WBsf_id )
               

=cut

# add / query specific feature type
sub feature
  {
    my $self = shift;
    my $feature = shift;
    my $data = shift; #@  182772  182773  WBsf01634

    # self=>feature=>SL   =>( x y id)
    #              =>polyA_site=>( x y id)

    if( $data ) { # adding new feature
      $self->{'feature'}->{ "$feature" } = $data;
    }
    else { 
      return $self->{'feature'}->{$feature};
    }
  }

=head2 SL1

    Usage   :   $seq->feature('SL',[$x, $y, $WBsf_id] );
                my $sl = $seq->SL1;
    Function:   get / set feature data
  
               

=cut
  
sub SL1
  {
    my ($self, $data) = @_;
    my $type = "SL";
    return $self->feature( $type, $data );
  }

=head2 SL2

   see SL1

=cut

sub SL2
  {
    my ($self, $data) = @_;
    my $type = "SL";
    return $self->feature( $type, $data );
  }


=head2 SL

   see SL1

=cut

sub SL
  {
    my ($self, $data) = @_;
    my $type = "SL";
    return $self->feature( $type, $data );
  }

=head2 polyA_site

   see SL1

=cut

sub polyA_site
  {
    my ($self, $data) = @_;
    my $type = "polyA_site";
    return $self->feature( $type, $data );
  }

=head2 polyA_signal_sequence

   see SL1

=cut

sub polyA_signal_sequence
  {
    my ($self, $data) = @_;
    my $type = "polyA_signal";
    return $self->feature( $type, $data );
  }

=head2 polyA_signal

   see SL1

=cut

# added this so i dont have to change loads of code where i've written it as this ;)
sub polyA_signal
  {
    my ($self, $data) = @_;
    my $type = "polyA_signal";
    return $self->feature( $type, $data );
  }

=head2 array_index

    Title   :   array_index
    Function:   get / set array_index - used to store position in ordered arrays in transcript_builder.pl

=cut

sub array_index
  {
    my $self = shift;
    my $index = shift;
    $self->{'index'} = $index if defined $index;
    return $self->{'index'};
  }

=head2 matching_cDNAs

    Returns :   array of SequenceObj s

=cut

sub matching_cDNAs
  {
    my $self = shift;
    return $self->{'matching_cdna'};
  }

=head2 paired_read

    Function:   get / set paired read for ESTs
    Returns :   SequenceObj
    Args    :   SequenceObj
               

=cut

sub paired_read
  {
    my $self = shift;
    my $pair = shift;

    $self->{'paired_read'} = $pair if $pair;
    return $self->{'paired_read'};
  }


1;
