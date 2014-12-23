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

use lib $ENV{'CVS_DIR'} ;
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

sub new {
  my ($class, 
      $name, 
      $exon_hash, 
      $strand) = @_;

  my $self = {};
  bless $self, $class;

  if ($name) {
    $self->exon_data($exon_hash);
    $self->sort_exons;

    $self->{name} = $name;
    $self->{strand} = $strand;
    $self->{probably_matching_cds} = [];

    my $prev_exon_end;
    my (%introns, @introns);
    
    foreach my $exon ( @{$self->sorted_exons}) {
      my $exon_start = $exon->[0];
      if (defined $prev_exon_end) {
        push @introns, [$prev_exon_end + 1, $exon_start - 1];
          $introns{$prev_exon_end + 1} = $exon_start - 1;
      }
      $prev_exon_end = $exon->[1];
    }
    $self->intron_data(\%introns);
    $self->{sorted_introns} = \@introns;
  }
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

sub sort_exons {
  my $self = shift;
  my ($start, $end);
  my $exon_data = $self->exon_data;
  
  my (%new_exons, @new_exons);

  foreach ( keys %{$exon_data} ) {
    if(not defined $start or $start > $_ ) {
      $start = $_;
    }
    if(not defined $end or $end < $exon_data->{$_} ) {
      $end = $exon_data->{$_};
    }
    $new_exons{$_} = $exon_data->{$_};
    push @new_exons, [$_,$exon_data->{$_}];
  }

  $self->{sorted_exons} = [sort { $a->[0] <=> $b->[0] } @new_exons];
  $self->exon_data(\%new_exons);

  $self->start( $start );
  $self->end( $end );
}


sub sort_introns {
  my $self = shift;

  my $intron_data = $self->intron_data;
  
  my (%new_introns, @new_introns);

  foreach ( keys %{$intron_data} ) {
    $new_introns{$_} = $intron_data->{$_};

    $self->{introns}->{$_} = $intron_data->{$_};
    push @new_introns, [$_,$intron_data->{$_}];
  }
  $self->{sorted_introns} = [sort { $a->[0] <=> $b->[0] } @new_introns];
  $self->intron_data(\%new_introns); 
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

sub transform_strand {
  my $self = shift;
  my $transformer = shift;
  my $direction = shift;
  
  my (%tmp_exons, %tmp_introns);

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
    else { 
      die "need a transformation direction\n";	   
    }
    
    $tmp_exons{$key} = $value;
  }
  $self->exon_data(\%tmp_exons); 
  $self->sort_exons;

  foreach ( keys %{$self->intron_data} ){
    my ($key,$value);
    if( $direction eq "transform" ){
      # swap exon start / end too
      $value = $transformer->transform_neg_coord($_);
      $key   = $transformer->transform_neg_coord($self->intron_data->{$_});
    } elsif ( $direction eq "revert" ) {
      # swap exon start / end too
      $value = $transformer->revert_to_neg($_);
      $key   = $transformer->revert_to_neg($self->intron_data->{$_});
    } else { 
      die "need a transformation direction\n";	   
    }
    
    $tmp_introns{$key} = $value;
  }
  $self->intron_data(\%tmp_introns); 
  $self->sort_introns;

  # transform feature data ( SL1 etc ).
  foreach my $ftype ( keys %{ $self->{feature}} ) {
    if ($direction eq "transform") {
      $self->{feature}->{$ftype} = [( $transformer->transform_neg_coord( $self->{feature}->{$ftype}->[1]),
                                      $transformer->transform_neg_coord( $self->{feature}->{$ftype}->[0]),
                                      $self->{feature}->{$ftype}->[2]
                                      )];
    } elsif ($direction eq "revert") {
      $self->{feature}->{$ftype} = [( $transformer->revert_to_neg( $self->{feature}->{$ftype}->[1]),
                                      $transformer->revert_to_neg( $self->{feature}->{$ftype}->[0]),
                                      $self->{feature}->{$ftype}->[2]
                                      )];
    }
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
	if ($cdna->{'exons'}->{"$cExonStart"} and
            $self->{'exons'}->{"$cExonStart"} == $cdna->{'exons'}->{"$cExonStart"} ) {
	  #exact match
	  print STDERR "SequenceObj::check_exon_match\tExact Match\n" if $debug;
	  $exon->[2] = 1;
	}
	#is this final gene exon
	elsif ( $cExonStart == $self->last_exon->[0] ) {
	  if( $cdna->exon_data->{"$cExonStart"} > $self->last_exon->[1] ) {
	    print STDERR "SequenceObj::check_exon_match\tMatch - last SeqObj exon\n" if $debug ;
	    $exon->[2] = 2;
	  }
	  elsif ( $cdna->exon_data->{"$cExonStart"} == $cdna->end ) {
	    print STDERR "SequenceObj::check_exon_match\tMatch - cDNA final exon ends in final gene exon\n"  if $debug;
	    $exon->[2] = 7;
	  }
	  else {
	    print STDERR "SequenceObj::check_exon_match\tMISS : cDNA splices in last SeqObj exon\n" if $debug ;
	    return 0;
	  }
	}
	# or final cDNA exon?
	elsif ( $cExonStart == $cdna->last_exon->[0] ) {
	  # . . must terminate within gene exon
	  if ( $cdna->{'exons'}->{"$cExonStart"} > $self->{'exons'}->{"$cExonStart"} ) {
	    print STDERR "SequenceObj::check_exon_match\tMISS - ",$cdna->name," $cExonStart => ",$cdna->{'exons'}->{$cExonStart}," extends over gene exon boundary\n" if $debug ;
	    return 0;
	  } else {
	    print STDERR "SequenceObj::check_exon_match\tMatch - last cDNA exon\n" if $debug ;
	    $exon->[2] = 3;
	  }
	}
	else {
	  print STDERR "SequenceObj::check_exon_match\tMISS -  ",$cdna->name," $cExonStart => ",$cdna->{'exons'}->{$cExonStart}," extends over gene exon boundary\n"  if $debug;
	  return 0;
	}
      }
      # do cDNA and gene share exon end position
      elsif ( ( $gExonS = $self->_exon_that_ends( $cdna->{'exons'}->{"$cExonStart"} ) and ($gExonS != 0) ) ) {
	#	# shared exon end
	
	if ( $gExonS == $self->first_exon->[0] ) { #is this the 1st gene exon 
	  if ( $cExonStart == $cdna->first_exon->[0] ) { # also cDNA start so always match
	    print STDERR "SequenceObj::check_exon_match\tMatch - 1st exons end in same place\n" if $debug ;
	    $exon->[2] = 4;
	  }
	  elsif ( $cExonStart < $self->first_exon->[0] ) { # cDNA exon overlap 1st gene exon
	    print STDERR "SequenceObj::check_exon_match\tMatch - cDNA exon covers 1st gene exon\n" if $debug ;
	    $exon->[2] = 5;
	    # extends 5'
	  }
	  else {
	    print STDERR "SequenceObj::check_exon_match\tMISS - cDNA exon splices in gene exon\n" if $debug ;
	    print STDERR "SequenceObj::check_exon_match\t\t",$cdna->name," $cExonStart => ",$cdna->exon_data->{$cExonStart},"\n" if $debug ;
	    print STDERR "SequenceObj::check_exon_match\t\t",$self->name," $gExonS => ",$self->exon_data->{$gExonS},"\n"  if $debug;
	    return 0;
	  }
	}
	# exon matched is not 1st of SeqObj
	elsif ( ($cExonStart == $cdna->first_exon->[0] ) and # start of cDNA
		($cExonStart >$gExonS ) ) { # . . . is in SeqObj exon
	  print STDERR "SequenceObj::check_exon_match\tMatch - 1st exon of cDNA starts in exon of SeqObj\n" if $debug ;
	  $exon->[2] = 6;
	} 
	else {
	  print STDERR "SequenceObj::check_exon_match\tMISS - exon ",$cdna->name," : $cExonStart => ",$cdna->exon_data->{$cExonStart}," overlaps start of gene exon : $gExonS => ",$self->exon_data->{$gExonS},"\n"  if $debug;
	  return 0;
	}
      }# cDNA_wholelyInExon
      elsif ( $self->_cDNA_wholelyInExon($cdna) ) {
	print STDERR "SequenceObj::check_exon_match\tMatch cDNA contained in exon\n"  if $debug;
	$exon->[2] = 7;
      }
      # single exon gene contained in cDNA
      elsif( ( $cExonStart < $self->first_exon->[0] ) and
	     ( $cdna->exon_data->{$cExonStart} > $self->last_exon->[1] ) and
	     ( scalar keys %{$self->exon_data} == 1 ) # single exon gene
	   ) {
	print STDERR "SequenceObj::check_exon_match\tMatch single exon gene contained in cDNA\n"  if $debug;
	$exon->[2] = 13;
      }
      # cDNA exon overlaps gene 1st exon start and terminate therein
      elsif( ( $cExonStart == $cdna->last_exon->[0] ) and #  last exon of cDNA
	     ( $cExonStart < $self->first_exon->[0] ) and 
             ( $cdna->last_exon->[1] > $self->first_exon->[0] and $cdna->last_exon->[1] <$self->first_exon->[1] )
	   ) {
	print STDERR "SequenceObj::check_exon_match\tcDNA final exon overlaps first exon of gene and ends therein\n" if $debug ;
	$exon->[2] = 8;
      }
      # cDNA exon starts in final gene exon and continues past end
      elsif( ($cdna->start > $self->last_exon->[0]) and 
	     ($cdna->start < $self->last_exon->[1]) and 
	     ($cdna->first_exon->[1] > $self->last_exon->[1] )
	   ) {
	# cdna has intron 3' of CDS end so final exon doesn't overlap. Needs to be added as exon rather than just extending.
	if( $cExonStart > $self->end ) {
	  # UTR exon
	  print STDERR "SequenceObj::check_exon_match\tMATCH : final cDNA exon overlaps gene end and cDNA has further splicing\n" if $debug;
	  $exon->[2] = 14;
	}
	else{
	  print STDERR "SequenceObj::check_exon_match\tMATCH : final cDNA exon starts in final gene exon and continues past end\n"  if $debug;
	  $exon->[2] = 9;
	}
      }
      # exon lies outside of CDS ( but other parts of cDNA overlap it )
      elsif( $exon->[1] < $self->start ) {
	print STDERR "SequenceObj::check_exon_match\tMATCH : 5\'UTR exon\n" if $debug ;
	$exon->[2] = 10;
      }
      elsif( $exon->[0] > $self->end ) {
	print STDERR "SequenceObj::check_exon_match\tMATCH : 3\'UTR exon\n"  if $debug;
	$exon->[2] = 11;
      }
      else {
	# doesnt match
	print STDERR "SequenceObj::check_exon_match\t" . $cdna->name . " doesnt match " . $self->name ."\n" if $debug;
	return 0;
      }
    }
    return 1;
  }

=head2 check_intron_match

    Title   :   check_intron_match
    Usage   :   $seq->check_intron_match( $cdna )
    Function:   count the number of consecutive introns that the cdna and SequenceObj structures have in common
    Returns :   number of consecutive matched introns
    Args    :   SequenceObj
               

=cut

sub check_intron_match {
  my $self = shift;
  my $cdna = shift;
  
  if ($debug) {
    printf(STDERR "SequenceObj::check_intron_match - comparing introns of %s and %s\n", $self->name, $cdna->name);
  }

  # check if cDNA introns fit with CDS introns
  my $max_cdna_contiguous_introns = 0;
  my $these_introns = 0; # count the introns in this contiguous series
  foreach my $intron ( @{$cdna->sorted_introns}) {
    my $intron_start = $intron->[0];
    # do cDNA and gene share exon start position
    if ( $self->{introns}->{$intron_start} ) { # does CDS intron start exist?
      if ($self->{introns}->{$intron_start} == $cdna->{introns}->{$intron_start} ) { # do intron ends match
	#exact match
	printf(STDERR "SequenceObj::check_intron_match\tExact Intron Match: %s %s\n", $self->name, $cdna->name) if $debug;
	$these_introns++; # count the number of contiguous introns
	if ($these_introns > $max_cdna_contiguous_introns) {$max_cdna_contiguous_introns = $these_introns}
      } else {
	$these_introns = 0;
      }
    }
  }

  # now do it again looking at all the CDS introns in case we have an
  # extra intron in the middle that doesn't match
  #
  # but we don't need to check again if there are no introns matching
  # or if there is only one as there can't be a missed intron in the
  # middle of one intron.
  if ($max_cdna_contiguous_introns <= 1) {return $max_cdna_contiguous_introns}

  my $max_cds_contiguous_introns = 0;
  # check if CDS introns fit with cDNA introns
  $these_introns = 0; # count the introns in this contiguous series
  foreach my $intron ( @{$self->sorted_introns}) {
    my $intron_start = $intron->[0];
    # do cDNA and gene share exon start position
    if ( $cdna->{'introns'}->{"$intron_start"} ) { # does cDNA intron start exist?
      if ($cdna->{'introns'}->{"$intron_start"} == $self->{'introns'}->{"$intron_start"} ) { # do intron ends match
	#exact match
	print STDERR "SequenceObj::check_intron_match\tExact Intron Match\n" if $debug;
	$these_introns++; # count the number of contiguous introns
	if ($these_introns > $max_cds_contiguous_introns) {$max_cds_contiguous_introns = $these_introns}
      } else {
	$these_introns = 0;
      }
    }
  }
  # return min value of the cdna or cds contiguous introns
  if ($max_cds_contiguous_introns > $max_cdna_contiguous_introns) {
    return $max_cds_contiguous_introns;
  } else {
    return $max_cdna_contiguous_introns;
  }
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


=head2 probably_matching_cds

    Function:   get / set the CDSs that probably match a cDNA together with the number of consecutive matching introns
    Returns / Args : ($cds object, number of introns) as a list

=cut

sub probably_matching_cds
  {
    my $self = shift;
    my $cds = shift;
    my $no_of_introns = shift;
    if ($no_of_introns) {
      push ( @{$self->{'probably_matching_cds'}}, [$cds, $no_of_introns] );
    }
    return $self->{'probably_matching_cds'};
  }


=head2 reset_probably_matching_cds

    Function:   reset the array of CDSs that probably match a cDNA
    Args:       none
    Returns:    none

=cut

sub reset_probably_matching_cds
  {
    my $self = shift;
    @{$self->{'probably_matching_cds'}} = ();
  }


=head2 list_of_matched_genes

    Title   :   list_of_matched_genes
    Usage   :   $self->intron_matched_genes
    Function:   returns the list of genes derived from the names of the CDSs in $cdna->probably_matching_cds
    Returns :   list of gene names
    Args    :   
               

=cut

sub list_of_matched_genes {
  my ($self, $cds2gene) = @_;

  my %genes;
  #print "Checking probably_matching_cds for ",$self->name,"\n";
  my @matches = @{$self->probably_matching_cds};
  if (! @matches) {return ()}
  foreach my $match (@matches ) {
    my $cds = $match->[0];

    if (exists $cds2gene->{$cds->name}) {
      $genes{$cds2gene->{$cds->name}} = 1;
    }
  }

  return keys %genes;
}



sub list_of_matched_genes_by_seqname {
  my ($self, $seq_name_regexp) = @_;

  my %genes;
  my @matches = @{$self->probably_matching_cds};

  return () if not @matches;

  foreach my $match (@matches) {
    my $cds = $match->[0];
    my ($gene) = ($cds->name =~ /($seq_name_regexp)/);
  
    if (defined $gene) {
      $genes{$gene} = 1;
    }
  }

  return keys %genes;
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

sub intron_data {
  my ($self, $val) = @_;

  if (defined $val) {
    $self->{introns} = $val;
  }
  return $self->{introns};
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

=head2 sorted_introns

    Title   :  sorted_introns
    Usage   :  my $second_intron = $seq->sorted_introns->[1]
               my $second_intron_start = $second_intron->[0]
    Function:  get sorted introns
    Returns :  sorted array of intron arrays ( sees synopsis )
    Args    :  none
               

=cut

sub sorted_introns
  {
    my $self = shift;
    return $self->{'sorted_introns'};
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

=head2 overlap

    Function:   check to see if the CDS object overlaps with a cDNA object
    Returns :   1 if the overlap, else 0
    Args    :   CDS object               

=cut

# check to see if the CDS object overlaps with a cDNA object in the same strand sense
# return 1 if the overlap, else 0
sub overlap {
  my $self = shift;
  my $cdna = shift;
  
  # check for overlap
  if ($self->strand ne $cdna->strand) {return 0;}
  if( $self->start > $cdna->end ) {
    return 0;
  } elsif( $cdna->start > $self->end ) {
    return 0;
  } else {
    return 1;
  }
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

=head2 coverage

    Function:   get / set coverage score for ESTs
    Returns :   coverage score
    Args    :   coverage score
               

=cut

sub coverage
  {
    my $self = shift;
    my $coverage = shift;

    $self->{'coverage'} = $coverage if $coverage;
    return $self->{'coverage'};
  }

sub check_features {
  my $self = shift;
  my $cdna = shift;

  if (my $SL = $cdna->SL ) {
    if ( $self->SL ) {
      unless( $SL->[0] == $self->SL->[0] ) { #same SL
        print STDERR "SequenceObj::check_features FAIL: ". $self->name . " " . $cdna->name . " : Conficting SLs\n" if $debug;
        return 0;
      }
    } else { 
      if ( $SL->[1] > $self->start ) {
        print STDERR "SequenceObj::check_features FAIL: ". $self->name . " " . $cdna->name . " : Splice leader of cdna within transcript\n" if $debug;
        return 0;
      }
    }
  }


  # if the cDNA extends past the end of the SL-spliced transcript, it could be
  # due to insufficient leader sequence clipping  at the 5' end. The result is that
  # the EST extends a few bps past the TSL site. To address this, we allow
  # the EST to match a transcript that has already been created using TSL data,
  # but do not (later) extend using it
  if ( $self->SL and $cdna->start < $self->start and $self->start - $cdna->start >= 10) {
    print STDERR "SequenceObj::check_features FAIL: ". $self->name . " " . $cdna->name . " : transcript has an SL, and cdna starts before it\n" if $debug;
    return 0;
  }
  
  # and polyA_Site
  my $polyA_site;
  if ( $polyA_site = $cdna->polyA_site ) {
    if ( $self->polyA_site ) {
      unless( $polyA_site->[0] == $self->polyA_site->[0] ) {
        print STDERR "SequenceObj::check_features FAIL: ". $self->name . " " . $cdna->name . " : conflicting polyA sites\n" if $debug;
        return 0;
      }
    } else {
      if ( $polyA_site->[0] < $self->end ) {
        print STDERR "SequenceObj::check_features FAIL: ". $self->name . " " . $cdna->name . " : polyA site of cdna within transcript\n" if $debug;
        return 0;
      }
    }
  }
  
  # . and polyA_signal
  my $polyA_sig;
  if ( $polyA_sig = $cdna->polyA_signal ) {
    if ( $self->polyA_signal ) {
      unless( $cdna->polyA_signal->[0] == $self->polyA_signal->[0] ) {
        print STDERR "SequenceObj::check_features FAIL: ". $self->name . " " . $cdna->name . " : conflicting polyA signals\n" if $debug;
        return 0;
      }
    } else {
      if ( $cdna->polyA_signal->[0] +  2500 > $self->end) {
        # NOT SURE ABOUT THIS ONE!
        print STDERR "SequenceObj::check_features FAIL: ". $self->name . " " . $cdna->name . " : polyA signal of cdna is too far past end of transcript\n" if $debug;    
        return 0;
      }
    }
  }	
  
  # transcript already has polyA and cDNA goes past this
  # if the cDNA extends past the end of the polyA-spliced transcript, it could be
  # due to insufficient polyA sequence clipping  at the 3' end. The result is that
  # the EST extends a few bps past the polyA site. To address this, we allow
  # the EST to match a transcript that has already been created using polyA data,
  # but do not (later) extend using it
  if ( $self->polyA_site and $cdna->end > $self->end and $cdna->end - $self->end >= 10) {
    print STDERR "SequenceObj::check_features: ". $self->name . " " . $cdna->name . " : transcript has a polyA site and cdna goes past the end of it\n" if $debug;
    return 0;
  }
  return 1;
}

1;
