package SequenceObj;

use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'} ;
use Carp;
use strict;

# new 
# expects name , %-ref of exon start/ends and strand ( +/- )

our $debug;

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
					   $transformer->transform_neg_coord( $self->{'feature'}->{"$feature"}->[0])
					 )];
    }
  }

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

sub exon_data
  {
    my $self = shift;
    my $new_exons = shift;
    $self->{'exons'} = $new_exons if $new_exons;
    return $self->{'exons'};
  }

sub first_exon 
  {
    my $self = shift;
    return $self->{'sorted_exons'}->[0];
  }

sub last_exon 
  {
    my $self = shift;
    return $self->{'sorted_exons'}->[-1];
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

sub name
  {
    my $self = shift;
    my $name = shift;
    $self->{'name'} = $name if $name;
    return $self->{'name'};
  }

sub strand
  {
    my $self = shift;
    return $self->{'strand'};
  }

sub sorted_exons
  {
    my $self = shift;
    return $self->{'sorted_exons'};
  }

#this is really a cDNA specific method but cant be bothered to create new class yet !
sub mapped
  {
    my $self = shift;
    my $state = shift;
    if( $state ){
      $self->{'mapped'} = $state;
    }
    return $self->{'mapped'};
  }

sub chromosome
  {
    my $self = shift;
    my $chromosome = shift;
    $self->{'chromosome'} = $chromosome if $chromosome;
    return $self->{'chromosome'};
  }

sub transformer 
  {
    my $self = shift;
    my $transformer = shift;
    $self->{'transformer'} = $transformer if $transformer;

    return $self->{'transformer'};
  }

sub debug
  {
    my $self = shift;
    my $set = shift;
    $debug = 1 if $set;
    return $debug;
  }

sub feature
  {
    my $self = shift;
    my $feature = shift;
    my $data = shift; #@  182772  182773  WBsf01634

    # self=>feature=>SL   =>( x y )
    #              =>polyA_site=>( x y )

    if( $data ) { # adding new feature
      $self->{'feature'}->{ "$feature" } = $data;
    }
    else { 
      return $self->{'feature'}->{$feature};
    }
  }
  
sub SL1
  {
    my ($self, $data) = @_;
    my $type = "SL";
    return $self->feature( $type, $data );
  }

sub SL2
  {
    my ($self, $data) = @_;
    my $type = "SL";
    return $self->feature( $type, $data );
  }


sub SL
  {
    my ($self, $data) = @_;
    my $type = "SL";
    return $self->feature( $type, $data );
  }

sub polyA_site
  {
    my ($self, $data) = @_;
    my $type = "polyA_site";
    return $self->feature( $type, $data );
  }

sub polyA_signal
  {
    my ($self, $data) = @_;
    my $type = "polyA_signal";
    return $self->feature( $type, $data );
  }



1;
