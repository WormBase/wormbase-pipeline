package SequenceObj;

use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'} ;
use Carp;

# new 
# expects name , %-ref of exon start/ends and strand ( +/- )

sub new
  {
    my $class = shift;
    my $name = shift;
    my $exon_data = shift; # \%
    my $strand = shift;

    my $self = {};
    $self->{'name'}   = $name;

    my ($start, $end);
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
    $self->{'start'} = $start;
    $self->{'end'}   = $end;
    $self->{'strand'}= $strand ;

    bless ( $self, $class );
    return $self;
  }


sub check_exon_match 
  {
    my $self = shift;
    my $cdna = shift;

    #use match_code to record match type eg exact exon match to make using this cdna in later steps easy

    #check if cDNA exon fits with gene model
    foreach my $cExonStart (keys %{$cdna->exon_data}) {
      my $gExonS;
      
      # do cDNA and gene share exon start position
      if ( $self->{'exons'}->{"$cExonStart"} ) {
	if ($self->{'exons'}->{"$cExonStart"} == $cdna->{'exons'}->{"$cExonStart"} ) {
	  #exact match
	  print "\tExact Match\n" if $verbose;
	  $match_code = 1;
	}
	#is this final gene exon
	elsif ( $cExonStart == $self->last_exon->[0] ) {
	  if( $cdna->{'exons'}->{"$cExonStart"} > $self->last_exon->[1] ) {
	    print "\tMatch - last SeqObj exon\n" if $verbose;
	    $match_code = 2;
	  }
	  else {
	    print STDERR "MISS : cDNA splices in last SeqObj exon\n" if $verbose;
	  }
	}
	# or final cDNA exon
	elsif ( $cExonStart == $cdna->last_exon->[0] ) {
	  # . . must terminate within gene exon
	  if ( $cdna->{'exons'}->{"$cExonStart"} > $self->{'exons'}->{"$cExonStart"} ) {
	    print STDERR "\tMISS - ",$cdna->name," $cExonStart => ",$cdna->{'exons'}->{$cExonStart}," extends over gene exon boundary\n" if $verbose;
	    return 0;
	  } else {
	    print "\tMatch - last cDNA exon\n" if $verbose;
	    $match_code = 3;
	  }
	}
      }
      # do cDNA and gene share exon end position
      elsif ( ( $gExonS = $self->_exon_that_ends( $cdna->{'exons'}->{"$cExonStart"} ) and ($gExonS != 0) ) ) {
	#	# shared exon end
	
	if ( $gExonS == $self->first_exon->[0] ) { #is this the 1st gene exon 
	  if ( $cExonStart == $cdna->first_exon->[0] ) { # also cDNA start so always match
	    print "\tMatch - 1st exons overlap\n" if $verbose;
	    $match_code = 4;
	  }
	  elsif ( $cExonStart < $self->first_exon->[0] ) { # cDNA exon overlap 1st gene exon
	    print "\tMatch - cDNA exon covers 1st gene exon\n" if $verbose;
	    $match_code = 5;
	    # extends 5'
	  }
	  else {
	    print STDERR "\tMISS - cDNA exon splices in gene exon\n" if $verbose;
	    print STDERR "\t\t",$cdna->name," $cExonStart => ",$cdna{'exons'}->{$cExonStart},"\n" if $verbose;
	    print STDERR "\t\t",$self->name," $gExonS => ",$self->{'exons'}->{$gExonS},"\n" if $verbose;
	    return 0;
	  }
	}
	# exon matched is not 1st of SeqObj
	elsif ( ($cExonStart == $cdna->first_exon->[0] ) and # start of cDNA
		($cExonStart >$gExonS ) ) { # . . . is in SeqObj exon
	  print"\tMatch - 1st exon of cDNA starts in exon of SeqObj\n" if $verbose;
	  $match_code = 6;
	} 
	else {
	  print STDERR "MISS - exon ",$cdna->name," : $cExonStart => ",$cdna{'exons'}->{$cExonStart}," overlaps start of gene exon : $gExonS => ",$self->{'exons'}->{$gExonS},"\n" if $verbose;
	  return 0;
	}
      }# cDNA_wholelyInExon
      elsif ( $self->_cDNA_wholelyInExon($cdna) ) {
	print "Match cDNA contained in exon\n" if $verbose;
	$match_code = 7;
      }
      # cDNA exon overlaps gene 1st exon start and terminate therein
      elsif( ( $cExonStart == $cdna->last_exon->[0] ) and #  last exon of cDNA
	     ( $cExonStart < $self->first_exon->[0] ) and 
	     ( $cdna->last_exon->[1] > $self->first_exon->[0] and $self->first_exon->[0] <$self->first_exon->[1] )
	   ) {
	print "\tcDNA final exon overlaps first exon of gene and end therein\n" if $verbose;
	$match_code = 8;
      }
      # cDNA exon starts in final gene exon and continues past end
      elsif( ($cdna->start > $self->last_exon->[0]) and 
	     ($cdna->start < $self->last_exon->[1]) and 
	     ($cdna->first_exon->[1] > $self->last_exon->[1] )
	   ) {
	print "MATCH : final cDNA exon starts in final gene exon and continues past end\n" if $verbose;
	$match_code = 9;

      }
      else {
	# doesnt match
	print STDERR $cdna->name," doesnt match ",$self->name,"\n";
	return 0;
      }
    }
    $self->add_matching_cDNA( $cdna, $match_code );
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
    return $self->{'start'};
  }

sub end 
  {
    my $self = shift;
    return $self->{'end'};
  }

sub name
  {
    my $self = shift;
    return $self->{'name'};
  }

sub strand
  {
    my $self = shift;
    return $self->{'strand'};
  }


1;
