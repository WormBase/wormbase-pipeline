package Transcript;

use lib -e "/wormsrv2/scripts"  ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'} ;
use Carp;
use Modules::SequenceObj;

@ISA = qw( SequenceObj );

sub new
  {
    my $class = shift;
    my $name = shift;
    my $CDS = shift;
    my $transformer = $CDS->transformer;

    my $self = SequenceObj->new($name, $CDS->exon_data, $CDS->strand);

    bless ( $self, $class );

    $self->transformer($transformer) if $transformer;
    return $self;
  }

sub map_cDNA
  {
    my $self = shift;
    my $cdna = shift;

    # check for overlap
    if( $self->start > $cdna->end ) {
      return 0;
    }
    elsif( $cdna->start > $self->end ) {
      return 0;
    }
    else {
      #this must overlap - check Splice Leader
      my $SL;
      if( $SL = $cdna->SL ){
	if( $self->SL ) {
	  unless( $SL->[0] == $self->SL->[0] ) {	    #same SL
	    print STDERR "Conficting SLs ",$self->name, "\t",$cdna->name,"\n" if $debug;
	    return 0;
	  }
	}
	else { 
	  if( $SL->[1] > $self->start ) {
	    print STDERR $cdna->name, " Splice Leader within ", $self->name, "\n" if $debug;
	    return 0;
	  }
	}
      }

      # reject cdna's that start before transcript SL
      if( $self->SL and $cdna->start < $self->start ) {
	return 0;
      }

      # and polyA_Site
      my $polyA_site;
      if( $polyA_site = $cdna->polyA_site ) {
	if( $self->polyA_site ) {
	  unless( $polyA_site->[0] == $self->polyA_site->[0] ) {
	    print STDERR "Conficting polyA_sites ",$self->name, "\t",$cdna->name,"\n" if $debug;
	    return 0;
	  }
	}
	else {
	  if( $polyA_site->[0] < $self-end ) {
	    print STDERR $cdna->name, " polyA_site within ", $self->name, "\n" if $debug;
	    return 0;
	  }
	}
      }

      # . and polyA_signal
      my $polyA_sig;
      if( $polyA_sig = $cdna->polyA_signal ) {
	if( $self->polyA_sig ) {
	  unless( $polyA_sig->[0] == $self->polyA_sig->[0] ) {
	    print STDERR "Conficting polyA_signals ",$self->name, "\t",$cdna->name,"\n" if $debug;
	    return 0;
	  }
	}
	else {
	  if( $polyA_sig->[0] + 30 > $self-end) {
	    print STDERR $cdna->name, " polyA_signal within ", $self->name, "\n" if $debug;
	    return 0;
	  }
	}
      }	

      # transcript already has polyA and cDNA goes past this
      if( $self->polyA_site and $cdna->end > $self->end ){
	return 0;
      }
      
      #check exon matching
      my $match = $self->check_exon_match( $cdna );
      if( $match == 1 ) {
	$match = $cdna->check_exon_match( $self );
	unless ($match == 0) { 
	  $self->add_matching_cDNA( $cdna );
	  $match = 1;
	}
      }
      $self->SL( $SL ) if $SL;
      $self->polyA_site( $polyA_site )  if $polyA_site;
      $self->polyA_signal( $polyA_sig ) if $polyA_sig;
      return $match;
    }
  }



sub add_matching_cDNA
  {
    my $self = shift;
    my $cdna = shift;
    #print STDERR $cdna->name," matches ",$self->name,"\n";
    push( @{$self->{'matching_cdna'}},$cdna);

###### modify current transcript structure to incorporate new cdna #####

#   match_codes
#  1 = Exact Match
# *2 = last SeqObj exon                   ( so may extend past exon end )
#  3 = last cDNA exon                     ( so may stop within exon )
# *4 = 1st exons end in same place
# *5 = cDNA exon covers 1st gene exon
#  6 = 1st exon of cDNA starts in exon of SeqObj
#  7 = Match cDNA contained in exon
# *8 = cDNA final exon overlaps first exon of gene and end therein
# *9 = final cDNA exon starts in final gene exon and continues past end
# *10 = 5'UTR exon
# *11 = 3'UTR exon

########################################################################

    
    foreach $exon ( @{$cdna->sorted_exons} ) {
      my $match_code = $exon->[2];
      next if  ($match_code == 1 or $match_code == 3 or $match_code == 6 or $match_code == 7);

      #extend 3'
      if( $match_code == 2 or $match_code == 9 ) {
	# cdna overlaps last exon so extend.
	my $last_exon_start = $self->last_exon->[0];
	$self->exon_data->{"$last_exon_start"} = $exon->[1];
      }
      #extend 5'
      elsif( $match_code == 4  or $match_code == 5 or $match_code == 8 ) {
	# 1st exons overlap so extend 5'
	my $curr_start = $self->start;

	next if( $curr_start < $exon->[0] );

	my $exon_end = $self->sorted_exons->[0]->[1];

	delete $self->exon_data->{$curr_start};
	$self->exon_data->{"$exon->[0]"} = $exon_end;
      }
      elsif( $match_code == 10  or $match_code == 11 ) {
	#add exon to UTR
	$self->exon_data->{"$exon->[0]"} = $exon->[1];
      }
    }
    # reset start end etc . . 
    $self->sort_exons;
  }

sub report
  {
    my $self = shift;
    my $fh = shift;
    my $coords = shift;
    
    my $real_start = $self->start;
    my $real_end = $self->end;

    if( $self->strand eq "-" ) {
      $real_start = $self->transformer->revert_to_neg( $self->start );
      $real_end = $self->transformer->revert_to_neg( $self->end );
    }

    my @clone_coords = $coords->LocateSpan($self->chromosome, $real_start,$real_end );

    # output S_Parent for transcript
    print $fh "\nSequence : $clone_coords[0]\n";
    print $fh "Transcript \"",$self->name,"\" $clone_coords[1] $clone_coords[2]\n";

    # . . and the transcript object
    print $fh "\nTranscript : \"", $self->name,"\"\n";
    our $start = $self->start - 1;
    our $end = $self->end + 1;
    my $exon_1 = 1;
    foreach $exon ( @{$self->sorted_exons} ) {    
      $exon->[0] = $exon->[0] - $start;
      $exon->[1] = $exon->[1] - $start;

      print $fh "Source_exons\t$exon->[0]\t$exon->[1]\n";
    }

    # . . and Matching_cDNA
    foreach (@{$self->{'matching_cdna'}}) {
      print $fh "Matching_cDNA \"",$_->name,"\"\n";
    }
    # .. and method
    print $fh "Method Coding_transcript\n";
  }









1;
