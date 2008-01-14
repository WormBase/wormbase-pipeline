##############
# Features.pm
# to contain the filters for use before alignment to the genome
# by Michael Han
# 
# Last updated by: $Author: pad $
# Last updated on: $Date: 2008-01-14 16:40:43 $

package Features;

use strict;
use warnings;

#############
# struct to hold the
# Feature data

use Class::Struct;

struct Feature =>{
    id         => '$',
    sequence   => '$',
    start      => '$',
    stop       => '$',
    psize      => '$',
    length     => '$',
    type       => '$',
    description=> '$',
    type2      => '$'
};


#############
# detect polyA the classic way
# trim_polyA(sequence)

sub trim_polyA{
  my($seq,$id)=@_;
  my $flag = "0";
  my $psize = length $seq;
  my ($maskedb,$tail)=0;
  $seq=~tr/acgt/ACGT/; #upcase sequence
  
  #Hidden polyAs and bad clipping.   AAAAAAAXYZPQR
  my ($probseq,$poorqual,$polyAlength,$polyAstart,$polyAend,$rubstart,$badseq,$rubbish_ratio)="";
  if ($seq=~/(AAAAAAAAAA+([^A]+\w+))/) {
    $probseq = length $1;
    $poorqual = length $2;
    $badseq = $2;
    
    my $poorqualscore = &check_poor_qual($badseq);#uses a small subroutine to assign a score to the post polyA sequence.
    if (($poorqualscore < 33) && ($poorqual > 50)) {
      print "//$id contains internal polyA run longer than 10bp but does not look like a real polyA!!\n\n";
    }
    elsif (($poorqualscore > 32) or ($poorqual < 25)) { #mask poor quality under 25bp scoring fails short seqs.
      $polyAlength = ($probseq - $poorqual);
      $polyAstart = ($psize - $probseq +1);
      $polyAend = ($polyAstart + $polyAlength -1);
      $rubstart = ($polyAend + 1);
      $flag = "1";
      
      if ($flag > 0) {
	my @features;
	my $feature1=Feature->new(
				  id          => $id,
				  sequence    => $seq,
				  start       => $polyAstart,
				  stop        => $polyAend,
				  psize       => $psize,
				  length      => $polyAlength,
				  type        => 'polyA',
				  description => 'polyA tail sequence',
				  type2       => 'polyA'
				 );
	my $feature2=Feature->new(
				  id          => $id,
				  sequence    => $seq,
				  start       => $rubstart,
				  stop        => $psize,
				  psize       => $psize,
				  length      => ($psize - $rubstart+1),
				  type        => 'low',
				  description => '['.$badseq.']'.'- low-complexity or poor quality/unclipped sequence(Score::'.$poorqualscore.'). Removed from BLAT pipeline.',
				  type2       => 'low'
				 );
	my $feature3=Feature->new(
				  id          => $id,
				  sequence    => $seq,
				  start       => ($polyAstart),
				  stop        => ($polyAstart+1),
				  psize       => $psize,
				  length      => 1,
				  type        => 'polyA_site',
				  description => 'polyA site',
				  type2       => 'polyA_site'
				 );
	return ($feature1,$feature2,$feature3);
      }
    }
      elsif (($poorqualscore > 32) && ($poorqual > 150)) {
	print "Warning, this sequence appears to be low quality/poly nucleotide runs and need investigating.\n"
      }
  }
  
  if (($seq=~/(A+)$/) && ($flag < 1)) {
    $tail=$1;
    $maskedb=length $1;
    }
  # generate 2 features polyA site and tail if size > 5
  if ($maskedb > 5) {
    
    # masking
    my $masked_tail= "n"x$maskedb;
    my @features;
    $seq=~s/$tail$/$masked_tail/i;
    
    my $feature1=Feature->new(
			      id          => $id,
			      sequence    => $seq,
			      start       => ($psize - $maskedb+1),
			      stop        => $psize,
			      psize       => $psize,
			      length      => $maskedb,
			      type        => 'polyA',
			      description => 'polyA tail sequence',
			      type2       => 'polyA'
			     );
    my $feature2=Feature->new(
			      id          => $id,
			      sequence    => $seq,
			      start       => ($psize - $maskedb),
			      stop        => ($psize - $maskedb+1),
			      psize       => $psize,
			      length      => 1,
			      type        => 'polyA_site',
			      description => 'polyA site',
			      type2       => 'polyA_site'
			     );
    ################
    # in case we migth look for polyA motifs in future
    #
    #
    #	if ($seq=~/(AATAAA)\w{5,30}?n+$/i){
    #	    my $offset=length $`;
    #	    my $motif=$1;
    #	    my $prettyseq=$seq;
    #	    $prettyseq=~s/(AATAAA)(\w{5,30}?n+$)/ \[$1\] $2/;
    #	   print '='x20,"\nfound polyA motif $motif at $offset\n".'='x20,"\nin\n$prettyseq\n";
    #	}
    return ($feature1,$feature2);
    }
  else {return undef}
}



##############
# idea borrowed from the CGI script
# at /nfs/WWWdev/INTWEB_docs/cgi-bin/Projects/C_elegans/tsl_search

sub get_tls{
  my %SL    = (
	       'SL1',  'GGTTTAATTACCCAAGTTTGAG',
	       'SL2',  'GGTTTTAACCCAGTTACTCAAG',
	       'SL2a', 'GGTTTATACCCAGTTAACCAAG',
	       'SL2b', 'GGTTTTAACCCAGTTTAACCAAG',
	       'SL2c', 'GGTTTTAACCCAGTTACCAAG',
	       'SL2d', 'GGTTTTTACCCAGTTAACCAAG',
	       'SL2e', 'GGTTTAAAACCCAGTTAACAAG',
	       'SL2f', 'GGTTTTAACCCAGTTAACCAAG',
	       'SL2g', 'GGTTTTAACCAGTTAACTAAG',
	       'SL2h', 'GGTTTTAACCCATATAACCAAG',
	       'SL2i', 'GGTTTTAACCCAAGTTAACCAAG',
	       'SL2j', 'GGTTTAAAACCCAGTTACCAAG',
	       'SL2k', 'GGTTTTAACCCAGTTAATTGAG',
          );
  my ($seq,$id)=@_;
  $seq=~tr/acgt/ACGT/; #upcase sequence
  my $psize=length $seq;
  
  # searches for n-mers of the TLS (max 23)
  for (my $i=23;$i>6;$i--){
    
    # loops through the TLS sequences 
    foreach my $slname (keys %SL){
      
      next if $i > length $SL{$slname};
      my $nmer=substr($SL{$slname},-$i);
      
      if ($seq=~ /^($nmer)/i) {
	my $stop=length $i;
	my $mask= 'n'x(length $nmer);
	
	$seq=~s/^$nmer/$mask/i; # trim sequence
	$slname=~/(SL\d)/; # extract type
	my $sltype=$1;
	my $feature=Feature->new(
				 id          => $id,
				 sequence    => $seq,
				 start       => 1,
				 stop        => (length $nmer),
				 psize       => $psize,
				 length      => (length $nmer),
				 type        => $sltype,
				 description => "$slname trans-spliced leader sequence",
				 type2       => 'TSL'
				);
	
	return $feature;
      }
	    }
  }
  return undef;
}

###################
# generate ACE from Feature

sub ace_feature{
  my ($feature)=@_;
  
  my $template="Sequence : \"%s\"\n".
    "Feature_data \"%s:%s\" 1 %d\n".
      "\n".
	"Feature_data : \"%s:%s\"\n".
	  "Feature %s  %d %d %d \"%s\"\n\n";
  
  my $id=$feature->id;
  my $type2=$feature->type2;
  
  return sprintf($template,$id,
		 $id,$type2,$feature->psize,
		 $id,$type2,
		 $feature->type,$feature->start,$feature->stop,$feature->length,$feature->description);
}


#####################################################################
# Sequence quality estimator.
# Subroutine provided by ar2 to give a score for a piece of sequence.
# Below 33 and the sequence is probably genomic.

sub check_poor_qual {
  my $dna = shift;
  my $dna_length = length ($dna);
  my $n = $dna =~ tr/N/N/ ;
  my $q = sprintf("%2d",$n/$dna_length * 100);
  
  my $n100 = substr($dna,0,100) =~ tr/N/N/ ;
  my $q100 = sprintf("%2d",$n100);
  
  my $rubbish =0;
  foreach ( qw( A C T G N ) ) {
    my @n = $dna =~ /$_+/g;
    my @sn = sort { length($b) <=> length($a) } @n;
    foreach (@sn) {
      if ( length($_) > 5 ) {
	$rubbish +=  length($_);
      }
    }
  }
  my $rubbish_ratio;
  if ( $rubbish == 0 or $dna_length == 0) {
    $rubbish_ratio = 0;
  } 
  else {
    $rubbish_ratio = sprintf("%2d",$rubbish/$dna_length *100);
  }
  return $rubbish_ratio;
}

#################
# create ACE from polyA and TSL

sub annot{
  my ($seq,$id)=@_;
  my @features=Features::trim_polyA($seq,$id);
  push(@features,Features::get_tls($seq,$id));
  my $ace='';
  foreach my $feature (@features){
    $ace=$ace.Features::ace_feature($feature) if $feature;
  }
  return $ace;
}


1;
__END__

=head1 NAME

Features.pm - Module for Feature_data processing

=head1 SYNOPSIS

use Features;

Features::annot($sequence,$id)


=head1 DESCRIPTION

Module containing functions to filter and detect polyA and trans-spliced leader (TSL) sequences.
It contains an ugly helper funtion to create ACE format which is used later
for the feature pipeline (which is also slightly misplaced in the module).

the Features_data is stored in a Feature structure which also contains the masked sequence.

=head3 struct Feature

    id         => '$', parent id ($id)
    sequence   => '$', masked sequence
    start      => '$',
    stop       => '$',
    psize      => '$', parent size
    length     => '$', feature length
    type       => '$', ACE method name
    description=> '$',
    type2      => '$'  ACE class name


=head3 trim_polyA($sequence,$id)

returns: upto 3 Features, one for polyA tail, one for polyA site if found else 
undef and a low-complexity feature if the read is poorly clipped.

minimum length is 6 bases

=head3 get_tls($sequence,$id)

returns: a Feature for TLS or undef if nothing found

searches for suffixes of known TSL sequences at the beginning of the sequence.
The minimum suffix length is 6.
In addition it will return only the longest match


=head3 ace_feature($feature)

returns: an ACE string

Generic method to handle Feature to ACE conversion

Note $feature is supposed to be an Feature struct

=head3 annot($sequence,$id)

returns an ace string of polyA sites/tails and TSL

=back

=head1 AUTHOR

$Author: pad $

=head1 VERSION

$Date: 2008-01-14 16:40:43 $

=cut
