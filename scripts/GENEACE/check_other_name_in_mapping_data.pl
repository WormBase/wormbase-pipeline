#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-03-19 15:27:03 $ 

use strict;
use lib "/wormsrv2/scripts"  ? "/wormsrv2/scripts" : glob("/wormmsrv1/chaokung/my-scripts");
use Wormbase;
use Ace;
use Getopt::Long;
use lib "/nfs/team71/worm/ck1/WORMBASE_CVS/scripts";
use GENEACE::GGeneace;

my $tace        = &tace;
my $ga          = init Geneace();
my $geneace_dir = $ga->geneace();
my $db          = Ace->connect(-path => $geneace_dir, -program =>$tace) || print Ace->error; 

push( my @multi_pt, $db->find("Find Multi_pt_data *") );

my %loci_2_multi = $ga->loci_have_multi_pt($db); # mutipt objs listed in each locus
my @other_names = $ga->other_name($db, "array");

my %Loci_multiPt; # key: locus, values: all loci found in a multipt obj


foreach (@multi_pt){
  my (@checks, $e);

  @checks = qw(Locus_A Locus_B Locus);
  foreach $e (@checks){
    if (defined $_ -> $e(1)){
      push( @{$Loci_multiPt{$_ -> $e(1)}}, $_ );
    }
  }
  
  @checks = qw(B_non_A A_non_B Combined);
  foreach $e (@checks){
    if (defined $_ -> $e){
      get_loci($_, $e);
    }
  }
}

foreach (keys %Loci_multiPt){ # loci grepped from each multi-pt 
#  print "$_ => @@@@\n"; 
  my @unique = $ga->get_unique_from_array( @{$Loci_multiPt{$_}} );

  my %unique;
  foreach (@eunique){$unique{$_}++};  # turn array element to hash key for quick look up

  my @diff = $ga->array_comp(\@unique, \@{$loci_2_multi{$_}}, "diff");

#  print "$_ -> @diff ###DIFFS\n" if @diff;
  foreach my $e (@diff){
    if (!exists $unique{$e}){
      print "$e is linked to $_, OK?\n";
    }
    else{
      print "$e should be linked to $_?\n";
    }
  }
 
 # print "   -> @{$loci_2_multi{$_}} ~~~~~~~~~~~~MULTI of each LOCUS\n";
}

#foreach (@other_names){
#  if ( exists $Loci_multiPt{$_} ){
#    print "Check $_ in Multi_Pt object @{$Loci_multiPt{$_}}\n";
#  }
#}  
 


###############################
#    s u b r o u t i n e s
###############################

sub get_loci {
  my ($obj, $check_point) = @_;
  
  # currently looking for up to 10 tiers of multi point mapping, should be really enough
  for(my $i=0; $i<10; $i=$i+3){
    if (defined $obj -> $check_point($i+2)){   
       push(@{$Loci_multiPt{$obj -> $check_point($i+2)}}, $obj );
    }
  }
}



__END__

=head2 NAME - geneace_check.pl  

=head3 <DESCRIPTION> 

B<>%loci_2_multi :  # mutipt objs listed in each locus
   (key: locus, values: multipt obj)

B<>%Loci_multiPt :  # loci listed in each multipt obj 
   (key: locus, values: multipt obj)

When values of these two hashes are compared, the differences tell you
what is wrong.

