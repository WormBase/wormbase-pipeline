#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-02-13 17:35:32 $ 

use strict;
use lib "/wormsrv2/scripts"  ? "/wormsrv2/scripts" : glob("/wormmsrv1/chaokung/my-scripts");
use Wormbase;
use Ace;
use Getopt::Long;
use Geneace;

my $tace        = &tace;
my $gobj        = init Geneace();
my $geneace_dir = $gobj->test_geneace();
my $db          = Ace->connect(-path => $geneace_dir, -program =>$tace) || print Ace->error; 

push( my @multi_pt, $db->find("Find Multi_pt_data *") );

print scalar @multi_pt;
push(my @other_names, $db->find("Find Gene_name * where Other_name_for AND !(Cb-* OR Cr-*)") );
print scalar @other_names;


__END__
my%Loci_multiPt;

foreach (@multi_pt){
  my $locus;
  if (defined $_ -> Locus_A(1)){
    # push(@{$multiPt_loci{$_}}, $_ -> Locus_A(1));
    push( @{$Loci_multiPt{$_ -> Locus_A(1)}}, $_ );
  }
  if (defined $_ -> Locus_B(1)){
    push(@{$Loci_multiPt{$_ -> Locus_B(1)}}, $_ );
  }
  if (defined $_ -> Locus(1)){
    push(@{$Loci_multiPt{$_ -> Locus(1)}}, $_ );
  }
  if (defined $_ -> B_non_A){
    get_loci($_, "B_non_A");
  }
  if (defined $_ -> A_non_B){
    get_loci($_, "A_non_B");
  }
  if (defined $_ -> Combined){
    get_loci($_, "Combined");
  }
}

sub get_loci {
  my ($obj, $check_point) = @_;
  
  if (defined $obj -> $check_point(2)){ 
    #  push(@{$multiPt_loci{$obj}}, $obj -> $check_point(2));
    push(@{$Loci_multiPt{$obj -> $check_point(2)}}, $obj );
  }
  if (defined $obj -> $check_point(5)){ 
    push(@{$Loci_multiPt{$obj -> $check_point(5)}}, $obj ); 
  }
  if (defined $obj -> $check_point(8)){ 
    push(@{$Loci_multiPt{$obj -> $check_point(8)}}, $obj );
  }
  if (defined $obj -> $check_point(11)){ 
    push(@{$Loci_multiPt{$obj -> $check_point(11)}}, $obj );
  }
  if (defined $obj -> $check_point(14)){ 
    push(@{$Loci_multiPt{$obj -> $check_point(14)}}, $obj );
  }
} 

my @other_names = $gobj->query($db);
print scalar @other_names;


foreach (@other_names){
  if ( exists $Loci_multiPt{$_} ){
    print "Check $_ in Multi_Pt object @{$Loci_multiPt{$_}}\n";
  }
}  
 




