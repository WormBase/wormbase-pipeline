#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-03-11 16:50:38 $ 

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Getopt::Long;
use lib "/nfs/team71/worm/ck1/WORMBASE_CVS/scripts";
use Geneace;

my $tace        = &tace;
my $ga          = init Geneace();
my $geneace_dir = $ga->geneace();
my $db          = Ace->connect(-path => $geneace_dir, -program =>$tace) || print Ace->error; 
my $rundate = rundate();

open(LOG, ">/wormsrv2/logs/dubious_multiPt_2_locus.$rundate") || die $!;
open(ACE, ">/wormsrv1/geneace/CHECKS/multiPt_2_locus.ace");

push( my @multi_pt, $db->find("Find Multi_pt_data *") );

my %Locus_2_multiPt_A = $ga->loci_have_multi_pt($db); # mutipt objs listed in each locus

my @other_names = $ga->other_name($db, "other");

my %Locus_2_multiPt_B; # key: locus, values: all loci found in a multipt obj
my %multi_geno;

foreach (@multi_pt){
  my (@checks, $e);
  $multi_geno{$_} = $_ -> Genotype(1) if defined $_ -> Genotype(1);	
  @checks = qw(Locus_A Locus_B Locus);
  foreach $e (@checks){
    if (defined $_ -> $e(1)){
      push( @{$Locus_2_multiPt_B{$_ -> $e(1)}}, $_ );
    }
  }
  
  @checks = qw(B_non_A A_non_B Combined);
  foreach $e (@checks){
    if (defined $_ -> $e){
      get_loci($_, $e);
    }
  }
}

foreach (keys %Locus_2_multiPt_A){ # multi-pt of a locus
  my @unique_B = ();
  my @unique_A = $ga->get_unique_from_array( @{$Locus_2_multiPt_A{$_}} ); # unique multi in a locus
  @unique_B = $ga->get_unique_from_array( @{$Locus_2_multiPt_B{$_}} ) if exists $Locus_2_multiPt_B{$_}; # unique multi in a locus

  if ( !exists $Locus_2_multiPt_B{$_} ){
    foreach my $ea ( @unique_A ) {
      print LOG "ERROR: $ea ($multi_geno{$ea}) should not be linked to $_?\n";
      print ACE "\n//ERROR: $ea ($multi_geno{$ea}) should not be linked to $_?\n";
      print ACE "Locus : \"$_\"\n";
      print ACE "-D Multi_point \"$ea\"\n";
    }
  }

  my (%unique_A, %unique_B, @diff);
  foreach (@unique_A){$unique_A{$_}++};  # turn array element to hash key for quick look up

  if (@unique_B){
    foreach (@unique_B){$unique_B{$_}++};
    @diff = $ga->array_comp(\@unique_A, \@unique_B, "diff");
    #print "$_ -> @diff ###DIFFS\n" if @diff;
    foreach my $e (@diff){
      $multi_geno{$e} = "NA" if !exists $multi_geno{$e};	
      if (!exists $unique_A{$e} ){
        print LOG "CHECK: $e ($multi_geno{$e}) should be linked to $_ ?\n";
      }
      if (!exists $unique_B{$e} ){
        print LOG "ERROR: $e ($multi_geno{$e}) should not be linked to $_?\n";
	print ACE "\n//ERROR: $e ($multi_geno{$e}) should not be linked to $_?\n";
	print ACE "Locus : \"$_\"\n";
	print ACE "-D Multi_point \"$e\"\n";
      }
    }
  }	
}    

print LOG "\n\nChecking locus names (other name) appear in multi-pt object\n";
print LOG "-----------------------------------------------------------\n\n";

foreach (@other_names){
  if ( exists $Locus_2_multiPt_B{$_} ){
    print LOG "CHECK: $_ found in Multi_Pt object @{$Locus_2_multiPt_B{$_}}. $_ an Other_name\n";
  }
}  
 
###############################
#    s u b r o u t i n e s
###############################

sub get_loci {
  my ($obj, $check_point) = @_;
  
  # currently looking for up to 10 tiers of multi point mapping, should be really enough
  for(my $i=0; $i<30; $i=$i+3){
    if (defined $obj -> $check_point($i+2)){   
       push(@{$Locus_2_multiPt_B{$obj -> $check_point($i+2)}}, $obj );
    }
  }
}



__END__

=head2 NAME - geneace_check.pl  

=head3 <DESCRIPTION> 

B<>%Locus_2_multiPt:  # mutipt objs listed in each locus
   (key: locus, values: multipt obj)

B<>%Locus_2_multiPt_B :  # loci listed in each multipt obj 
   (key: locus, values: multipt obj)

When values of these two hashes are compared, the differences tell you
what is wrong.

