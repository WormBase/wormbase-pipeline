#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-05-12 16:07:41 $ 

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;
use Getopt::Long;
use GENEACE::Geneace;

my $tace        = &tace;
my $ga          = init Geneace();
my $geneace_dir = $ga->geneace();
my $db          = Ace->connect(-path => $geneace_dir, -program =>$tace) || print Ace->error; 
my $rundate = rundate();

open(LOG, ">/wormsrv2/logs/dubious_multiPt_2_locus.$rundate") || die $!;
open(ACE, ">/wormsrv1/geneace/CHECKS/multiPt_2_locus.ace");

push( my @multi_pt, $db->find("Find Multi_pt_data *") );

my %Gene_info = $ga->gene_info();

my %Gene_id_2_multiPt_A = $ga->gene_id_has_multi_pt($db); # mutipt objs listed in each locus

my %Gene_id_2_multiPt_B; # key: locus, values: all loci found in a multipt obj
my %multi_geno;

my $counter = 0;

foreach (@multi_pt){
  my (@checks, $e);
  $multi_geno{$_} = $_ -> Genotype(1) if defined $_ -> Genotype(1);
  $multi_geno{$_} = "NA" if !defined $_ -> Genotype(1);

  @checks = qw(Gene_A Gene_B Gene);
  foreach $e (@checks){
    if (defined $_ -> $e(1)){
      push( @{$Gene_id_2_multiPt_B{$_ -> $e(1)}}, $_ ); # key is locus
    }
  }

  @checks = qw(B_non_A A_non_B Combined);
  foreach $e (@checks){
    if (defined $_ -> $e){
      get_loci($_, $e);
    }
  }
}

foreach (keys %Gene_id_2_multiPt_A){ # key : locus obj, value: multi-pt of that locus obj
  my @unique_B = ();
  my @unique_A = $ga->get_unique_from_array( @{$Gene_id_2_multiPt_A{$_}} ); # unique multi in a locus obj
  @unique_B = $ga->get_unique_from_array( @{$Gene_id_2_multiPt_B{$_}} ) if exists $Gene_id_2_multiPt_B{$_}; # unique multi in a locus obtained from a multi-pt obj

  if ( !exists $Gene_id_2_multiPt_B{$_} ){
    foreach my $ea ( @unique_A ) {
      $multi_geno{$ea} = "NA" if !exists $multi_geno{$ea};
      $counter++;
      print LOG "$counter. ERROR: $ea ($multi_geno{$ea}) should not be linked to $_?\n";
      print ACE "\n//ERROR: $ea ($multi_geno{$ea}) should not be linked to $_ ($Gene_info{$_}{'CGC_name'}) ? ===(1)\n";
      print ACE "Gene : \"$_\"\n";
      print ACE "-D Multi_point \"$ea\"\n";
    }
  }

  my (%unique_A, %unique_B, @diff);
  foreach (@unique_A){$unique_A{$_}++};  # turn array element to hash key for quick look up

  if (@unique_B){

    foreach (@unique_B){$unique_B{$_}++};
    @diff = $ga->array_comp(\@unique_A, \@unique_B, "diff");

    foreach my $e (@diff){
      $multi_geno{$e} = "NA" if !exists $multi_geno{$e};	
      if (!exists $unique_A{$e} ){
	$counter++;
        print LOG "$counter. CHECK: $e ($multi_geno{$e}) should be linked to $_ ($Gene_info{$_}{'CGC_name'}) ?\n";
      }
      if (!exists $unique_B{$e} ){
	$counter++;
        print LOG "$counter. ERROR: $e ($multi_geno{$e}) should not be linked to $_ ($Gene_info{$_}{'CGC_name'}) ?\n";
	
	print ACE "\n//ERROR: $e ($multi_geno{$e}) should not be linked to $_ ($Gene_info{$_}{'CGC_name'}) ? ===(2)\n";
	print ACE "Gene : \"$_\"\n";
	print ACE "-D Multi_point \"$e\"\n";
      }
    }	
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
       push(@{$Gene_id_2_multiPt_B{$obj -> $check_point($i+2)}}, $obj );
    }
  }
}



__END__

=head2 NAME - check_dubious_multi_pt_2_locus.pl

=head3 <DESCRIPTION> 

B<>%Gene_ids_2_multiPt_A:  # mutipt objs listed in each locus
   (key: locus, values: multipt obj)

B<>%Gene_id_2_multiPt_B :  # loci listed in each multipt obj 
   (key: locus, values: multipt obj)

When values of these two hashes are compared, the differences tell you
what is wrong.

