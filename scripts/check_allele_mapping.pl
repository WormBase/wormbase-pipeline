#!/usr/local/bin/perl5.6.1 -w
# 
# check_allele_mapping.pl
#
# by Chao-Kung Chen
#
# Script to check if allele sequence from allele mapping script is the same as current 
#
# Last updated by: $Author: ck1 $
# Last updated on: $Date: 2003-07-09 11:36:52 $

use strict;
use lib "/wormsrv2/scripts/"; 
use Wormbase;


my $tace = &tace;
my $curr_db = "/nfs/disk100/wormpub/DATABASES/current_DB";
my $ga_dir = "/wormsrv1/geneace/";
my (@mapping, $allele, $seq, %allele_seq, %allele_seq_map);

my $current = `grep "NAME WS" $curr_db/wspec/database.wrm`; $current =~ s/NAME WS//; chomp $current;
$current += 1;
print $current, "\n";


@mapping = `echo "table-maker -p /wormsrv1/geneace/wquery/allele_has_seq.def" | $tace $ga_dir`;

foreach (@mapping){
  chomp;
  if ($_ =~ /^\"(.+)\"\s+\"(.+)\"/){ # $1 = allele, $2 = seq
   $allele_seq{$1} = $2;
  }
}
foreach (keys %allele_seq){
 # print "$_ -> $allele_seq{$_}\n";
}

my $update = "/wormsrv2/autoace/MAPPINGS/map_alleles_geneace_update".$current.".ace";

open(IN, "/wormsrv2/autoace/MAPPINGS/map_alleles_geneace_update".$current.".ace") || die $!;

while(<IN>){
  if ($_ =~ /^Allele\s+:\s+\"(.+)\"/){
    $allele = $1;
  }
  if ($_ =~ /^Sequence\s+\"(.+)\"/){
    $seq = $1;
  }  
  $allele_seq_map{$allele} = $seq;
}

foreach (keys %allele_seq_map){
  print "$_ from mapping has a different seq than the one in geneace\n" 
  if (exists $allele_seq{$_} && $allele_seq_map{$_} ne $allele_seq{$_})

}
    
  

