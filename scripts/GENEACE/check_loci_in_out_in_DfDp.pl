#!/usr/local/bin/perl5.8.0 -w
# 
# check_loci_in_out_in_DfDp.pl
#
# by Chao-Kung Chen
#
# Script to check for the loci no longer as inside/outside locus of a Df/Dp based on modified gmap via interpolation
#
# Last updated by: $Author: mt3 $
# Last updated on: $Date: 2005-12-09 13:25:34 $


use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Ace;

######################
# global variables   # 
######################

my $tace = &tace;  
my $geneace_dir = "/nfs/disk100/wormpub/DATABASES/geneace";

my $locus_2_gmap=<<EOF;
  Table-maker -p "/nfs/disk100/wormpub/DATABASES/geneace/wquery/locus_2_gmap.def" quit
EOF

open (FH, "echo '$locus_2_gmap' | $tace $geneace_dir | ") || die "Couldn't access geneace\n";

my %locus_chrom_pos;
while (<FH>){
  chomp $_;
  if ( $_ =~ /^\"/ ) {	
    my ($locus, $chrom, $pos) = split(/\s+/, $_) if $_ =~ /^\"/;
    $locus =~ s/\"//g;
    $chrom =~ s/\"//g;
    $pos =~ s/\"//g;		
    push (@{$locus_chrom_pos{$locus}}, $chrom, $pos);
  }
}
close FH;

my $rearr_ends=<<EOF;
  Table-maker -p "/nfs/disk100/wormpub/DATABASES/geneace/wquery/rearr_L_R_ends.def" quit
EOF

open (FH, "echo '$rearr_ends' | $tace $geneace_dir | ") || die "Couldn't access geneace\n";

my %Rearr_map_ends;
while (<FH>){
  chomp $_;
  my ($rearr, $map, $L_end, $R_end) = split(/\s+/, $_) if $_ =~ /^\"/;
  if ($R_end){
    $rearr =~ s/\"//g;
    $map   =~ s/\"//g;
    $L_end =~ s/\"//g;
    $R_end =~ s/\"//g;
    push(@{$Rearr_map_ends{$rearr}}, $map, $L_end, $R_end);
   # print "$rearr, $map, $L_end, $R_end\n";
  }
}
close FH;

my $db = Ace->connect(-path  => $geneace_dir,
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();};

my $rearr_query = "Find Rearrangement *df* where Map";
push(my @rearr, $db->find( $rearr_query ) );

my %seen_rearr=();
my $in_count = 0; 
my $out_count = 0;
my $has_map = 0;

foreach (@rearr){
  if ( defined $_ ->Locus_inside(1) ){
    foreach my $e ($_ ->Locus_inside(1) ){
      $in_count++;
      $has_map++ if $locus_chrom_pos{$e}->[0];
      if ( exists $locus_chrom_pos{$e} && $locus_chrom_pos{$e}->[0] eq $Rearr_map_ends{$_}->[0] && 
         !($locus_chrom_pos{$e}->[1] > $Rearr_map_ends{$_}->[1] && $locus_chrom_pos{$e}->[1] < $Rearr_map_ends{$_}->[2]) ){
	print "$e (new pos $locus_chrom_pos{$e}->[1]) is NOT inside of $_ (pos range @{$Rearr_map_ends{$_}}[1..2]) \n";
	push(@{$seen_rearr{$_}},$_);
      }
    }
  }
  if ( defined $_ ->Locus_outside(1) ){
    foreach my $e ( $_ ->Locus_outside(1) ){
      $out_count++;
      $has_map++ if $locus_chrom_pos{$e}->[0];
      if ( exists $locus_chrom_pos{$e} && $locus_chrom_pos{$e}->[0] eq $Rearr_map_ends{$_}->[0] && 
	   ($locus_chrom_pos{$e}->[1] > $Rearr_map_ends{$_}->[1] && $locus_chrom_pos{$e}->[1] < $Rearr_map_ends{$_}->[2]) ){
	print "$e (new pos $locus_chrom_pos{$e}->[1]) is NOT outside of $_ (pos range @{$Rearr_map_ends{$_}}[1..2]) \n";
	push(@{$seen_rearr{$_}}, $_);
      }
    }
  }
}   
 
print "\n\nThe following ", scalar keys %seen_rearr, " Df objects have loci map positions not matching the Locus_inside/Locus_outside information\n\n";

my $count = 0;
foreach (keys %seen_rearr){
  print "$_ has ", scalar @{$seen_rearr{$_}}, " unmatched\n";
  $count += scalar @{$seen_rearr{$_}};
}
print "\n-------------------------------\n";
print "Total: $count unmatched out of ", $in_count + $out_count, " loci ($has_map have map position)  in ", scalar @rearr, " Df objects\n";
