#!/usr/local/bin/perl5.6.1 -w
# 
# parse_predicted_tRNA.pl
#
# by Chao-Kung Chen
#
# Script to parse columns of predicted tRNA by tRNASCAN-SE
#
# Last updated by: $Author: ck1 $
# Last updated on: $Date: 2003-09-02 15:33:33 $

use strict;
use lib "/wormsrv2/scripts/"; 
use Wormbase;
use Getopt::Long;

###########
# variables
###########

my $tace = &tace;  
my $curr_db = "/nfs/disk100/wormpub/DATABASES/current_DB";
my ($help, $file);

GetOptions ("h|help"        => \$help,
            "f|file=s"     => \$file,
	   );

if ($help){
  system("perldoc $0");
  exit(0);
}

#############################################
# data structure: 2 hashes for tRNA gff files
# %gff_cols:     key is start coord 
# %gff_cols_end: key is end coord
#############################################

my ($chrom, $id, $start, $end, $type, $anticodon, $istart, $iend, $strand, $junk, 
    $transcript, %gff_cols, %gff_cols_end, %predict_cols, %predict_cols_end);

# use latest WS release version
my $current = `grep "NAME WS" $curr_db/wspec/database.wrm`; $current =~ s/NAME //; chomp $current;

# use also same version of gff file for tRNA
my $gfffile = "/wormsrv2/autoace/GFF_SPLITS/$current/CHROMOSOME_*.exon_tRNA.gff";


# fetch only some columns in tRNA gff file
my @gffcols = `cut -f 1,4,5,7,9 $gfffile`;
foreach(@gffcols){
  chomp;
 # print $_, "\n";
  if ($_ =~ /^CHR.+/){
    ($chrom, $start, $end, $strand, $junk, $transcript) = split(/\s+/, $_);
    $chrom =~ s/CH.+_//;
    $transcript =~ s/Transcript|\"//g; 
 
    # swap start with end coord in gff file if found in "-" strand
    if ($strand eq "-"){
      push(@{$gff_cols{$end}}, $start, $transcript, $chrom, $strand);
      push(@{$gff_cols_end{$start}}, $end, $transcript, $chrom, $strand);
    }
    else {
      push(@{$gff_cols{$start}}, $end, $transcript, $chrom, $strand);
      push(@{$gff_cols_end{$end}}, $start, $transcript, $chrom, $strand);
     }
  } 
}

#foreach (sort keys %gff_cols){
#    print "$_ -> @{$gff_cols{$_}}\n";
#}

####################################################
# data structure: 2 hashes for predicted tRNA output
# %predict_cols:     key is start coord 
# %predict_cols_end: key is end coord
####################################################

my @predict_cols = `cut -f 1,2,3,4,5,6,7,8,9 $file`;

# fetch columns from tRNAScan output
foreach (@predict_cols){
  chomp;
  if ($_ =~ /^CHR.+/){
    ($chrom, $id, $start, $end, $type, $anticodon, $istart, $iend, ) = split(/\s+/, $_);
    $chrom =~ s/CH.+_//;

    if ($istart !=0 || $iend != 0){
      push(@{$predict_cols{$start}}, $istart-1, $type, $anticodon, $chrom, $id);
      push(@{$predict_cols{$iend+1}}, $end, $type, $anticodon, $chrom, $id);
      push(@{$predict_cols_end{$end}}, $iend+1, $type, $anticodon, $chrom, $id);
    }
    else {
      push(@{$predict_cols{$start}}, $end, $type, $anticodon, $chrom, $id);
      push(@{$predict_cols_end{$end}}, $start, $type, $anticodon, $chrom, $id);
    }
  }
}

#foreach (sort keys %predict_cols){
#  print "$_ -> @{$predict_cols{$_}}\n";
#}

#################################################################
# comparing tRNA coordinates of existing ones with predicted ones
#################################################################

my $end_diff = 0;
my $start_diff = 0;
my $same_start_end = 0;
my $same_end = 0;
my $pred_tRNA_id = 0;
my $diff_start = 0;
my $diff_end = 0;
my $diff_start_end = 0;
my $newp = 0;
my $pseudo = 0;
my ($g_chrom, $g_end, $p_chrom, $p_end); 

# loop through each predicted tRNA to compare their start and end coords with those of existing ones

foreach (sort keys %predict_cols){

  # ${@{$gff_cols{$_}}}[2]:     chrom from %gff_col
  # ${@{$predict_cols{$_}}}[3]: chrom from %predict_col

  # predicted tRNA with same start coord as in existing one on the same chromosome
  if (exists $gff_cols{$_} && exists ${@{$gff_cols{$_}}}[2] && exists ${@{$predict_cols{$_}}}[3] 
    && ${@{$gff_cols{$_}}}[2] eq ${@{$predict_cols{$_}}}[3] ){
    
    $g_end =   ${@{$gff_cols{$_}}}[0];       # end coord of existing tRNA
    $pred_tRNA_id = ${@{$predict_cols{$_}}}[4];
    $p_end =   ${@{$predict_cols{$_}}}[0];   # end coord of predicted tRNA
    #print "$g_end -> $p_end\n";
    
    # identical start and end coords
    if ($g_end == $p_end){
      $same_start_end++;
     # print "$transcript exists\n";
    }
    # same start coord, diff end coord
    if ($g_end != $p_end){
      $diff_end++;
     # print "Predicted tRNA #$pred_tRNA_id exists with diff end coord: $g_end -> $p_end\n";
    }
  }
  # predicted tRNAs with diff start coord
  else {

    $p_end = ${@{$predict_cols{$_}}}[0]; 

    # ${@{$gff_cols_end{$p_end}}}[2]: chrom from %gff_col
    # ${@{$predict_cols{$_}}}[3]: chrom from %predict_col

    $pred_tRNA_id = ${@{$predict_cols{$_}}}[4];
    $type = ${@{$predict_cols{$_}}}[1];
    

    # diff start coord and end coord, but not pseudogene
    if (!exists $gff_cols_end{$p_end} && $type ne "Pseudo"){ 
      $diff_start_end++;
      print "(P:$diff_start_end) Diff start & end codon: $_ -> @{$predict_cols{$_}}\n";
    }
    # diff start coord and end coord, and found to be pseudogene
    if (!exists $gff_cols_end{$p_end} && $type eq "Pseudo"){
      $pseudo++;
      print "(Pp:$pseudo) Diff start & end codon: $_ -> @{$predict_cols{$_}}\n";
    }  
    # diff start coord, but same end coord on the same chromosome
    if (exists $gff_cols_end{$p_end} && exists  ${@{$gff_cols_end{$p_end}}}[2] && exists ${@{$predict_cols{$_}}}[3]
        && ${@{$gff_cols_end{$p_end}}}[2] eq ${@{$predict_cols{$_}}}[3] ){ 
      $same_end++;
      print "Diff start codon: $pred_tRNA_id($_) -> E(${@{$gff_cols_end{$p_end}}}[0]), SAME end codon $p_end\n"; 
    
    }
  }
}

###########################################
# Existing tRNA not found in new prediction
###########################################

my $problem = 0;

print "\nExisting tRNAs having diff start AND end coords from predicted\n";
foreach (sort keys %gff_cols){
  $start = $_;
  $end =  ${@{$gff_cols{$_}}}[0];
  $chrom =  ${@{$gff_cols{$_}}}[2];
  $strand = ${@{$gff_cols{$_}}}[3];
  
  if (!exists $predict_cols{$start} && !exists $predict_cols_end{$end}){
    $problem++;
    print "   $start -> @{$gff_cols{$start}}\n"; 
  }
}

######################################################
# simple comp stats about predicted and existing tRNAs
######################################################

print "\nComparing existing tRNA: ", scalar (keys %gff_cols), " with predicted tRNA: ", scalar (keys %predict_cols), "\n";
print "\nIn prediction:\n$same_start_end tRNA with identical coords as existing ones already exist\n";
print "$diff_end with diff end coord but identical start coord: update existing tRNA\n";
print "$same_end with diff start coord, but identical end coord: update existing tRNA\n";
print $diff_start_end+$pseudo, " with diff start/end coords: among them $pseudo are Pseudos\n";
print "\nIn existing tRNAs:\n$problem have start & end coords diff from new prediction: check by hand\n";

__END__


=head2 NAME - parse_predicted_tRNA.pl  

=head2 DESCRIPTION

            This script parses the output of tRNA prediction suite tRNASCAN SE and write ace file 
            for predicted tRNA or existing tRNA objects.
         
=head3 <USAGE> 
 

=head2 Options: [h or help] [f or file]


B<-help:>     
            Displays documentation

B<-file:>   
            specify file generated by tRNASCAN SE


