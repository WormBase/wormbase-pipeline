#!/usr/local/bin/perl5.8.0 -w

# get_allele_flank_seq.pl

# Author: Chao-Kung Chen
 

# Retrieve 30 bp to the right and left of a specified amino acid

use strict;
use lib "/wormsrv2/scripts/"; 
use Wormbase;
use Getopt::Long;

###################################################
# variables and command-line options with aliases # 
###################################################

my ($aa, $dna, $cds);

GetOptions (
	    "cds=s" =>  \$cds,
	    "aa=s"  =>  \$aa,
            "dna=s" =>  \$dna
           );

open(IN1, "/nfs/disk100/wormpub/WORMPEP/wormpep_current") || die $!; 
#open(IN2, "/nfs/disk100/wormpub/WORMPEP/wormpep.dna_current") || die $!;
my $dna_file = "/wormsrv1/chaokung/my-scripts/$cds.dna";
open(IN2, $dna_file) || die $!;

my ($prot_seq, $DNA_seq, @DNA, @prot);  

if ($aa){
  $prot_seq = get_seq($cds, *IN1);
#  print "$cds\n$prot_seq\n\n";   
  @prot = split(//, $prot_seq);
  print "aa $aa = $prot[$aa-1] [length = ", scalar @prot, "]\n";
}

if ($dna){
  $DNA_seq = get_seq($cds, *IN2);
#  print "$cds\n$DNA_seq\n\n";   
  @DNA = split(//, $DNA_seq);
  print "DNA [length = ", scalar @DNA, "]\n";
  print "position $dna = $DNA[$dna-1]\n";
}

#####################################
# check for latest exon table version
#####################################

my (@current, $current, $archive);
my @version = `ls /wormsrv2/autoace/GFF_SPLITS/`;
while (<@version>){chomp; if ($_ =~ /^WS(\d+)/){push (@current, $1)}}
@current = sort {$a<=>$b} @current;

my @archive = `ls /wormsrv1/geneace/ALLELE_DATA/EXON_TABLES/`;
while (<@archive>){chomp; if ($_ =~ /.+(WS\d+)/){$archive = $1}}
if ("WS$current[-1]" ne $archive) {
  print "Update exon tables of all CDS/Transcripts\n";
  `echo "table-maker -o /wormsrv1/geneace/ALLELE_DATA/EXON_TABLES/CDS_exons_WS$current[-1] -p /wormsrv1/geneace/wquery/get_CDS_source_exons.def" |
   tace /nfs/disk100/wormpub/DATABASES/current_DB`;
}

my @exons = `grep $cds /wormsrv1/geneace/ALLELE_DATA/EXON_TABLES/*`;  


exons_to_codons($cds, \@exons, \@DNA, \@prot);

#######################
# s u b r o u t i n e s
#######################

sub exons_to_codons {
  my ($cds, $exons) = @_;
  my ($i, $j, $start, $end, @exon_start_end, $num_aa, $remainder, %codon_seq, $codon_seq, $ref, %remainder_hash, $start_bp, @return, @all);
  my $aa_length = 0;
  my $aa_codon = 0;
    
  foreach (@$exons){
    chomp;
    my ($cds, $start, $end) = split (/\s+/, $_);
    push (@exon_start_end, $start, $end);	         
  }

  # verify sequence 
  for ($i = 0; $i < scalar @exon_start_end; $i=$i+2){
    if ($i == 0){
  #    print $exon_start_end[$i],">>>>\n";
      $remainder = ($exon_start_end[$i+1]-$exon_start_end[$i]+1)%3;
      $remainder_hash{$i} = $remainder;
      $num_aa = ($exon_start_end[$i+1]-$exon_start_end[$i]+1-$remainder)/3;
  #    print "$num_aa -> $remainder #\n";
    }
    if ($i > 0 && $remainder_hash{$i-2} == 1){
      $remainder = ($exon_start_end[$i+1]-$exon_start_end[$i]-1)%3;    
      $remainder_hash{$i} = $remainder;
      $num_aa = ($exon_start_end[$i+1]-$exon_start_end[$i]-1-$remainder)/3;
    }
    if ($i > 0 && $remainder_hash{$i-2} == 2){
      $remainder = ($exon_start_end[$i+1]-$exon_start_end[$i])%3;  
      $remainder_hash{$i} = $remainder;
      $num_aa = ($exon_start_end[$i+1]-$exon_start_end[$i]-$remainder)/3;
    }
    if ($i > 0 && $remainder_hash{$i-2} == 0){
      $remainder = ($exon_start_end[$i+1]-$exon_start_end[$i]+1)%3;  
      $remainder_hash{$i} = $remainder;
      $num_aa = ($exon_start_end[$i+1]-$exon_start_end[$i]+1-$remainder)/3;
    }
       
    $aa_length += $num_aa;
#    print $aa_length, "\n";
    $start_bp = $exon_start_end[$i];
#    print $remainder_hash{$i},"\n";

   ########################################################################################################################
    if ($remainder == 0){
      if ($i == 0){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]+1;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "A", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq);
      }
      if ($i > 0 && $remainder_hash{$i-2} == 1 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-1;
         @return = codon_to_seq($start_bp, $num_bp, \@DNA, "B", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq);
      }
      if ($i > 0 && $remainder_hash{$i-2} == 2 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i];
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "C", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq);
      }   
      if ($i > 0 && $remainder_hash{$i-2} == 0 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]+1;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "A", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq);
      }
    }
    #########################################################################################################################
    if ($remainder == 1){
      if ($i == 0){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i];
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "D", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq);
      }
      if ($i > 0 && $remainder_hash{$i-2} == 1 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-2;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "H", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq);
      }
      if ($i > 0 && $remainder_hash{$i-2} == 2 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-1;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "C", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq);
      }
      if ($i > 0 && $remainder_hash{$i-2} == 0 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i];
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "D", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq);
      }
    }
    ######################################################################################################################
    if ($remainder == 2){
      if ($i == 0){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-1;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "E", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq);
      }
      if ($i > 0 && $remainder_hash{$i-2} == 1 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-3;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "F", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq);
      }
      if ($i > 0 && $remainder_hash{$i-2} == 2 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-2;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "G", $aa_codon, $i, @exon_start_end);
       %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq);
      }
      if ($i > 0 && $remainder_hash{$i-2} == 0 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-1;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "E", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq);
      }
    }
    ######################################################################################################################
  }

   my $codons = scalar @all;
#   print $codons, ">>>\n";
   for($i=0; $i < scalar @all; $i=$i+2){
     push(@{$codon_seq{$all[$i]}}, ${@{$all[$i+1]}}[0], ${@{$all[$i+1]}}[1], ${@{$all[$i+1]}}[2], ${@{$all[$i+1]}}[3], ${@{$all[$i+1]}}[4], ${@{$all[$i+1]}}[5]);
#     print "${@{$all[$i+1]}}[0], ${@{$all[$i+1]}}[1], ${@{$all[$i+1]}}[2], ${@{$all[$i+1]}}[3], ${@{$all[$i+1]}}[4], ${@{$all[$i+1]}}[5]\n";
   }  

  print "$prot[$aa-1] = @{$codon_seq{$aa}}\n";
  print "$prot[$aa-1]($aa) = @{$codon_seq{$aa}}->[0] @{$codon_seq{$aa}}->[2] @{$codon_seq{$aa}}->[4]\n";


  if ($aa > 11 ){
    print @DNA[@{$codon_seq{$aa}}->[1]-31..@{$codon_seq{$aa}}->[1]-2], " --- ", @DNA[@{$codon_seq{$aa}}->[5]..@{$codon_seq{$aa}}->[5]+29],"\n";
  }
  else {print "Left flanking sequences out of range -- ", @DNA[@{$codon_seq{$aa}}->[5]+1..@{$codon_seq{$aa}}->[5]+30],"\n"}
}

sub codon_to_seq {
  my ($start_bp, $num_bp, $DNA, $option, $aa_codon, $i, @exon_start_end) = @_;
  my %codon_seq;
 # print "### $aa_codon ###\n";
  #print "\$i = $i\n";
  for (my $j=0; $j < $num_bp-1; $j=$j+3){
    my $pos = $j + $start_bp;
    $aa_codon++;
    if ($option eq "A" || $option eq "D" || $option eq "E"){
#      print "Codon $aa_codon: ${@$DNA}[$pos-1] ", $pos, " ${@$DNA}[$pos] ", $pos+1, " ${@$DNA}[$pos+1] ", $pos+2,"\n";
      push(@{$codon_seq{$aa_codon}}, ${@$DNA}[$pos-1], $pos, ${@$DNA}[$pos], $pos+1, ${@$DNA}[$pos+1], $pos+2);
    }
    if ($option eq "B" || $option eq "H" || $option eq "F"){
#      print "Codon $aa_codon#: ${@$DNA}[$pos+1] ", $pos+2, " ${@$DNA}[$pos+2] ", $pos+3, " ${@$DNA}[$pos+3] ", $pos+4,"\n";
      push(@{$codon_seq{$aa_codon}}, ${@$DNA}[$pos+1], $pos+2, ${@$DNA}[$pos+2], $pos+3, ${@$DNA}[$pos+3], $pos+4);
    }
    if ($option eq "C" || $option eq "G"){
#      print "Codon $aa_codon: ${@$DNA}[$pos] ", $pos+1, " ${@$DNA}[$pos+1] ", $pos+2, " ${@$DNA}[$pos+2] ", $pos+3,"\n";
      push(@{$codon_seq{$aa_codon}}, ${@$DNA}[$pos], $pos+1, ${@$DNA}[$pos+1], $pos+2, ${@$DNA}[$pos+2], $pos+3);
    }  
  }
  
  if ($option eq "D" || $option eq "H" || $option eq "C"){
    $aa_codon++;
    my $pos_up = $exon_start_end[$i+1];
    my $pos_down = $exon_start_end[$i+2];
   # print "\$i = $i\n";
#    print "Codon $aa_codon: ${@$DNA}[$pos_up-1], ", $pos_up, " ${@$DNA}[$pos_down-1], ", $pos_down, " ${@$DNA}[$pos_down], ", $pos_down+1, "\n";
    push(@{$codon_seq{$aa_codon}}, ${@$DNA}[$pos_up-1], $pos_up, ${@$DNA}[$pos_down-1], $pos_down, ${@$DNA}[$pos_down], $pos_down+1);
  }
  if ($option eq "E" || $option eq "F" || $option eq "G"){
    $aa_codon++;
    my $pos_up = $exon_start_end[$i+1];
    my $pos_down = $exon_start_end[$i+2];
#    print "Codon $aa_codon: ${@$DNA}[$pos_up-2], ", $pos_up-1, " ${@$DNA}[$pos_up-1], ", $pos_up, " ${@$DNA}[$pos_down-1], ", $pos_down, "\n";
    push(@{$codon_seq{$aa_codon}}, ${@$DNA}[$pos_up-2], $pos_up-1, ${@$DNA}[$pos_up-1], $pos_up, ${@$DNA}[$pos_down-1], $pos_down);
  }
  
  return \%codon_seq, $aa_codon;
}

sub get_seq {

  my ($cds, $FH) = @_;
  my $seq =();
  my $line_count = 0;

  while(<$FH>){
    chomp;
    if ($_ =~ /^>$cds/){$line_count++}
    if ($_ !~ /^>/ && $line_count > 0){$seq .= $_; $line_count++}
    if ($_ =~ /^>/ && $line_count > 1){last}	
  }
  return $seq; 
}



