#!/usr/local/bin/perl5.8.0 -w

# parse_allele_info_from_paper.pl 

# Date: 2004-02-27

# Author: Chao-Kung Chen

# Last updated by: $Author: ar2 $
# Last updated on: $Date: 2005-02-14 14:46:17 $

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'}; 
use Wormbase;
use Getopt::Long;

###################################################
# variables and command-line options with aliases # 
###################################################

my ($aa, $dna, $help, $verbose, $infile, $cds, $locus, $allele, $paper, $output);

GetOptions (
	    "i|input=s"  =>  \$infile,
	    "cds=s"      =>  \$cds,
	    "locus=s"    =>  \$locus,
	    "allele=s"   =>  \$allele,
	    "paper=s"    =>  \$paper,
	    "o|output=s" =>  \$output,
           );


my %RA = (
	  'a' => "Greenwald and Horvitz (1986)",
	  'b' => "[cgc1800]",
	  'c' => " Greenwald and Horvitz (1980)",
	  'd' => " I. Perez de la Cruz and H. R. Horvitz, unpublished results",
	  'e' => "[cgc1525]",
	  'f' => "J. Levin and H. R. Horvitz, unpublished results",
	  'g' => " I. Greenwald and H. R. Horvitz, unpublished results",

	 );


open(ACE, ">$output") || die $!;


# reading input file
my ($Allele, $ref, $WT_codon, $MU_codon, $WT_aa, $aa_pos, $MU_aa, $mutagen);

if ($infile){
  my @input = `cat $infile`;
  foreach (@input){
    chomp;

    # n345e TTA to TGA 	L20Ochre unc-93(e1500) 	NTG (format used in this paper)
    my($col1, $WT_codon, $col3, $MU_codon, $col5, $col6, $mutagen) = split(/\s+/, $_);
  
    # n345e (e is ref. code in this paper)
    if ($col1 =~ /(\w+\d+)(\w)/ ){
      $Allele = $1;
      $ref = $2;
      
    }
    if ($col5 =~ /([A-Z])(\d+)([a-zA-Z]+)/ ){
      $WT_aa  = uc($1);
      $aa_pos = $2;
      $MU_aa  = $3;
    }
    if ($col6 =~ /(\w{3,3}-\d+)\(\w+\d+\)/ ){
      $locus = $1;
    }					   

    print "$Allele, $cds, $ref, $WT_codon, $MU_codon, $WT_aa, $aa_pos, $MU_aa, $mutagen\n";
    my $MU_mis =();
    $MU_mis = "x" if $MU_aa eq "Amber" || $MU_aa eq "Ochre" || $MU_aa eq "Opal";
    my $param = "$aa_pos-$MU_aa-$WT_aa-$MU_codon-$WT_codon-$output-$Allele" if !$MU_mis;
    $param = "$aa_pos-x-$WT_aa-$MU_codon-$WT_codon-$output-$Allele" if $MU_mis;

   # $cds = "C46F11.1" if $locus eq "unc-93" or "sup-18";
#    $cds = "R09G11.1" if $locus eq "sup-10";
#    $cds = "F34D6.3" if $locus eq "sup-9";

    my $seq = $cds;
    $seq =~ s/\..+//;
    
    
    print ACE "\n\nVariation : \"$Allele\"\n";
    print ACE "Sequence \"$seq\"\n";
    print ACE "Mutagen \"$mutagen\"\n";
    
    print ACE "Amber_UAG\n" if $MU_aa eq "Amber";
    print ACE "Ochre_UAA\n" if $MU_aa eq "Ochre";
    print ACE "Opal_UGA\n"  if $MU_aa eq "Opal";

    print ACE "Missense\n" if $MU_mis ne "x";

    print ACE "Predicted_CDS \"$cds\"\n";
    print ACE "Species \"Caenorhabditis elegans\"\n";
    
    $MU_aa = uc($MU_aa) if length($MU_aa) == 1;

    if ($ref eq "b" or $ref eq "e"){
      print ACE "Remark \"$Allele is a $WT_aa($aa_pos) to $MU_aa mutation\" Paper_evidence \"$RA{$ref}\"\n" if $ref eq "b" or $ref eq "e";
    }
    else {
      print ACE "Remark \"$Allele is a $WT_aa($aa_pos) to $MU_aa mutation, $RA{$ref}\"\n"; 
    }
    print ACE "Method \"Substitution_allele\"\n";
    system("perl5.8.0 /nfs/team71/worm/ck1/WORMBASE_CVS/scripts/get_allele_flank_seq_multiple.pl -cds $cds -aa $param "); 
    
  }
}
