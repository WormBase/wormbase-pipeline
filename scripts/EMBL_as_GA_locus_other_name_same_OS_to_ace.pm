#!/usr/local/bin/perl5.6.0 -w

# What it does: Write .ace file (including NDB info) for EMBL names matching GA Locus or Other_name
# 

# Last updated by $Author: ck1 $
# Last updated on: $Date: 2002-09-12 13:59:09 $ 

package EMBL_as_GA_locus_other_name_same_OS_to_ace;

use Exporter;

@ISA = qw(Exporter);

@EXPORT = qw(Locus_Other_name); 

use lib "/wormsrv1/chaokung/my-scripts";
use EMBL_GeneAce;


###################################
# subroutine for writing .ace files
###################################

sub Locus_Other_name {

  my ($Embl_name, $Embl_OS, $DB_OS, $DB_locus, $DB_other, $peptide_status, $AC, $file, $DNA)=@_;
 # print "$DB_locus\n";

  my @Acc=split(/ /, $AC);       
  my @ACs;
  %seen=();
  foreach $e (@Acc){
    push(@ACs, $e) unless $seen{$e}++;
  }
  
  my (@LT, @EMBL_entry);

  if ($Embl_OS eq "C.elegans"){
    $Embl_OS = "Caenorhabditis elegans";
  }
  if ($DB_OS eq "C.elegans"){
    $DB_OS = "Caenorhabditis elegans";
  }
  if ($peptide_status eq "Live"){
    if ($DB_locus eq "NA"){
      ACE_NDB_C($file, $Embl_name, $DNA, $Embl_OS, @ACs);
    }
    else {
      ACE_NDB_A($file, $DB_locus, @ACs);
    }
  }  
}  
1;
 sub ACE_NDB_A {

    my ($file, $DB_locus, @ACs) = @_;

    #test
    #print $file, "\n";
    #print $DB_locus, "\n";
    #print "@ACs\n";

    my $acefile=$file."_EMBL_Gene_to_GeneAce.ace";
    open(OUT_EG, ">>/wormsrv1/chaokung/wormdata/ACE/$acefile") || die "Can't write to file!";
    print OUT_EG "Locus : \"$DB_locus\"\n";
    foreach $e (@ACs){
      print OUT_EG "Other_sequence\t\"$e\"\n";
    }
    foreach $e (@ACs){
      print OUT_EG "\nSequence : $e\n";
      print OUT_EG "DB_annotation\tEMBL\t\"$e\"\n\n";
    }
	
    foreach $e (@ACs){
      print OUT_EG "LongText : \"$e\"\n";

      ###########################################################
      # Subroutine @/wormsrv1/chaokung/my-scripts/EMBL_GeneAce.pm 
      ###########################################################

      my @LT=Delete_EMBL_DNA_Seq($e);
      foreach (@LT){
	print OUT_EG "$_\n";
      }
      print OUT_EG "***LongTextEnd***\n";
      print OUT_EG "\n";

      #####################################################################################
      # Get NDB info from subroutine get_NDB @/wormsrv1/chaokung/my-scripts/EMBL_GeneAce.pm
      #####################################################################################

      my @EMBL_entry=get_NDB(@LT);

      foreach (@EMBL_entry){
	print OUT_EG "$_";
      }
    } 
  }  
1;
 sub ACE_NDB_C {

    my ($file, $Embl_name, $DNA, $Embl_OS, @ACs) = @_;

    #test
    #print $file, "\n";
    #print $DB_locus, "\n";
    #print "@ACs\n";

    my $acefile=$file."_EMBL_Gene_to_GeneAce.ace";
    open(OUT_EG, ">>/wormsrv1/chaokung/wormdata/ACE/$acefile") || die "Can't write to file!";
    print OUT_EG "Locus : \"$Embl_name\"\n";
    print OUT_EG "Gene\n";
    print OUT_EG "Genomic_sequence\t$DNA\n"; 
    foreach $e (@ACs){
      print OUT_EG "Other_sequence\t$e\n";
    }   
    print OUT_EG "Species\t\"$Embl_OS\"\n";
    foreach $e (@ACs){
       print OUT_EG "$e, ";
    }
    print OUT_EG "[100802 ck1]\"\n";
    foreach $e (@ACs){
      print OUT_EG "\nSequence : $e\n";
      print OUT_EG "DB_annotation\tEMBL\t\"$e\"\n\n";
    }
    
    foreach $e (@ACs){
      print OUT_EG "LongText : \"$e\"\n";
      @LT=Delete_EMBL_DNA_Seq($e);
      foreach (@LT){
	print OUT_EG "$_\n";
      }
      print OUT_EG "***LongTextEnd***\n";
      print OUT_EG "\n";

      ######################################
      # Get NDB info from subroutine get_NDB
      ######################################

      @EMBL_entry=get_NDB(@LT);

      foreach (@EMBL_entry){
        print OUT_EG "$_";
      }
    }
  }
1;
