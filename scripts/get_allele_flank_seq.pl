#!/usr/local/bin/perl5.8.0 -w

# get_allele_flank_seq.pl 

# Date: 2003-05-28

# Author: Chao-Kung Chen

# Last updated by: $Author: ck1 $
# Last updated on: $Date: 2003-06-03 08:28:43 $

use strict;
use lib "/wormsrv2/scripts/"; 
use Wormbase;
use Term::ANSIColor;
use Getopt::Long;

###################################################
# variables and command-line options with aliases # 
###################################################

my ($aa, $dna, $cds, $help);

GetOptions (
	    "cds=s"  =>  \$cds,
	    "aa=s"   =>  \$aa,
            "dna"    =>  \$dna,
            "h|help" =>  \$help 
           );

if (!$cds || $help){system ("perldoc $0"); exit(0)}

my $mutation = $aa;  
$mutation =~ s/\d+//; $mutation = uc($mutation); # aa mutation
$aa =~ s/\D//;         
$cds = uc($cds);

#####################################
# check for latest exon table version
#####################################

my $tace = &tace;  
my $curr_db = "/nfs/disk100/wormpub/DATABASES/current_DB";
my $exon_tbl = "/wormsrv1/geneace/ALLELE_DATA/EXON_TABLES";

my ($current, @archive, $archive);
$current = `grep "NAME WS" $curr_db/wspec/database.wrm`; $current =~ s/NAME //; chomp $current;
@archive = `ls $exon_tbl`; foreach (@archive){chomp; if ($_ =~ /CDS.+(WS\d+)/){$archive = $1}} 

if ("$current" ne "$archive") {	
  print "New WS release available . .\nFetching latest source exons of all CDS/Transcripts and 6 chromosomal DNA sequencs\n\n";
  system ("rm -f $exon_tbl/*");
  `echo "table-maker -o $exon_tbl/CDS_exons_$current -p /wormsrv1/geneace/wquery/get_CDS_source_exons.def" | $tace $curr_db`;  
   system ("cp -f /nfs/disk69/ftp/pub/wormbase/WS102/CHROMOSOMES/CHR*.dna.gz $exon_tbl/");
  system ("chmod 775 $exon_tbl/*; gunzip $exon_tbl/CHR*");
}

####################################################
# get DNA sequence (exon/intron) of a CDS/Transcript
####################################################
  
my $DNA;
chdir "/wormsrv2/autoace/GFF_SPLITS/$current/";

my $CDS_coords = `grep $cds *.genes.gff *.rna.gff | cut -f 1,4,5,7,9`; 
my ($chrom, $left, $right, $strand, $CDS)= split(/\s+/, $CDS_coords);
$chrom =~ s/.+CHROMOSOME_//;

# 30 bp extension beyond 1st/last nucleotide
$left -= 30;
$right += 30;

#my $dna_file = "/nfs/team71/worm/ck1/DNA/CHROMOSOME_".$chrom.".dna";
my $dna_file = "$exon_tbl/CHROMOSOME_".$chrom.".dna";
my @line = `egrep "[atcg]" $dna_file`;
my $line;
foreach (@line){chomp; $line .= $_}

if ($strand eq "-"){
  $DNA = substr($line, $left-1, $right-$left+1);
  $DNA = reverse $DNA; $DNA =~ tr/atcg/tagc/;
}
if ($strand eq "+"){
  $DNA = substr($line, $left-1, $right-$left+1);
}

my @DNA = split(//, $DNA);

########################################################
# get protein sequence (exon/intron) of a CDS/Transcript
########################################################

open(IN1, "/nfs/disk100/wormpub/WORMPEP/wormpep_current") || die $!; 

my ($prot_seq, $DNA_seq, @prot);  

$prot_seq = get_seq($cds, *IN1);
#  print "$cds\n$prot_seq\n\n";   
@prot = split(//, $prot_seq);
#  print "aa $aa = $prot[$aa-1] [length = ", scalar @prot, "]\n";

########################################
# fetch source exons of a CDS/Transcript
########################################

my $src_exon = "$exon_tbl/CDS*";
my @exons = `grep $cds $src_exon`;  
#print "@exons\n";

########################################################################################
# retrieving flank seq of a specified codon or mutation site via exons_to_codons routine
########################################################################################

exons_to_codons($cds, \@exons, \@DNA, \@prot);

#######################
# s u b r o u t i n e s
#######################

#############################################################################################
# This chunk does several things:
# 1. process soruce exon coods to figure out frame shift 
# 2. the result of 1 is passed into codon_to_seq routine to retrieves 30 bp flanks of a codon 
#############################################################################################

sub exons_to_codons {
  my ($cds, $exons) = @_;
  my ($i, $j, $start, $end, @exon_start_end, $num_aa, $remainder, %codon_seq, $codon_seq, $ref, %remainder_hash, $start_bp, @return, @all);
  my $aa_length = 0;
  my $aa_codon = 0;
    
  foreach (@$exons){
    chomp;
    my ($cds, $start, $end) = split (/\s+/, $_);
    $start += 30; $end += 30; # 30 bp extension to get flank seq. of 1st/last amino acid
    push (@exon_start_end, $start, $end);
    $start=(); $end=();	         
  }

  for ($i = 0; $i < scalar @exon_start_end; $i=$i+2){
    if ($i == 0 || ($i > 0 && $remainder_hash{$i-2} == 0)){
      $remainder = ($exon_start_end[$i+1]-$exon_start_end[$i]+1)%3;
      $remainder_hash{$i} = $remainder;
    }
    if ($i > 0 && $remainder_hash{$i-2} == 1){
      $remainder = ($exon_start_end[$i+1]-$exon_start_end[$i]-1)%3;    
      $remainder_hash{$i} = $remainder;
    }
    if ($i > 0 && $remainder_hash{$i-2} == 2){
      $remainder = ($exon_start_end[$i+1]-$exon_start_end[$i])%3;  
      $remainder_hash{$i} = $remainder;
    }
       
    $start_bp = $exon_start_end[$i];

    ########################################################################################################################
    if ($remainder == 0){
      if ($i == 0 || ($i > 0 && $remainder_hash{$i-2} == 0) ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]+1;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "A", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq); next;
      }
      if ($i > 0 && $remainder_hash{$i-2} == 1 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-1;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "B", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq); next;
      }
      if ($i > 0 && $remainder_hash{$i-2} == 2 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i];	
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "C", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq); next;
      }   
    }
    #########################################################################################################################
    if ($remainder == 1){
      if ($i == 0 || ($i > 0 && $remainder_hash{$i-2} == 0)){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i];
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "D", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq); next;
      }
      if ($i > 0 && $remainder_hash{$i-2} == 1 ){
	my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-2;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "H", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq); next;
      }
      if ($i > 0 && $remainder_hash{$i-2} == 2 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-1;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "C", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq); next;
      }
    }
    ######################################################################################################################
    if ($remainder == 2){
      if ($i == 0 || ($i > 0 && $remainder_hash{$i-2} == 0) ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-1;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "E", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq); next;
      }
      if ($i > 0 && $remainder_hash{$i-2} == 1 ){
        my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-3;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "F", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq); next;
      }
      if ($i > 0 && $remainder_hash{$i-2} == 2 ){
     	my $num_bp = $exon_start_end[$i+1]-$exon_start_end[$i]-2;
        @return = codon_to_seq($start_bp, $num_bp, \@DNA, "G", $aa_codon, $i, @exon_start_end);
        %codon_seq = %{$return[0]}; $aa_codon = $return[1];
        push(@all, %codon_seq); next;
      }
    }
    ######################################################################################################################
  }

  my $codons = scalar @all;
  for($i=0; $i < scalar @all; $i=$i+2){
    push(@{$codon_seq{$all[$i]}}, ${@{$all[$i+1]}}[0], ${@{$all[$i+1]}}[1], ${@{$all[$i+1]}}[2], ${@{$all[$i+1]}}[3], ${@{$all[$i+1]}}[4], ${@{$all[$i+1]}}[5]);
#     print "${@{$all[$i+1]}}[0], ${@{$all[$i+1]}}[1], ${@{$all[$i+1]}}[2], ${@{$all[$i+1]}}[3], ${@{$all[$i+1]}}[4], ${@{$all[$i+1]}}[5]\n";
  }  

  print "\n$prot[$aa-1]($aa) = "; 
  print "@{$codon_seq{$aa}}->[0] (", @{$codon_seq{$aa}}->[1]-30, ") @{$codon_seq{$aa}}->[2] (";
  print @{$codon_seq{$aa}}->[3]-30, ") @{$codon_seq{$aa}}->[4] (", @{$codon_seq{$aa}}->[5]-30, ") [full-length aa of this gene: ", scalar @prot, "]\n\n"; 
  my $codon = "@{$codon_seq{$aa}}->[0]"."@{$codon_seq{$aa}}->[2]"."@{$codon_seq{$aa}}->[4]";

  print "Triplet color coding: mutated site in ", color("red"), "RED", color("reset"), ", the other two sites behind or before it in the triplet in ", color("blue"), "BLUE", color("reset"), "\n"; 
  print "Flanking bp coding when frame shift:  upsteam of mutated site in ", color("magenta"), "MAGENTA", color("reset"), ", downstream in ", color("green"), "GREEN", color("reset"), "\n";
  print "This color coding does not apply to dinucleotide mutation sites. But flanking seqs are there anyway\n\n"; 

  ################################
  # output 30 bp flanks of a codon
  ################################

  # look up genetic code of aa 
  my %code = (
	  'A' => ["gct", "gcc", "gca", "gcg"],                     
	  'B' => ["aat", "aac", "gat", "gac"], 
	  'C' => ["tgt", "tgc"],
	  'D' => ["gat", "gac"], 
	  'E' => ["gaa", "gag"], 
	  'F' => ["ttt", "ttc"],
	  'G' => ["ggt", "ggc", "gga", "ggg"],
	  'H' => ["cat", "cac"],
	  'I' => ["att", "atc", "ata"],
	  'K' => ["aaa", "aag"],
	  'L' => ["tta", "ttg", "ctt", "ctc", "cta", "ctg"], 
	  'M' => ["atg"],
	  'N' => ["aat", "aac"], 
	  'P' => ["cct", "ccc", "cca", "ccg"], 
	  'Q' => ["caa", "cag"],
	  'R' => ["cgt", "cgc", "cga", "cgg", "aga", "agg"], 
	  'S' => ["tct", "tcc", "tca", "tcg", "agt", "agc"],
	  'T' => ["act", "acc", "aca", "acg"], 
	  'V' => ["gtt", "gtc", "gta", "gtg"],
	  'W' => ["tgg"], 
	  'Y' => ["tat", "tac"], 
	  'Z' => ["gaa", "gag", "caa", "cag"],
	  'X' => ["taa", "tag", "tga"], #stop codons
  );

  print "-------------------------------------\n";
  print "   	Codon      ($prot[$aa-1]):\t\t$codon\n";
  for ($i=0; $i < scalar @{$code{$mutation}}; $i++){
    print "	Mutated to \($mutation\):\t\t@{$code{$mutation}}->[$i-1]\n" if $mutation ne "X";
    print "	Mutated to \(STOP\):\t@{$code{$mutation}}->[$i-1]\n" if $mutation eq "X";
  }
  print "-------------------------------------\n";    

  my ($first_bp, $second_bp, $third_bp, $first_site, $second_site, $third_site);
  $first_bp = @{$codon_seq{$aa}}->[1];
  $second_bp = @{$codon_seq{$aa}}->[3];
  $third_bp = @{$codon_seq{$aa}}->[5];
  $first_site = @{$codon_seq{$aa}}->[0];
  $second_site = @{$codon_seq{$aa}}->[2];
  $third_site = @{$codon_seq{$aa}}->[4];

  #######################################################################
  # output 30 bp flank seq under frame shift or no frame shift situations
  #######################################################################

  ################
  # no frame shift
  ################

  if ($first_bp == $second_bp-1 && $second_bp == $third_bp-1){
    for($i=1; $i<4; $i++){
      print "\n\# 1st site mutation:\n" if $i == 1; print "\n\# 2nd site mutation:\n" if $i == 2; print "\n\# 3rd site mutation:\n" if $i == 3;
      print @DNA[$first_bp-32+$i..$first_bp-2];
      print color("red"), " $first_site ", color("reset"), color("blue"), $second_site.$third_site, color("reset") if $i == 1;
      print color("blue"), $first_site, color("reset"), color("red"), " $second_site ", color("red"), color("reset"), color("blue"), $third_site, color("reset") if $i == 2;
      print color("blue"), $first_site.$second_site, color("reset"), color("red"), " $third_site ", color("reset") if $i == 3;
      print @DNA[$third_bp..$third_bp+26+$i], "\n";
    }
  }

  ######################  
  # second frame shifted
  ######################

  if ($first_bp != $second_bp-1 && $second_bp == $third_bp-1){

    print "\n\# 1st site mutation:\n";
    print color("green"), " $DNA[$first_bp]: ", $first_bp+1-30, color("reset"), "\n";
    print @DNA[$first_bp-31..$first_bp-2];
    print color("red"), " $first_site ", color("reset");
    print color("green"), $DNA[$first_bp], color("reset"), @DNA[$first_bp+1..$first_bp+29], "\n";   

    print "\n\# 2nd site mutation:\n"; 
    print color("magenta"), " $DNA[$second_bp-2]: ", $second_bp-1-30, color("reset"), "\n";;
    print @DNA[$second_bp-31..$second_bp-3], color("magenta"), $DNA[$second_bp-2], color("reset");
    print color("red"), " $second_site ", color("reset"), color("blue"), $third_site, color("reset"); 
    print @DNA[$third_bp..$third_bp+28], "\n";  

    print "\n\# 3rd site mutation:\n";
    print color("magenta"), " $DNA[$second_bp-2]: ", $second_bp-1-30, color("reset"), "\n";;
    print @DNA[$second_bp-30..$second_bp-3], color("magenta"), $DNA[$second_bp-2], color("reset");;   
    print color("blue"), $second_site, color("reset"), color("red"), " $third_site ", color("reset");                                                      
    print @DNA[$third_bp..$third_bp+29], "\n"; 
  }
 
 ######################
  # third frame shifted
  #####################

  if ($first_bp == $second_bp-1 && $second_bp != $third_bp-1){

    print "\n\# 1st site mutation:\n";
    print color("green"), " $DNA[$second_bp]: ", $second_bp+1-30, color("reset"), "\n";
    print @DNA[$first_bp-31..$first_bp-2];
    print color("red"), " $first_site ", color("reset"), color("blue"), $second_site, color("reset");
    print color("green"), $DNA[$second_bp], color("reset"), @DNA[$second_bp+1..$second_bp+28], "\n";
  
    print "\n\# 2nd site mutation:\n";
    print color("green"), " $DNA[$second_bp]: ", $second_bp+1-30, color("reset"), "\n";
    print @DNA[$first_bp-30..$first_bp-2];
    print color("blue"), $first_site, color("reset"), color("red"), " $second_site ", color("reset");  
    print color("green"), $DNA[$second_bp], color("reset"), @DNA[$second_bp+1..$second_bp+29], "\n";

    print "\n\# 3rd site mutation:\n";
    print color("magenta"), " $DNA[$third_bp-2]: ", $third_bp-1-30, color("reset"), "\n"; 
    print @DNA[$third_bp-31..$third_bp-3], color("magenta"), $DNA[$third_bp-2], color("reset");
    print color("red"), " $third_site ", color("reset");                                             
    print @DNA[$third_bp..$third_bp+29], "\n";
  }
  print "\n";
}

################################################################################################
# get DNA triplet of a specified amino acid based on source exons processed in the above routine
################################################################################################

sub codon_to_seq {
  my ($start_bp, $num_bp, $DNA, $option, $aa_codon, $i, @exon_start_end) = @_;
  my %codon_seq;

 # print "\$i = $i\n";
 # print "start: $start_bp end: $exon_start_end[$i+1] = ", $exon_start_end[$i+1] - $exon_start_end[$i]+1, " Bp = $num_bp","\n";

  for (my $j=0; $j < $num_bp; $j=$j+3){
    my $pos = $j + $start_bp;
    $aa_codon++;
    if ($option eq "A" || $option eq "D" || $option eq "E"){
     # print "Codon $aa_codon(ADE): ${@$DNA}[$pos-1] ", $pos, " ${@$DNA}[$pos] ", $pos+1, " ${@$DNA}[$pos+1] ", $pos+2,"\n";
      push(@{$codon_seq{$aa_codon}}, ${@$DNA}[$pos-1], $pos, ${@$DNA}[$pos], $pos+1, ${@$DNA}[$pos+1], $pos+2);
    }
    if ($option eq "B" || $option eq "H" || $option eq "F"){
     #  print "Codon $aa_codon(FH): ${@$DNA}[$pos+1] ", $pos+2, " ${@$DNA}[$pos+2] ", $pos+3, " ${@$DNA}[$pos+3] ", $pos+4,"\n";
       push(@{$codon_seq{$aa_codon}}, ${@$DNA}[$pos+1], $pos+2, ${@$DNA}[$pos+2], $pos+3, ${@$DNA}[$pos+3], $pos+4);
    }
    if ($option eq "C" || $option eq "G"){
     #print "Codon $aa_codon(CG): ${@$DNA}[$pos] ", $pos+1, " ${@$DNA}[$pos+1] ", $pos+2, " ${@$DNA}[$pos+2] ", $pos+3,"\n";
      push(@{$codon_seq{$aa_codon}}, ${@$DNA}[$pos], $pos+1, ${@$DNA}[$pos+1], $pos+2, ${@$DNA}[$pos+2], $pos+3);
    }  
  }

  if ($option eq "D" || $option eq "H"){
    $aa_codon++;
    my $pos_up = $exon_start_end[$i+1];
    my $pos_down = $exon_start_end[$i+2];
    # print "\$i = $i\n";
    #print "Codon $aa_codon(CDH): ${@$DNA}[$pos_up-1], ", $pos_up, " ${@$DNA}[$pos_down-1], ", $pos_down, " ${@$DNA}[$pos_down], ", $pos_down+1, "\n";
    push(@{$codon_seq{$aa_codon}}, ${@$DNA}[$pos_up-1], $pos_up, ${@$DNA}[$pos_down-1], $pos_down, ${@$DNA}[$pos_down], $pos_down+1);
  }
  if ($option eq "E" || $option eq "F" || $option eq "G"){
    $aa_codon++;
    my $pos_up = $exon_start_end[$i+1];
    my $pos_down = $exon_start_end[$i+2];
   # print "Codon $aa_codon(EFG): ${@$DNA}[$pos_up-2], ", $pos_up-1, " ${@$DNA}[$pos_up-1], ", $pos_up, " ${@$DNA}[$pos_down-1], ", $pos_down, "\n";
    push(@{$codon_seq{$aa_codon}}, ${@$DNA}[$pos_up-2], $pos_up-1, ${@$DNA}[$pos_up-1], $pos_up, ${@$DNA}[$pos_down-1], $pos_down);
  }
  
  return \%codon_seq, $aa_codon;
}

##########################################
# get protein seq via wormpep.current file
##########################################

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


__END__


=head2 NAME - get_allele_flank_seq.pl  

=head2 DESCRIPTION

            This script is suitable for curating allele flanking sequences described in paper such as "
            .... amino acid Q at position 235 is mutated to G ..."

            This script can output 30 bp flanking DNA seqs on two sides of any mutated site of a triplet of a defined gene (including the 
            first and last codon, although such alleles rarely occur). 

            You need to supply an amino acid coordinate followed by a mutation in single-letter code (case-insensitive)
            and a CDS/Transcript name (case-insensitive) as arguments (see USAGE) to retrieve the flanking sequences.  
            
	    Scenario:
                        An amino acid G has DNA triplet ggg at positions 590, 591 and 640 (3rd frame shifted),
                        and has a missense mutation to E. 
                        You want to retrieve the flanking seq. of the mutated sites in a triplet.
            
            Output:

                        (1) G(45) = g (590) g (591) g (640) [full-length aa of this gene = 255] 
            
			(2) Triplet color coding: mutated site in RED, the other two sites behind or before it in the triplet in BLUE
			    Flanking bp coding when frame shift:  upsteam of mutated site in MAGENTA, downstream in GREEN
			    This color coding does not apply to dinucleotide mutation sites. But flanking seqs are there anyway.

                        (3)
			-------------------------------------
			           Codon (G):    ggg
			     Mutation to (E):    gag
			     Mutation to (E):    gaa
			-------------------------------------

                        (4)
                         
                        # 1st site mutation:
			g: 592
			acagactacttaaacattgtaaaaggatat g ggtaagaatatatatcttatacaaccctta

			# 2nd site mutation:
			g: 592
			cagactacttaaacattgtaaaaggatatg g gtaagaatatatatcttatacaacccttac

			# 3rd site mutation:
			 g: 639 t: 641
			tacaacccttactgaattttaatttttcag g tgctactctcaagttggacgaactggagga


            Comments on output:
	
                       (1) Verification. This simply tells you that G(45) described in paper (scenario) is the same as current WS dataset. 
                           So should be OK to run script for this allele, . . . usually.  

                       (2) Color coding (cannot be seen here, but when you run the script) helps you quickly identify the flanking 
                           sequences (4) of a mutated site (the one in between spaces), especially in cases where frame shift occur so that 
                           the immediate flanking nucleotide maybe in the intron between two sites of a codon.  
                          
                       (3) As three potential single-site mutation can occur in a codon, (3) gives you genetic codes of the amino acid 
                           resulted in mutation and allows you a quick look up of bp substitution. This table is also fine 
                           for a dinucleotide mutation.
                           In the scenario, eg, if codon (G) has a missense mutation and changed to E, the codon table conveniently tells you 
                           that ggg has been mutated to gag. So, this would be a 2nd site mutation or a [g/a] substitution.
                           You should then choose the matching flanking sequences in (4).       

=head3 <USAGE> 
         
            Obligatory arguments -cds, -aa for retrieving 30 bp flank seqs
            Sample query:  perl get_allele_flank_seq.pl -cds 4R79.1 -aa 45E 


=head2 Options: [h or help] [cds] [aa]
             
          
B<-help:>     
            Displays this POD.

B<-cds:>   
            Specifies a CDS or Transcript, eg. -cds 4R79.1 
            (case insensitive)

B<-aa:>  
            Specifies an amino acid coordinate and mutation
            (letter, case insensitive). eg. -aa 45E 
            Type "X" for nonsense mutation, eg. -aa 45X
 
             

