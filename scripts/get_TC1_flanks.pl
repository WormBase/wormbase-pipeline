#!/usr/local/bin/perl5.8.0 -w

# Author: Chao-Kung Chen
# Last updated by $Author: ck1 $
# Last updated on: $Date: 2004-02-13 17:35:32 $ 

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'}; 
use Wormbase;
use Ace;
use lib "/nfs/team71/worm/ck1/WORMBASE_CVS/scripts/";
use Geneace;
use Coords_converter;



# ------global variables

my (%WB_TC1_seqs, @tc1_BG, %Tc1_PSL, %Tc1_chrom, %WB_TC1, %Tc1_OK, %WB_TC1_chrom);

my $curr_db = "/nfs/disk100/wormpub/DATABASES/current_DB";
my $coords = Coords_converter->invoke("$curr_db");
my $output_dir = "/nfs/team71/worm/ck1/TMP";

# -----start working

# parse TC1 sequence
my @tc1_seqs = `cat $output_dir/tc1.dna`;
my ($Tc1,$seq);

foreach (@tc1_seqs){
  chomp;
  if ($_ =~ /^\>(TC1:.+)/ ){$Tc1 = $1}
  if ($_ =~ /([actgn]+)/){$seq = $1; $WB_TC1_seqs{$Tc1} .= $seq}
}

open(TMP, ">$output_dir/tc1_temp.fasta") || die $!;;
open(ACE, ">$output_dir/tc1.ace") || die $!;
&check_unique_genomic_from_Bargmann_lab;
&get_equal_blat_chrom_span;


#######################################################
#               s u b r o u t i n e s
#######################################################

sub check_unique_genomic_from_Bargmann_lab {

  #  ----- parse Bargmann seq. from http://bargmann.ucsf.edu/~smccarroll/tc1locations.html
  #  ----- check that genomic seq. (200bp) after TC1 insertion site from Bargmann lab 
  #  ----- matches 100% to genome without ambiguity
  #  ----- create allele obj. based on them (some TC1 seq. are present in Wormbase)

  my %TC1_chrom_genomic;
  my @TC1_file = `cut -f 1,4,6 $output_dir/Bargmann_TC1`;
  
  foreach (@TC1_file){
    chomp;
    my ($chrom, $tc1, $genomic) = split(/\s+/, $_);
    push (@{$TC1_chrom_genomic{$tc1}}, $chrom, $genomic);
  }
  
#  # ----- convert Bargmann file into FASTA format for BLAT
#  my $tc1_query = "$output_dir/Bargmann_TC1_blat_query";
#  open(BLAT, ">$tc1_query") || die $!;

#  foreach (keys %TC1_chrom_genomic){
#    print  BLAT ">$_ -> $TC1_chrom_genomic{$_}->[0]\n";
#    print  BLAT "$TC1_chrom_genomic{$_}->[1]\n\n";
#  }
#  close BLAT;
 
#  # -----BLATting
                              
#  my $blat = "/nfs/disk100/wormpub/blat/blat";                          # blat binary
#  my $DNAs = "/nfs/team71/worm/ck1/SEQUENCES/CHROMOSOME_all.dna";       # sequence input for blat
  my $psl_file = "$output_dir/Bargman_TC1.psl";  # BLAT psl output file 

#  # run blat and output a psl file as $output
#  #`$blat -noHead $DNAs $tc1_query $psl_file`;
#  `$blat $DNAs $tc1_query $psl_file`;


  # -----Check BLAT result of Bargmann sequences
  my $count = 0;
  my @psl = `cut -f 1,2,9,10,14,16,17 $psl_file`;
  foreach (@psl){
    chomp;
    if ($_ =~ /^\d/){ # omit header
      my ($match, $mismatch, $strand, $tc1, $chrom, $start, $end) = split(/\s+/, $_);
      if ($mismatch == 0 && $match == length $TC1_chrom_genomic{$tc1}->[1]){
        $count++;
        $chrom =~ s/CHROMOSOME_//;

       # print "$match, $mismatch, $tc1, $chrom\t";
       # print length $TC1_chrom_genomic{$tc1}->[1], "***\n";
        if ( $tc1 =~ /pk(\d+)_1_0_Tc1/ ){$tc1 = "TC1:PK".$1."_L"}
        if ( $tc1 =~ /pk(\d+)_0_1_Tc1/ ){$tc1 = "TC1:PK".$1."_R"}
       # print $tc1, "#####\n";
        push(@tc1_BG, $tc1);
        push(@{$Tc1_PSL{$tc1}}, $strand, $chrom, $start, $end);
        push(@{$Tc1_chrom{$chrom}}, $tc1);
      }
    }
  }
  print $count, "\n";
  &get_30_bp_flanks();  # for all Bargmann TC1 sequences matched 100% to genomic sequence
}

sub get_equal_blat_chrom_span {  
  
  my $ga = init Geneace();
  my $tace = &tace;

  #-------------------------------------
  #  grep the length of TC1 insertion 
  #-------------------------------------
  my $db = Ace->connect(-path  => $curr_db,
			-program =>$tace) || do { print "Connection failure: ",Ace->error; die();};
  
  my $tc1_query = "find Sequence TC1*";
  my %tc1_leng;
  
  push( my @TC1s, $db->find($tc1_query) );
  foreach (@TC1s){
    my $tc1_seq = $_ -> DNA(1);
    my $tc1_length = $_ -> DNA(2);
    $tc1_leng{$tc1_seq} = $tc1_length;
  }
  
  my $gff_dir = "/wormsrv2/autoace/GFF_SPLITS/WS118";
  
  open(IN, "$output_dir/tc1_BLAT_BEST") || die $!;
  
  my (%TC1_blat_LR, %BC, %tc1_BG, %Tc1_to_check, %equal_span);

  foreach (@tc1_BG){$tc1_BG{$_}++} # array lookup
  my %seen;
  
  while(<IN>){
    chomp;
    if ($_ =~ /.+:CHROMOSOME_(.+)\s+(\d+)\s+(\d+)\s+Target\s+\".+:(TC1:PK\d+_.+)\"\s+(\d+)\s(\d+)$/){
      my $chrom = $1;
      my $L = $2;
      my $R = $3;
      my $TC1 = $4;
      my $blat_L = $5;
      my $blat_R = $6;
#    print "$chrom $L $R $blat_L $blat_R  #####\n";
      if (exists $tc1_BG{$TC1}){
	$Tc1_OK{$TC1} = $TC1;
      }
      else {
	$Tc1_to_check{$TC1} = $TC1;
      }

      push(@{$WB_TC1{$chrom}},$TC1) if !exists $seen{$TC1}; # record same TC1 only once

      $WB_TC1_chrom{$TC1} = $chrom;
      
      $seen{$TC1}=$TC1; 

      push(@{$TC1_blat_LR{$TC1}{'blat_L'}}, "$blat_L:$L");
      push(@{$TC1_blat_LR{$TC1}{'blat_R'}}, "$blat_R:$R");
    }
  }
 
  my $count_OK = 0;
  my $count_match =0;
  my $count_error =0; 
 
  foreach my $tc1 (keys %TC1_blat_LR){
    my ($chrom_coord, $blat_coord, @blats);
    
    push(my @coords, @{$TC1_blat_LR{$tc1}{'blat_L'}}, @{$TC1_blat_LR{$tc1}{'blat_R'}});
    foreach (@coords){
      ($blat_coord, $chrom_coord) = split(/:/, $_);
      $BC{$blat_coord} = $chrom_coord;   # %BC key: blat coord, value: corresponding chrom. coord
      push(@blats, $blat_coord);
    }
    
    @blats = sort {$a<=>$b} @blats;
    my $blat_length = $BC{$blats[-1]} - $BC{$blats[0]}  + 1 if $BC{$blats[-1]} >  $BC{$blats[0]};
       $blat_length = $BC{$blats[0]}  - $BC{$blats[-1]} + 1 if $BC{$blats[0]}  >  $BC{$blats[-1]};

    #if ( ($tc1_leng{$tc1} == $blat_length or $tc1_leng{$tc1} < $blat_length) && exists $Tc1_to_check{$tc1} ){
    if ( ($tc1_leng{$tc1} == $blat_length) && exists $Tc1_to_check{$tc1} ){  
      push(@{$equal_span{$tc1}}, $BC{$blats[0]}, $BC{$blats[-1]} ); # values are chrom. coords
     # print $tc1, "////\n";  OK
    }
    else {
      blast_TC1_seq($tc1) if $tc1 eq "TC1:PK5000_L";
    } 
  }

  undef %seen;# =(); # undef
  get_30_bp_flanks(%equal_span);
  
  #-----------------------------------------------------------------------
  #   little stats comparing Wormbase TC1 and those from Bargmann lab
  #-----------------------------------------------------------------------
  my @new_Tc1;
  foreach (keys %tc1_BG){
    if (!exists  $Tc1_to_check{$_} && !exists  $Tc1_OK{$_} ){
      push(@new_Tc1, $_);
    }
  }
  
  print scalar (keys %Tc1_OK), " OK, ", scalar (keys %Tc1_to_check), " to check in WB:\n";
  my @checks = keys %Tc1_to_check;
  print "@checks\n";
  print scalar @new_Tc1, " new TC1 from Bargmann lab:\n  @new_Tc1\n";
}

sub blast_TC1_seq {
  my $tc1 = shift;
  my $input = (); $input = ">$tc1\n".$WB_TC1_seqs{$tc1};
#  print $input, "###\n";
  print TMP $input;
  my $chrom = $WB_TC1_chrom{$tc1};
  
  `bl2seq -i $output_dir/tc1_temp.fasta -g -j /nfs/disk100/wormpub/DATABASES/current_DB/CHROMOSOMES/CHROMOSOME_$chrom.dna -p blastn -o $output_dir/BLAST/tc1_blast_output.$tc1`;
}

sub get_30_bp_flanks {
  
  my (%WBTc1OK) = @_;
  my @LGs = qw (I II III IV V X);
  
  foreach my $chrom (@LGs){
    
    # -- prepare DNA files
    my $dna_file = "$curr_db/CHROMOSOMES/CHROMOSOME_".$chrom.".dna";
    my @line = `egrep "[atcg]" $dna_file`;
    my $line;

    foreach (@line){chomp; $line .= $_}
 
    if (!%WBTc1OK){
      foreach ( @{$Tc1_chrom{$chrom}} ){
	
	my $strand = $Tc1_PSL{$_} -> [0];
	my $start  = $Tc1_PSL{$_} -> [2];
	my $end    = $Tc1_PSL{$_} -> [3];

	ace_output($_, $strand, $start, $end, $line, $chrom); 
      }
    }
    else {
      foreach ( @{$WB_TC1{$chrom}} ){
	if ( exists $WBTc1OK{$_} ){
	  
	  my $start = $WBTc1OK{$_}->[0];
	  my $end   = $WBTc1OK{$_}->[1];
	  my $strand;
	 # print "$_, $start, $end (3)\n"; # if $R < $L; 
	  if ($start > $end){
	    $strand ="-";
	  }
	  else{
	    $strand = "+";
	  }
	  ace_output($_, $strand, $start, $end, $line, $chrom); 
	}
      }
    }
  }
}


sub ace_output {
  my ($tc1, $strand, $start, $end, $line, $chrom) =@_;

  $tc1 =~ /TC1:PK(\d+)_\w/;
  my $allele = "pkP".$1;

  if ($strand eq "+"){

    my $DNA_L = substr($line, $start-31, 30);
    my $DNA_R = substr($line, $start+1, 30);
    my @clone = $coords->LocateSpan("CHROMOSOME_$chrom",$start, $start);
    
    print ACE "\nAllele : \"$allele\"\n";
    print ACE "Evidence PMID_evidence \"8962114\"\n";
    print ACE "Insertion\n";
    print ACE "Location \"NL\"\n";
    print ACE "Species \"Caenorhabditis elegans\"\n";
    print ACE "Sequence \"$clone[0]\"\n";
    print ACE "Flanking_sequences \"$DNA_L\" \"$DNA_R\"\n";
    print ACE "Transposon_insertion \"Tc1\"\n";
    print ACE "Method Transposon_insertion\n";
    print ACE "Remark \"This insertion site was identified using sequence from a shotgun library of Tc1 flanks.  It should be confirmed before it is studied further\" PMID_evidence \"8962114\"\n"; 
  }
  else {
    
    my $DNA_L = substr($line, $end+1, 30);
    $DNA_L = reverse $DNA_L; $DNA_L =~ tr/atcg/tagc/;
    my $DNA_R = substr($line, $end-31, 30);
    $DNA_R = reverse $DNA_R; $DNA_R =~ tr/atcg/tagc/;
    my @clone = $coords->LocateSpan("CHROMOSOME_$chrom",$end, $end);

    print ACE "\nAllele : \"$allele\"  // REV\n"; 
    print ACE "Evidence PMID_evidence \"8962114\"\n";
    print ACE "Insertion\n";
    print ACE "Location \"NL\"\n";
    print ACE "Species \"Caenorhabditis elegans\"\n";
    print ACE "Sequence \"$clone[0]\"\n";
    print ACE "Flanking_sequences \"$DNA_L\" \"$DNA_R\"\n";
    print ACE "Transposon_insertion \"Tc1\"\n";
    print ACE "Method Transposon_insertion\n";
    print ACE "Remark \"This insertion site was identified using sequence from a shotgun library of Tc1 flanks.  It should be confirmed before it is studied further\" PMID_evidence \"8962114\"\n";
  }
}
