#!/usr/local/bin/perl5.6.1 -w
# 
# parse_predicted_tRNA.pl
#
# by Chao-Kung Chen
#
# Script to parse columns of predicted tRNA by tRNASCAN-SE
#
# Last updated by: $Author: ck1 $
# Last updated on: $Date: 2003-09-26 13:41:34 $

use strict;
use lib "/wormsrv2/scripts/"; 
use Wormbase;
use Getopt::Long;
use Ace;
use Coords_converter;


###########
# variables
###########

my $tace = &tace;  
my ($help, $file, $database);

GetOptions ("h|help"          => \$help,
            "f|file=s"        => \$file,
	    "db|database=s"   => \$database
	   );

if ($help){
  system("perldoc $0");
  exit(0);
}

# database used; default is current_db

if (!$database){$database = "/nfs/disk100/wormpub/DATABASES/current_DB"};

#############################################
# data structure: 2 hashes for tRNA gff files
# %gff_cols:     key is start coord 
# %gff_cols_end: key is end coord
#############################################

open(tRNA, ">tRNA_dataset") || die $!;

my ($chrom, $id, $start, $end, $strand, $junk, $transcript, %gff_cols, %gff_cols_end);

# use latest WS release version
my $current = `grep "NAME WS" $database/wspec/database.wrm`; $current =~ s/NAME //; chomp $current;

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

foreach (sort keys %gff_cols){
  print tRNA "$_ -> @{$gff_cols{$_}}\n";
}

####################################################
# data structure: 2 hashes for predicted tRNA output
# %predict_cols:     key is start coord 
# %predict_cols_end: key is end coord
####################################################

my ($type, $anticodon, $istart, $iend, %predict_cols, %predict_cols_end, %predict_cols_exon, %predict_cols_span);

my @predict_cols = `cut -f 1,2,3,4,5,6,7,8,9 $file`;

# fetch columns from tRNAScan-SE output
foreach (@predict_cols){
  chomp;
  if ($_ =~ /^CHR.+/){
    ($chrom, $id, $start, $end, $type, $anticodon, $istart, $iend, ) = split(/\s+/, $_);
    $chrom =~ s/CH.+_//;
    
    # tRNA with > 1 exons
    if ($istart !=0 || $iend != 0){
      # + strain
      if ($start < $end){
	push(@{$predict_cols_exon{$start}}, 1, $istart-$start);                   # source_exon 1 
	push(@{$predict_cols_exon{$iend+1}}, $iend+1-$start+1, $end-$start+1);    # source_exon 2 
	push(@{$predict_cols{$start}}, $istart-1, $type, $anticodon, $chrom, $id, $istart); 
	push(@{$predict_cols{$iend+1}}, $end, $type, $anticodon, $chrom, $id, $istart); 

        # start and end coords for clone coords of tRNA with multiple exons
	push(@{$predict_cols_span{$start}}, $start, $end); 
	push(@{$predict_cols_span{$iend+1}}, $start, $end); 
      }
      # - strain
      else {
	push(@{$predict_cols_exon{$start}}, 1, $start-$istart+2);                 # source_exon 1 
	push(@{$predict_cols_exon{$iend-1}}, $start-$iend, $start-$end+1);        # source_exon 2 
	push(@{$predict_cols{$start}}, $istart+1, $type, $anticodon, $chrom, $id, $istart); 
	push(@{$predict_cols{$iend-1}}, $end, $type, $anticodon, $chrom, $id, $istart); 

        # start and end coords for clone coords of tRNA with multiple exons
        push(@{$predict_cols_span{$start}}, $start, $end);
        push(@{$predict_cols_span{$iend-1}}, $start, $end);
      }
    }
    # tRNA with single exon
    else {
      push(@{$predict_cols{$start}}, $end, $type, $anticodon, $chrom, $id, $istart);
    }
  }
}

foreach (sort keys %predict_cols){
  print tRNA "$_ -> @{$predict_cols{$_}}\n";
}

foreach (keys %predict_cols_span){
  print tRNA "$_ -> @{$predict_cols_span{$_}}\n";
}

####################################################################
# retrieve all existing tRNAs types for comparing with predicted one
####################################################################

my @tRNA_type = `echo "table-maker -p /nfs/disk100/wormpub/DATABASES/current_DB/wquery/tRNA_type.def" | $tace $database`; 

my %tRNA_type;

foreach (@tRNA_type){
  chomp;
  if ($_ =~ /^\"(.+)\"\s+\"(.+)\"\s+\"(.+)\"$/){
    my $tRNA = $1;
    my $type = $2;
    my $lab  = $3;
    my ($codon, $aa_3, $aa_1) = split(/\s/, $type);

    # some tRNA do not have 2nd, 3rd fields
    if (defined $codon && !defined $aa_3 && !defined $aa_1){
      $aa_3 = $codon;
    }
    if (defined $aa_3 && ($aa_3 ne "???")){
      push(@{$tRNA_type{$tRNA}}, $aa_3, $lab);
    }
    else {
      push(@{$tRNA_type{$tRNA}}, "NA", $lab);
    }
  }
}

foreach (keys %tRNA_type){
  print tRNA "$_ -> @{$tRNA_type{$_}}###\n";
}

##################################
# hash key: clone, value: HX or RW
##################################

my @clone_to_lab = `echo "table-maker -p /nfs/disk100/wormpub/DATABASES/current_DB/wquery/genome_seq_to_lab.def" | $tace $database`;
my %clone_lab;
        
foreach (@clone_to_lab){
  chomp;
  if ($_ =~ /^\"(.+)\"\s+\"(.+)\"$/){
    my $clone = $1; my $lab = $2;
    $clone_lab{$clone} = $lab;
  }
}
foreach (keys %clone_lab){
  print tRNA "$_ -> $clone_lab{$_}    <clone_lab>\n";
}

##########################################################
# a hash of clones with last t# number from Transcrpt obj.
# used to assign next available seq. name for new tRNAs
##########################################################

my $curr_db = Ace->connect(-path => $database) || die "Connection failure: ", Ace->error;

my $query = "find Transcript where *.t* & !cb*"; 
my @seq_names;
push(@seq_names, $curr_db->find($query));
my %clone_last;

foreach (@seq_names){
  $_ =~ /(.+)\.t(\d+)/;
  my $clone = $1; my $last = $2;
  $clone_last{$clone} = $last;
}

foreach (keys %clone_last){
  print tRNA "$_ -> $clone_last{$_}\n";
}

##########################################################
# output acefile of existing tRNA with timestamp in a hash
# This is used when tRNA needs to be moved to pseudogene
# and is way faster than aceperl
##########################################################

my $command=<<EOF;
  find Transcript *.t* where Species = "C*elegans"; 
  show -a -T -f /tmp/trnas.ace
EOF

open (DUMP, "| $tace $database") || die "Failed to connect to $database";
print DUMP $command;

my @pseudo = `less /tmp/trnas.ace`;
my (%tRNA_ace, $tRNA);

foreach (@pseudo) {
  chomp;
  if ($_ =~ /^Transcript : \"(.+)\" .+/){
    $tRNA = $1;
  }
  else {
    push(@{$tRNA_ace{$tRNA}}, $_);
  }
}

#################################################################
# comparing tRNA coordinates of existing ones with predicted ones
# write acefile and output a stats list
#################################################################

my $HX_a = 0; my $HX_c = 0; my $HX_e = 0; my $HX_g = 0; my $HX_i = 0;
my $RW_b = 0; my $RW_d = 0; my $RW_f = 0; my $RW_h = 0; my $RW_j = 0;
my $HX_k = 0; my $HX_m = 0; my $HX_o = 0; my $HX_q = 0; my $HX_s = 0;
my $RW_l = 0; my $RW_n = 0; my $RW_p = 0; my $RW_r = 0; my $RW_t = 0;

my $same_start_same_end = 0;
my $same_start_same_end_same_type = 0;
my $same_start_end_pseudo = 0;
my $same_start_end_not_pseudo = 0;

my $same_start_diff_end = 0;
my $same_start_diff_end_same_type = 0;
my $same_start_diff_end_pseudo = 0;
my $same_start_diff_end_not_pseudo = 0;

my $diff_start_diff_end = 0;
my $diff_start_diff_end_pseudo = 0;
my $diff_start_diff_end_not_pseudo = 0;

my $diff_start_same_end = 0;
my $diff_start_same_end_same_type = 0;
my $diff_start_same_end_pseudo = 0;
my $diff_start_same_end_not_pseudo = 0;

my ($codon, $g_chrom, $g_end, $p_chrom, $p_end, @clone_coords); 
my %clone_name =();

####################################################################################################
# loop through each predicted tRNA to compare their start and end coords with those of existing ones
####################################################################################################

open(ACE, ">tRNA_prediction.ace") || die $!;

foreach $start (sort keys %predict_cols){

  # ${@{$gff_cols{$start}}}[2]:     chrom from %gff_col
  # ${@{$predict_cols{$start}}}[3]: chrom from %predict_col
  # for clarity  
  $p_end       = ${@{$predict_cols{$start}}}[0]; # end coord of predicted tRNA
  $anticodon   = ${@{$predict_cols{$start}}}[2]; $anticodon =~ tr/T/U/; 
  $type        = ${@{$predict_cols{$start}}}[1];

  # predicted tRNA with same start coord as in existing one on the same chromosome
  if (exists $gff_cols{$start} && exists ${@{$gff_cols{$start}}}[2] && exists ${@{$predict_cols{$start}}}[3] 
    && ${@{$gff_cols{$start}}}[2] eq ${@{$predict_cols{$start}}}[3] ){
    
    $g_end        = ${@{$gff_cols{$start}}}[0];       # end coord of existing tRNA
    $transcript   = ${@{$gff_cols{$start}}}[1];

    # identical start and end coords
    if ($g_end == $p_end){
      $same_start_same_end++; 

      if ($type eq ${@{$tRNA_type{$transcript}}}[0]){
	$same_start_same_end_same_type++; 
       
	print ACE "\n\/\/(A1:) Same start, same end, type identical: ${@{$tRNA_type{$transcript}}}[0]<->$type\n";
	print ACE "\nTranscript : \"$transcript\"\n";
	print ACE "-D Method\n";
	print ACE "-D Property\n";
	print ACE "\nTranscript : \"$transcript\"\n";
	print ACE "Type\t\"$type\"\n";
	print ACE "Anticodon\t\"$anticodon\"\n";
	print ACE "Method \"tRNAscan-SE-1.23\"\n";
      }
      if ($type ne ${@{$tRNA_type{$transcript}}}[0] && $type eq "Pseudo"){
	$same_start_end_pseudo++;

	$HX_a++ if ${@{$tRNA_type{$transcript}}}[1] eq "HX";
        $RW_b++ if ${@{$tRNA_type{$transcript}}}[1] eq "RW";

	print ACE "\n\/\/(A2:) same S, same E, $transcript is ${@{$tRNA_type{$transcript}}}[0], prediction is $type\n";
        print ACE "\/\/Move $transcript to Pseudogene Class\n\n";
        move_to_pseudo($start, $transcript, "pseudo");
      } 
      if ($type ne ${@{$tRNA_type{$transcript}}}[0] && $type ne "Pseudo"){
        $same_start_end_not_pseudo++;

	$HX_c++ if ${@{$tRNA_type{$transcript}}}[1] eq "HX";
        $RW_d++ if ${@{$tRNA_type{$transcript}}}[1] eq "RW";

	print ACE "\n\/\/(A3:) same S, same E, $transcript is ${@{$tRNA_type{$transcript}}}[0], prediction is $type\n";
	write_ace($start, $transcript, $type, $anticodon, "exon");
      }
    }
    # same start coord, diff end coord
    if ($g_end != $p_end){
      $same_start_diff_end++;

      # same start coord, diff end coord, same type -> change end coord
      if ($type eq ${@{$tRNA_type{$transcript}}}[0]){
	$same_start_diff_end_same_type++; 
	
	$HX_e++ if ${@{$tRNA_type{$transcript}}}[1] eq "HX";
        $RW_f++ if ${@{$tRNA_type{$transcript}}}[1] eq "RW";

	print ACE "\n\/\/(B1:) Same start, diff end, type identical: ${@{$tRNA_type{$transcript}}}[0]<->$type\n";
        print ACE "\/\/start: $start, end: $p_end\n";
        write_ace($start, $transcript, $type, $anticodon);
      }
      if ($type ne ${@{$tRNA_type{$transcript}}}[0] && $type eq "Pseudo"){
        $same_start_diff_end_pseudo++;
	
	$HX_g++ if ${@{$tRNA_type{$transcript}}}[1] eq "HX";
        $RW_h++ if ${@{$tRNA_type{$transcript}}}[1] eq "RW";

        print ACE "\n\/\/(B2:) Same start, diff end, type diff: ${@{$tRNA_type{$transcript}}}[0], prediction is $type\n";
        print ACE "\/\/start: $start, end: $p_end\n";
        move_to_pseudo($start, $transcript, "pseudo");
      }
      if ($type ne ${@{$tRNA_type{$transcript}}}[0] && $type ne "Pseudo"){
	print ACE "\n\/\/(B3:) $transcript is $tRNA_type{$transcript}, prediction is $type\n";
        $same_start_diff_end_not_pseudo++;
	
	$HX_i++ if ${@{$tRNA_type{$transcript}}}[1] eq "HX";
        $RW_j++ if ${@{$tRNA_type{$transcript}}}[1] eq "RW";

        write_ace($start, $transcript, $type, $anticodon);
      }
    }
  }

  ## predicted tRNAs with diff start coord
  if (!exists $gff_cols{$start}) {
    
    # ${@{$gff_cols_end{$p_end}}}[2]: chrom from %gff_col
    # ${@{$predict_cols{$start}}}[3]: chrom from %predict_col
    
    # get transcript with identical end coord
    $transcript = ${@{$gff_cols_end{$p_end}}}[1] if exists $gff_cols_end{$p_end};
    $transcript = "NA" if !exists $gff_cols_end{$p_end};

    # diff start coord AND end coord 
    if (!exists $gff_cols_end{$p_end}){
     $diff_start_diff_end++;
    } 

    # diff start coord AND end coord and pred. type is pseudo
    if (!exists $gff_cols_end{$p_end} && $type eq "Pseudo"){
      print ACE "\n\/\/(C1:)Diff start & end coord: $start : ${@{$predict_cols{$start}}}[0], pred. type: $type\n";
      $diff_start_diff_end_pseudo++;
      my $LAB = move_to_pseudo($start, ${@{$predict_cols{$start}}}[0], "new");  # start and end coords of new tRNA
      $HX_k++ if $LAB eq "HX";
      $RW_l++ if $LAB eq "RW";
    }
    
    # diff start coord AND end coord and pred. type is NOT pseudo
    if (!exists $gff_cols_end{$p_end} && $type ne "Pseudo"){
      print ACE "\n\/\/(C2:)Diff start & end coord: $start : $p_end, pred. type $type\n";  

      $diff_start_diff_end_not_pseudo++;
      my $LAB = write_ace($start, $p_end, $type, $anticodon);
      $HX_m++ if $LAB eq "HX";
      $RW_n++ if $LAB eq "RW";
    }

    # diff start coord, but same end coord on the same chromosome
    if (exists $gff_cols_end{$p_end} && exists  ${@{$gff_cols_end{$p_end}}}[2] && exists ${@{$predict_cols{$start}}}[3]
        && ${@{$gff_cols_end{$p_end}}}[2] eq ${@{$predict_cols{$start}}}[3]){
        $diff_start_same_end++;
    }

    # diff start coord, but same end coord on the same chromosome, same type 
    if (exists $gff_cols_end{$p_end} && exists  ${@{$gff_cols_end{$p_end}}}[2] && exists ${@{$predict_cols{$start}}}[3]
        && ${@{$gff_cols_end{$p_end}}}[2] eq ${@{$predict_cols{$start}}}[3]
        && $type eq ${@{$tRNA_type{$transcript}}}[0]){
      $diff_start_same_end_same_type++;

      $HX_o++ if ${@{$tRNA_type{$transcript}}}[1] eq "HX";
      $RW_p++ if ${@{$tRNA_type{$transcript}}}[1] eq "RW";   

      print ACE "\n\/\/(C3:)Diff start & same end coord, same type: $start -> $p_end, $type<->${@{$tRNA_type{$transcript}}}[0]\n"; 
      write_ace($start, $transcript, $type, $anticodon);
    }

    # diff start coord, but same end coord on the same chromosome, diff type and pred. is pseudo
    if (exists $gff_cols_end{$p_end} && exists  ${@{$gff_cols_end{$p_end}}}[2] && exists ${@{$predict_cols{$start}}}[3]
        && ${@{$gff_cols_end{$p_end}}}[2] eq ${@{$predict_cols{$start}}}[3] 
        && $type ne ${@{$tRNA_type{$transcript}}}[0] && $type eq "Pseudo" ){ 

      $diff_start_same_end_pseudo++;
      $HX_q++ if ${@{$tRNA_type{$transcript}}}[1] eq "HX";
      $RW_r++ if ${@{$tRNA_type{$transcript}}}[1] eq "RW";

      print ACE "\n\/\/(C4:)Diff start & same end coord, diff type: $start -> ${@{$gff_cols_end{$p_end}}}[0], SAME end coord: $p_end, ";
      print ACE "\n\/\/pred. type: $type, current: ${@{$tRNA_type{$transcript}}}[0]\n"; 
      print ACE "\n\/\/(C4a:) $transcript: update end codon from ${@{$gff_cols_end{$p_end}}}[0] -> $start (", $start-${@{$gff_cols_end{$p_end}}}[0], ")\n"; 
      move_to_pseudo($start, $transcript, "pseudo"); 
    }

    # diff start coord, but same end coord on the same chromosome, diff type and pred. is NOT pseudo
    if (exists $gff_cols_end{$p_end} && exists  ${@{$gff_cols_end{$p_end}}}[2] && exists ${@{$predict_cols{$start}}}[3]
        && ${@{$gff_cols_end{$p_end}}}[2] eq ${@{$predict_cols{$start}}}[3] 
        && $type ne ${@{$tRNA_type{$transcript}}}[0] && $type ne "Pseudo" ){

      $diff_start_same_end_not_pseudo++;
      $HX_s++ if ${@{$tRNA_type{$transcript}}}[1] eq "HX"; 
      $RW_t++ if ${@{$tRNA_type{$transcript}}}[1] eq "RW";

      print ACE "\/\/(C5:)Diff start & same end coord: $start -> ${@{$gff_cols_end{$p_end}}}[0], SAME end coord: $p_end, ";
      print ACE "\/\/pred. type: $type, current: ${@{$tRNA_type{$transcript}}}[0]\n";
      print ACE "\/\/(C5a:) $transcript: update end codon from ${@{$gff_cols_end{$p_end}}}[0] -> $start (", $start-${@{$gff_cols_end{$p_end}}}[0], ")\n";
      write_ace($start, $transcript, $type, $anticodon);
    }
  }
}

$curr_db->close;
system("rm -f /tmp/trna.ace");

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

print "\nComparing existing tRNA (", scalar (keys %gff_cols), ") with predicted tRNA (", scalar (keys %predict_cols), ") based on $current\n\n";
print "Legends:\n--------\n";
print "<Start>: start coord, <End>: end coord, <Type>: tRNA type, <#>: number\n";
print "<+>: same, <->: diff, <P>: pred. type is pseudo, <N>: pred. type is not pseudo\n\n";

print "Start\tEnd\tType\t  #\tHX\tRW\t\n";
print "-------------------------------------------\n";
print "  +  \t + \t   \t $same_start_same_end\n";
print "  +  \t + \t + \t $same_start_same_end_same_type\n";
print "  +  \t + \t P \t $same_start_end_pseudo\t$HX_a\t$RW_b\n";
print "  +  \t + \t N \t $same_start_end_not_pseudo\t$HX_c\t$RW_d\n";
print "\n";
print "  +  \t - \t   \t $same_start_diff_end\n";
print "  +  \t - \t + \t $same_start_diff_end_same_type\t$HX_e\t$RW_f\n";
print "  +  \t - \t P \t $same_start_diff_end_pseudo\t$HX_g\t$RW_h\n";
print "  +  \t - \t N \t $same_start_diff_end_not_pseudo\t$HX_i\t$RW_j\n";
print "\n";
print "  -  \t + \t   \t $diff_start_same_end\n";
print "  -  \t + \t + \t $diff_start_same_end_same_type\t$HX_o\t$RW_p\n";
print "  -  \t + \t P \t $diff_start_same_end_pseudo\t$HX_q\t$RW_r\n";
print "  -  \t + \t N \t $diff_start_same_end_not_pseudo\t$HX_s\t$RW_t\n";
print "\n";
print "  -  \t - \t   \t $diff_start_diff_end (NEW)\n";
print "  -  \t - \t P \t $diff_start_diff_end_pseudo\t$HX_k\t$RW_l\n";
print "  -  \t - \t N \t $diff_start_diff_end_not_pseudo\t$HX_m\t$RW_n\n";

#######################
# s u b r o u t i n e s
#######################

sub write_ace {
  my ($coord, $transcript, $type, $anticodon, $exon) = @_;
  
  my $chrom = "CHROMOSOME_".${@{$predict_cols{$coord}}}[3];
  my ($version, @seq_names);
  my $HX = (); my $RW =();  

  if ($transcript !~ /^[0-9]{1,}$/){
    print ACE "\nTranscript : \"$transcript\"\n";
    print ACE "-D Source_Exons\n" if !$exon;
    print ACE "-D Brief_identification \"tRNA-$type\"\n";
    print ACE "-D DB_remark\n";
    print ACE "-D Property\n";
    print ACE "-D Method\n";
   
    # watch out for multiple exons: the transcript obj should use the same clone coords
    my $start_coord = $coord;
    my $end_coord;
    if (exists $predict_cols_span{$coord}){
      $start_coord = ${@{$predict_cols_span{$start_coord}}}[0];
      $end_coord   = ${@{$predict_cols_span{$start_coord}}}[1];
    } 
    else {  
      $end_coord = ${@{$predict_cols{$coord}}}[0];
    }
    
    my $coords  = Coords_converter->invoke($database);
    @clone_coords = $coords->LocateSpan($chrom, $start_coord, $end_coord);
   
    print ACE "\nSequence : \"$clone_coords[0]\"\n";
    print ACE "Transcript_child $transcript $clone_coords[1] $clone_coords[2]\n";
  }
  
  if ($transcript =~ /^[0-9]{1,}$/){
    print ACE "\/\/$transcript  B4\n";

    # watch out for multiple exons: the transcript obj should use the same clone coords
    my $start_coord = $coord;
    my $end_coord;
    if (exists $predict_cols_span{$coord}){
      $start_coord = ${@{$predict_cols_span{$start_coord}}}[0];
      $end_coord   = ${@{$predict_cols_span{$start_coord}}}[1];     
    } 
    else {
      $end_coord = $transcript;
    }

    my $coords  = Coords_converter->invoke($database);
    @clone_coords = $coords->LocateSpan($chrom, $start_coord, $end_coord);
    print ACE "\/\/@clone_coords\n";
    @seq_names = ();  
 
    $HX = "HX" if $clone_lab{$clone_coords[0]} eq "HX";
    $RW = "RW" if $clone_lab{$clone_coords[0]} eq "RW";

    # assign seq. name to tRNA/Pseudogene
    if (exists $clone_last{$clone_coords[0]}){
      my $last = $clone_last{$clone_coords[0]};
      $last = $last+1 if !(exists $predict_cols_exon{$coord} && ${@{$predict_cols_exon{$coord}}}[0] != 1);
      $transcript = "$clone_coords[0]\.t".$last;
      $clone_last{$clone_coords[0]} = $last;
    }
    else {
      $transcript = "$clone_coords[0]\.t1";
      $clone_last{$clone_coords[0]} = 1;
    }
    print ACE "\nSequence : \"$clone_coords[0]\"\n";
    print ACE "Transcript_child $transcript $clone_coords[1] $clone_coords[2]\n";
    print ACE "\nSequence : \"$transcript\"\n";
    print ACE "Corresponding_transcript $transcript\n";
  }
  print ACE "\nTranscript : \"$transcript\"\n";
  print ACE "Sequence $clone_coords[0]\n" if $clone_coords[0];
  print ACE "From_Laboratory \"$clone_lab{$clone_coords[0]}\"\n" if $clone_coords[0];
  print ACE "DB_remark \"Predicted tRNA\"\n";
  print ACE "Species \"Caenorhabditis elegans\"\n";
  print ACE "Brief_identification \"tRNA-$type\"\n";
  print ACE "Type\t\"$type\"\n";
  print ACE "Anticodon\t\"$anticodon\"\n";
  print ACE "Method \"tRNAscan-SE-1.23\"\n"; 

  if (exists $predict_cols_exon{$coord} ){
    print ACE "t1-Source_Exons\t${@{$predict_cols_exon{$coord}}}[0]\t${@{$predict_cols_exon{$coord}}}[1]\n";          # if has > 1 exons
  }
  if ($coord < ${@{$predict_cols{$coord}}}[0]){
    print ACE "t2-Source_Exons 1 ", ${@{$predict_cols{$coord}}}[0]-$coord+1, "\n" if ${@{$predict_cols{$coord}}}[5] == 0;   # if only 1 exon
  } 
  else {
    print ACE "t3-Source_Exons 1 ", $coord-${@{$predict_cols{$coord}}}[0]+1, "\n" if ${@{$predict_cols{$coord}}}[5] == 0;   # if only $
  }
  return $HX if $HX;
  return $RW if $RW;
} 

sub move_to_pseudo {
  my ($coord, $transcript, $pseudo) = @_;
  my $HX =(); my $RW =();
  my $chrom = "CHROMOSOME_".${@{$predict_cols{$coord}}}[3];
  my $parent;

  print ACE "\n-D Transcript \"$transcript\"\n" if $pseudo && $pseudo ne "new";
  if ($transcript !~ /^[0-9]+$/){
    print ACE "\nPseudogene \"$transcript\"\n"; 
    print ACE "-D DB_info\n";

    # watch out for multiple exons: the transcript obj should use the same clone coords
    my $start_coord = $coord;
    my $end_coord;
    if (exists $predict_cols_span{$start_coord}){
      $start_coord = ${@{$predict_cols_span{$start_coord}}}[0];
      $end_coord   = ${@{$predict_cols_span{$start_coord}}}[1];
    }
    else {
      $end_coord = ${@{$predict_cols{$coord}}}[0];
    }

    my $coords  = Coords_converter->invoke($database);
    @clone_coords = $coords->LocateSpan($chrom, $start_coord, $end_coord);
 
    print ACE "\nSequence \"$clone_coords[0]\"\n";
    print ACE "Pseudogene \"$transcript\" $clone_coords[1] $clone_coords[2]\n";

    print ACE "\nPseudogene : \"$transcript\"\n";
  }
  if ($transcript =~ /^[0-9]+$/){
    print ACE "\/\/$transcript   Before\n";

    # watch out for multiple exons: the transcript obj should use the same clone coords
    my $start_coord = $coord;
    my $end_coord;
    if (exists $predict_cols_span{$start_coord}){
      $start_coord = ${@{$predict_cols_span{$start_coord}}}[0];
      $end_coord   = ${@{$predict_cols_span{$start_coord}}}[1];
    }
    else {
      $end_coord = $transcript;
    }
     
    my $coords  = Coords_converter->invoke($database);
    @clone_coords = $coords->LocateSpan($chrom, $start_coord, $end_coord); 
    print ACE "\/\/@clone_coords\n";
    
    print  "$coord : $clone_coords[0] : $clone_lab{$clone_coords[0]}\n";   
    $HX = "HX" if $clone_lab{$clone_coords[0]} eq "HX";
    $HX = "HX" if $clone_coords[0] =~ /SUPERLINK_HX.+/;
    
    $RW = "RW" if $clone_lab{$clone_coords[0]} eq "RW";      
    $RW = "RW" if $clone_coords[0] =~ /SUPERLINK_RW.+/; 
      
    # assign seq. name to tRNA/Pseudogene
    if (exists $clone_last{$clone_coords[0]}){
      my $last = $clone_last{$clone_coords[0]};
      $last = $last+1 if !(exists $predict_cols_exon{$coord} && ${@{$predict_cols_exon{$coord}}}[0] != 1);
      $transcript = "$clone_coords[0]\.t".$last;    
      $clone_last{$clone_coords[0]} = $last; 
    }
    else { 
      $transcript = "$clone_coords[0]\.t1";
      $clone_last{$clone_coords[0]} = 1;
    }
    print ACE "\nSequence \"$clone_coords[0]\"\n";
    print ACE "Pseudogene \"$transcript\" $clone_coords[1] $clone_coords[2]\n";
    print ACE "\nPseudogene : \"$transcript\"\n"; 
    print ACE "Sequence \"$clone_coords[0]\"\n";
    print ACE "From_Laboratory \"HX\"\n" if $HX;
    print ACE "From_Laboratory \"RW\"\n" if $RW;	
  }
  print ACE "Type \"RNA_pseudogene\"\n";
  print ACE "Remark \"Predicted tRNA by tRNAscan-SE-1.11 but predicted as pseudogene by tRNAscan-SE-1.23\"\n" if $pseudo ne "new";
  print ACE "Remark \"Predicted pseudogene by tRNASCAN-SE-1.23\"\n" if $pseudo eq "new";
  print ACE "Species \"Caenorhabditis elegans\"\n";
  print ACE "Method \"Pseudogene\"\n";

  if (exists $predict_cols_exon{$coord} ){
    print ACE "t4-Source_Exons\t${@{$predict_cols_exon{$coord}}}[0]\t${@{$predict_cols_exon{$coord}}}[1]\n";               # if has > 1 exons
  }
  # + strain
  if ($coord < ${@{$predict_cols{$coord}}}[0]){
    print ACE "t5-Source_Exons 1 ", ${@{$predict_cols{$coord}}}[0]-$coord+1, "\n" if ${@{$predict_cols{$coord}}}[5] == 0;   # if only 1 exon
  }
  # - strain
  else {
    print ACE "t6-Source_Exons 1 ", $coord-${@{$predict_cols{$coord}}}[0]+1, "\n" if ${@{$predict_cols{$coord}}}[5] == 0;   # if only 1 exon
  }
  if (exists $tRNA_ace{$transcript}){
    foreach (@{$tRNA_ace{$transcript}}) {
      chomp;
      if ($_ =~ /Remark/ || $_ =~ /^Method/ || $_ =~ /^Structure/ ||
          $_ =~ /Brief_identification/ || $_ =~ /^Interpolated_map_position/ ||
          $_ =~ /^Properties/ ) {}
      else{
        print ACE "$_\n";
      }
    }
  }
  return $HX if $HX && $pseudo eq "new";
  return $RW if $RW && $pseudo eq "new";
}

__END__


=head2 NAME - parse_predicted_tRNA.pl  

=head2 DESCRIPTION

            This script parses the output of tRNA prediction suite tRNASCAN-SE 1.23 (current) and write ace file 
            for predicted tRNA or existing tRNA objects.
         
=head3 <USAGE> 
 

=head2 Options: [h or help] [f or file]


B<-help:>     
            Displays documentation

B<-file:>   
            specify file generated by tRNASCAN-SE


