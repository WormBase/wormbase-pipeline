#!/usr/local/bin/perl5.6.1 -w

# get_interpolated_gmap.pl

# by Chao-Kung Chen [030113]

# This script calculates interpolated gmap for CDS/transcripts lying between as well as outside genetic markers.
# Output ace file of such information

# Last updated on: $Date: 2003-03-31 10:37:47 $
# Last updated by: $Author: ck1 $

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Cwd 'chdir';
use Getopt::Long;

my $start=`date +%H:%M:%S`; chomp $start;
my $rundate = `date +%y%m%d`; chomp $rundate;

###################################################
# variables and command-line options with aliases #
###################################################

my ($diff, $reverse);

GetOptions ("diff" => \$diff,
            "rev|reverse" => \$reverse,
           );

my $curr_db_dir = "/wormsrv2/current_DB/";
#my $autoace_dir = "/wormsrv2/autoace/";

#################################################
# look for latest WS release in GFF_SPLITS folder
#################################################

my $gff_dir = "/wormsrv2/autoace/GFF_SPLITS/";

my $ga_dir="/wormsrv1/geneace";
#my $ga_dir="/wormsrv1/chaokung/wormdata/CK1_GENEACE/";

my @versions=dataset($gff_dir, "folder");

my @order = sort {$a <=> $b} @versions;

my $gff_location = "/wormsrv2/autoace/GFF_SPLITS/WS"."$order[-1]";

#my $cmp_file = glob "/wormsrv2/logs/cmp_gmap_with_coord_order.".$rundate;
my $cmp_file = glob "/wormsrv1/chaokung/DOCS/cmp_gmap_with_coord_order.".$rundate;
if($cmp_file){system("rm -f $cmp_file")}

open(OUT, ">>$cmp_file") || die $!;
system ("chmod 755 $cmp_file");
print OUT "\nGFF file version from $gff_location\n\n";
             
chdir $gff_location;

my (@data, %CDS_mapping, %CDS_isoforms_mapping, %CDS_variants, %predicted_gene_to_locus);
my (%chrom_meanCoord, %genetics_mapping, %I_gmap_cds_locus, %II_gmap_cds_locus, %III_gmap_cds_locus,
    %IV_gmap_cds_locus, %V_gmap_cds_locus, %X_gmap_cds_locus);
my ($cds, $parts, @coords, $i, $mean_coords, %cds_mean, %mean_coord_cds,
    %mean_coord_cds_I,  %mean_coord_cds_II,  %mean_coord_cds_III,
    %mean_coord_cds_IV,  %mean_coord_cds_V,  %mean_coord_cds_X,
   );

if ($reverse){
  print "\nChecking reverse physicals for genetic markers . . . .\n";
  #my $logfile="/wormsrv2/logs/reverse_physicals.".$rundate;
  my $logfile="/wormsrv1/chaokung/DOCS/reverse_physicals.".$rundate;
  open (LOG, ">$logfile") || die $!;
  system("chmod 755 $logfile");
  print LOG "Checking linearity of gmap markers (reverse physicals) ....\n";
  print LOG "\n";
}

if ($diff){print "Checking for chrom. mapping discrepancies by genetics / coordinates......\n"}

if(!$diff){

  #my $acefile="/wormsrv2/logs/interpolated_gmap.ace.".$rundate;
  my $acefile="/wormsrv1/chaokung/DOCS/interpolated_gmap.ace.".$rundate;
  open (ACE, ">$acefile") || die "Can't output file!\n";
  system("chmod 755 $acefile");

  #############################################
  # get list of predicted CDSes linked to locus
  #############################################

  my $predicted_CDS_linked_to_locus=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/get_predicted_gene_to_locus.def" quit
EOF
    
  open (FH, "echo '$predicted_CDS_linked_to_locus' | tace $curr_db_dir | ") || die "Couldn't access current_DB\n";
  #open (FH, "echo '$predicted_CDS_linked_to_locus' | tace $autoace_dir | ") || die "Couldn't access autoace\n";
  while (<FH>){
    chomp($_);
    if ($_ =~ /^\".+\"$/){
      my ($predict, $locus)= split(/\s+/, $_);
      $predict =~ s/\"//g;
      $locus =~ s/\"//g;
      $predicted_gene_to_locus{$predict} = $locus;
    }
  }
  close FH;
}

##########################################################
# retrieve data from gff file: CHROMOSOME_number.genes.gff
##########################################################

my @gff_files_cds=dataset($gff_location, "genes");
my @gff_files_rna=dataset($gff_location, "rna");
  
my $cds_count=0;
my $rna_count=0;
  
foreach (@gff_files_cds){
  @data = `less $_ | cut -f 1,4,5,9`;
  foreach (@data){
    if ($_ =~ /^CHROM.+/){
      my ($chrom, $left, $right, $junk1, $CDS)= split(/\s+/,$_);
      $CDS =~ s/\"//g;
      my $ori = $CDS;
      $cds_count++;
      if ($CDS =~ /(.+\.\d+)\D+/){
        $CDS = $1;
        push(@{$CDS_variants{$CDS}}, $ori); #parent of isoforms
        push(@{$CDS_isoforms_mapping{$CDS}}, $chrom, $left, $right);
      }
      else {
        $chrom =~ s/CHROMOSOME_//;
        push(@{$CDS_mapping{$CDS}}, $chrom, $left, $right);
      }
    }
  }
}

print "all cds(with isoforms): ",$cds_count, "\n";

foreach (@gff_files_rna){
  @data = `less $_ | cut -f 1,4,5,9,10`;
  foreach (@data){
    if ($_ =~ /^CHR.+/){
      my ($chrom, $left, $right, $junk, $RNA)= split(/\s/,$_);
      $RNA =~ s/\"//g;
      my $ori = $RNA;
      $rna_count++;
      $chrom =~ s/CHROMOSOME_//;
      push(@{$CDS_mapping{$RNA}}, $chrom, $left, $right);
    }
  }
}
        
print "all transcripts: ",$rna_count, "\n";
        
##############################################################
# from GFF files: genes.gff / rna.gff
# get mean value of chrom coords of each CDS/Transcript,
# mean of isoforms is the mean of the lowest and highest value
##############################################################
     
my $count=0;

foreach (sort keys %CDS_mapping){

  ####################
  # cds has no isoform   
  ####################
  
  # print "$_ -> ${@{$CDS_mapping{$_}}}[1] --- ${@{$CDS_mapping{$_}}}[2]\n";
  $mean_coords = (${@{$CDS_mapping{$_}}}[1] + ${@{$CDS_mapping{$_}}}[2]) / 2;
      
  $cds_mean{$_}= $mean_coords; # key: cds(w/0 isoform) value: mean coords
      
  if (${@{$CDS_mapping{$_}}}[0] eq "I")  {$mean_coord_cds_I  {$mean_coords}=$_}
  if (${@{$CDS_mapping{$_}}}[0] eq "II") {$mean_coord_cds_II {$mean_coords}=$_}
  if (${@{$CDS_mapping{$_}}}[0] eq "III"){$mean_coord_cds_III{$mean_coords}=$_}
  if (${@{$CDS_mapping{$_}}}[0] eq "IV") {$mean_coord_cds_IV {$mean_coords}=$_}
  if (${@{$CDS_mapping{$_}}}[0] eq "V")  {$mean_coord_cds_V  {$mean_coords}=$_}
  if (${@{$CDS_mapping{$_}}}[0] eq "X")  {$mean_coord_cds_X  {$mean_coords}=$_}
}

foreach (sort keys %CDS_isoforms_mapping){
  $parts = scalar @{$CDS_isoforms_mapping{$_}};

  ##################
  # cds has isoforms
  ##################

  for ($i=0; $i < $parts; $i=$i+3 ){
    push(@coords, ${@{$CDS_isoforms_mapping{$_}}}[$i+1], ${@{$CDS_isoforms_mapping{$_}}}[$i+2]);  # get coords of all isoforms
  }
  # print "$_ -> @coords\n";
  @coords = sort {$a <=> $b} @coords;
  $mean_coords = ($coords[-1] + $coords[0]) / 2;

  $cds_mean{$_}= $mean_coords;

  if (${@{$CDS_isoforms_mapping{$_}}}[0] eq "I")  {$mean_coord_cds_I  {$mean_coords}=$_}
  if (${@{$CDS_isoforms_mapping{$_}}}[0] eq "II") {$mean_coord_cds_II {$mean_coords}=$_}
  if (${@{$CDS_isoforms_mapping{$_}}}[0] eq "III"){$mean_coord_cds_III{$mean_coords}=$_}
  if (${@{$CDS_isoforms_mapping{$_}}}[0] eq "IV") {$mean_coord_cds_IV {$mean_coords}=$_}
  if (${@{$CDS_isoforms_mapping{$_}}}[0] eq "V")  {$mean_coord_cds_V  {$mean_coords}=$_}
  if (${@{$CDS_isoforms_mapping{$_}}}[0] eq "X")  {$mean_coord_cds_X  {$mean_coords}=$_}
  @coords =();
}

################################################################################################
# retrieve locus data: chromosome, map position, genomic_sequence/transcript from gmap (Geneace)
################################################################################################
 
my $marker_gmap_of_each_chrom=<<EOF;
Table-maker -p "/wormsrv1/geneace/wquery/marker_gmap_of_each_chrom.def" quit
EOF
  
my ($mean_coord, %chrom_pos);
 
open (FH, "echo '$marker_gmap_of_each_chrom' | tace $ga_dir | ") || die "Couldn't access Geneace\n";
while (<FH>){
  chomp($_);
  # $4=cds/transcript, $2=chrom, $1=locus, $3=gmap position
  if ($_ =~ /^\"(.+)\"\s+\"(.+)\"\s+(-\d+\.\d+|\d+\.\d+)\s+\"(.+)\"/){   # postive and negative values
    my $locus = $1;
    my $chrom = $2;
    my $gmap = $3;
    my $cds = $4;
    #print "$1 -> $2 -> $3 -> $4\n";
    if ($2 eq "I"){push (@{$I_gmap_cds_locus{$3}}, $4, $1)}
    if ($2 eq "II"){push (@{$II_gmap_cds_locus{$3}}, $4, $1)}
    if ($2 eq "III"){push (@{$III_gmap_cds_locus{$3}}, $4, $1)}
    if ($2 eq "IV"){push (@{$IV_gmap_cds_locus{$3}}, $4, $1)}
    if ($2 eq "V"){push (@{$V_gmap_cds_locus{$3}}, $4, $1)}
    if ($2 eq "X"){push (@{$X_gmap_cds_locus{$3}}, $4, $1)}
    if ($cds =~ /(.+\.\d+)\D+/){
      my $cds = $1;
      $mean_coord=$cds_mean{$cds};
      push(@{$genetics_mapping{$cds}}, $chrom, $locus, $gmap, $mean_coord);
      push(@{$chrom_pos{$chrom}}, $gmap, $mean_coord, $locus, $cds);
    }
    if ($cds !~ /(.+\.\d+)\D+/) { 
      $mean_coord=$cds_mean{$cds};
      push(@{$genetics_mapping{$cds}}, $chrom, $locus, $gmap, $mean_coord);
      push(@{$chrom_pos{$chrom}}, $gmap, $mean_coord, $locus, $cds);
    }
  }
}
close FH;
  
#######################################################################
# output discrepancy of chromosomal mapping by genetics and coordinates 
# this is called from geneace_check.pl
#######################################################################

$count = 0;
if ($diff){
  #open(LOG, ">/wormsrv2/logs/mapping_diff.".$rundate) || die $!;
  open(LOG, ">/wormsrv1/chaokung/DOCS/mapping_diff.".$rundate) || die $!;
  foreach (sort keys %genetics_mapping){
    
    if (exists ${@{$genetics_mapping{$_}}}[0] && exists ${@{$CDS_isoforms_mapping{$_}}}[0] ){
      ${@{$CDS_isoforms_mapping{$_}}}[0] =~ s/CHROMOSOME_//;
      if (${@{$genetics_mapping{$_}}}[0] ne ${@{$CDS_isoforms_mapping{$_}}}[0]){
        $count++;
        print  LOG "ERROR: $_ (${@{$genetics_mapping{$_}}}[1]) is mapped to ${@{$genetics_mapping{$_}}}[0] by genetics\n";
        print  LOG "ERROR: $_ (${@{$genetics_mapping{$_}}}[1]) is mapped to ${@{$CDS_isoforms_mapping{$_}}}[0] by coordinates.\n\n";
      }
    }
    if (exists ${@{$genetics_mapping{$_}}}[0] && exists ${@{$CDS_mapping{$_}}}[0] ){
      ${@{$CDS_mapping{$_}}}[0] =~ s/CHROMOSOME_//;
      if (${@{$genetics_mapping{$_}}}[0] ne ${@{$CDS_mapping{$_}}}[0]){
        $count++;
        print  LOG "ERROR: $_ (${@{$genetics_mapping{$_}}}[1]) is mapped to ${@{$genetics_mapping{$_}}}[0] by genetics\n";
        print  LOG "ERROR: $_ (${@{$genetics_mapping{$_}}}[1]) is mapped to ${@{$CDS_mapping{$_}}}[0] by coordinates.\n\n";
      }
    }
  }
  print LOG "There are $count discrepancies in genetics\/chrom. coords mapping\n";
  exit(0);
}

 ##################################################################
  # hard coded chrom length - as change in the future is very seldom
  ##################################################################
    
  my %chrom_length = (
                       I   => 15080471,
                       II  => 15279300,
                       III => 13783267,
                       IV  => 17493790,
                       V   => 20922239,
                       X   => 17705013,
                     );

#################################################################
# get end coords of each chromosom -> 
# for calculating int. gmap of CDSes lying outside of tip markers
#################################################################

my (@all_coords, %All_Coords);

@all_coords= keys %mean_coord_cds_I;
my @all_coords_I = sort {$a <=>$b} @all_coords;

@all_coords= keys %mean_coord_cds_II;
my @all_coords_II = sort {$a <=>$b} @all_coords;

@all_coords= keys %mean_coord_cds_III;
my @all_coords_III = sort {$a <=>$b} @all_coords;

@all_coords= keys %mean_coord_cds_IV;
my @all_coords_IV = sort {$a <=>$b} @all_coords;

@all_coords= keys %mean_coord_cds_V;
my @all_coords_V = sort {$a <=>$b} @all_coords;

@all_coords= keys %mean_coord_cds_X;
my @all_coords_X = sort {$a <=>$b} @all_coords;

%All_Coords=(
              I   => [$all_coords_I[0], $all_coords_I[-1]],
              II  => [$all_coords_II[0], $all_coords_II[-1]],
              III => [$all_coords_III[0], $all_coords_III[-1]],
              IV  => [$all_coords_IV[0], $all_coords_IV[-1]],
              V   => [$all_coords_V[0], $all_coords_V[-1]],
              X   => [$all_coords_X[0], $all_coords_X[-1]],
            );

#########################################################
#  
# 
#########################################################

my (@chroms, $ea, @pos_order, $part, @unique_pos_order, %pos_order_to_mean_coord_locus_cds, 
    @mean_coords_order, $chrom, @chrom_I, @chrom_II, @chrom_III, @chrom_IV, @chrom_V, @chrom_X, 
    $units, $bp_length_per_unit,
);

my $rev_phys =0;

@chroms=qw(I II III IV V X);
#@chroms=qw(I);
foreach $chrom (@chroms){
  #print "$chrom -> @{$chrom_pos{$chrom}}\n";
  $parts=scalar @{$chrom_pos{$chrom}};
  #print $parts, "\n";
  for (my $i=0; $i< $parts; $i=$i+4){
    push(@pos_order, @{$chrom_pos{$chrom}}->[$i]); #gmap
    my $pos = @{$chrom_pos{$chrom}}->[$i];  # gmap
    my $mean_coord = @{$chrom_pos{$chrom}}->[$i+1];
    my $locus = @{$chrom_pos{$chrom}}->[$i+2];
    $cds = @{$chrom_pos{$chrom}}->[$i+3];
    push(@{$pos_order_to_mean_coord_locus_cds{$pos}}, $mean_coord, $locus, $cds); # key is gmap
    @pos_order = sort {$a <=>$b} @pos_order;  # sorting gmap positions of markers from left tip to right tip
  }      
  #print "@pos_order\n";
  my %seen=();
  foreach my $e (@pos_order){
    push(@unique_pos_order, $e) unless $seen{$e}++;
  }
  #print "@unique_pos_order#####\n";

  ###########################################
  # Left/Right tip gmap
  # Total units of a chromosom
  ###########################################

  $units = ($unique_pos_order[-1] - $unique_pos_order[0]);
  #print "Chromosome $chrom:\tLEFT tip\t$unique_pos_order[0] - RIGHT tip\t$unique_pos_order[-1]\t-> Total: $units Units\n";

  #####################################################################################
  # divide physical length of each chromosome by gmap units of corresponding chromosome
  # $bp_length_per_unit is the baseline for interpolated gmap outside tip markers 
  #####################################################################################
   
  $bp_length_per_unit = $chrom_length{$chrom} / $units;
  #print "Chrom_".$chrom.": $bp_length_per_unit\tbp/U\n";

  ###############################################
  # get gmap for CDSes identified as gmap markers
  ###############################################

  foreach my $gmap (@unique_pos_order){
    $mean_coord = @{$pos_order_to_mean_coord_locus_cds{$gmap}}->[0];
    $cds = @{$pos_order_to_mean_coord_locus_cds{$gmap}}->[2];
    if (exists $cds_mean{$cds} && $cds_mean{$cds} == $mean_coord){
      if (exists $CDS_variants{$cds}){
        foreach (@{$CDS_variants{$cds}}){
         # print  "\nSequence : \"$_\"\n";
         #print  "Interpolated_gMap\t$chrom\t$gmap\n";
         #print  "\nLocus : \"$predicted_gene_to_locus{$_}\"\n";
        #  print  "Map\t\"$chrom\"\tPosition\t$gmap\n";
        }
      }
      else {
       # print  "\nSequence : \"$cds\"\n";  
        #print  "Interpolated_gMap\t$chrom\t$gmap\n";
        #print  "\nLocus : \"$predicted_gene_to_locus{$_}\"\n";
      #  print  "Map\t\"$chrom\"\tPosition\t$gmap\n";
      }
    }
  }
  
  ##############################################################
  # get interpolated gmap for CDSes lying outside of tip markers
  ##############################################################
  
  my $outside=0;
  my $L_tip = $unique_pos_order[0];  
  my $R_tip = $unique_pos_order[-1];
  my $L_mean_coord = @{$pos_order_to_mean_coord_locus_cds{$L_tip}}->[0];
  my $R_mean_coord = @{$pos_order_to_mean_coord_locus_cds{$R_tip}}->[0];
  my $L_end = @{$All_Coords{$chrom}}->[0];
  my $R_end = @{$All_Coords{$chrom}}->[1];
  print "Chrom $chrom: $L_tip -> $R_tip -> $L_mean_coord -> $R_mean_coord\n";

  foreach $ea (sort keys %cds_mean){
    if ($cds_mean{$ea} < $L_mean_coord){
      #print "Chrom $chrom: $L_tip -> $R_tip -> $L_mean_coord -> $R_mean_coord -> LE: $L_end -> RE: $R_end\n";
      #print "$cds_mean{$ea} < $L_mean_coord\n"; 
      $outside++;
      my $length_diff = $cds_mean{$ea} - $L_end;
      #print "LD: $length_diff\n";
      my $position = $length_diff / $bp_length_per_unit;
      $position = $L_tip - $position;
      ace_output($ea, $chrom, $position);
    }

    if ($cds_mean{$ea} > $R_mean_coord){ 
      #print "$cds_mean{$ea} > $R_mean_coord\n";
      $outside++;
      my $length_diff = $cds_mean{$ea} - $R_end;  
      my $position = $length_diff / $bp_length_per_unit;
      $position = $R_tip + $position;
      ace_output($ea, $chrom, $position);
    }
  }
  print "Outside: $outside\n";
    
  ############################################################
  # get interpolated gmap for CDSes lying in between 2 markers
  ############################################################

  my ($baseline, $length_diff, $down_mean_coord, $up_mean_coord, $gmap_down, $gmap_up);
   
  foreach $ea (sort keys %cds_mean){
    for ($i=0; $i < (scalar @unique_pos_order) - 1; $i++){

      $gmap_down = $unique_pos_order[$i];  
      $gmap_up =   $unique_pos_order[$i+1];
      $down_mean_coord = @{$pos_order_to_mean_coord_locus_cds{$gmap_down}}->[0]; 
      $up_mean_coord =   @{$pos_order_to_mean_coord_locus_cds{$gmap_up}}->[0];
      $gmap_down = $unique_pos_order[$i];
      $gmap_up =   $unique_pos_order[$i+1];
     # print "$down_mean_coord  ->  $up_mean_coord\n";
    
      if ($cds_mean{$ea} > $down_mean_coord && $cds_mean{$ea} < $up_mean_coord){
        $length_diff = $cds_mean{$ea} - $down_mean_coord;

        ################################
        # negative gmap VS negative gmap
        ################################
        
        if ($gmap_down < 0 && $gmap_up <= 0){
          $baseline=($up_mean_coord - $down_mean_coord) / (-$gmap_down + $gmap_up);
          
        } 
        ################################
        # positive gmap VS positive gmap
        ################################
  
        if ($gmap_down >= 0 && $gmap_up > 0){
          $baseline=($up_mean_coord - $down_mean_coord) / ($gmap_up - $gmap_down);
        }
        ################################
        # negative gmap VS positive gmap
        ################################
    
        if ($gmap_down < 0 && $gmap_up > 0){
          $baseline=($up_mean_coord - $down_mean_coord) / ($gmap_up - $gmap_down);
        }
  
        my $position = $length_diff / $baseline;
        $position = $gmap_down + $position;
        ace_output($ea, $chrom, $position);       
      }
    }
  }  

  #####################################################
  # output list of gmap position and corresponding info
  #####################################################
        
  for ($i=0; $i<scalar @unique_pos_order; $i++){
    my $mean_coord_down = @{$pos_order_to_mean_coord_locus_cds{$unique_pos_order[$i]}}->[0];
    my $mean_coord_up;
    if ($i != (scalar @unique_pos_order)-1){
       $mean_coord_up =   @{$pos_order_to_mean_coord_locus_cds{$unique_pos_order[$i+1]}}->[0];
    }
    my $locus =           @{$pos_order_to_mean_coord_locus_cds{$unique_pos_order[$i]}}->[1];
    my $cds =             @{$pos_order_to_mean_coord_locus_cds{$unique_pos_order[$i]}}->[2];
    my $gmap =            $unique_pos_order[$i];
    print OUT "$chrom\t$gmap\t$locus\t$cds\t$mean_coord_down\n";
          
    #################################################
    # command line option to output reverse physicals
    #################################################
        
    if ($reverse){
      if ($i != (scalar @unique_pos_order)-1){
        #print "$mean_coord_down > $mean_coord_up\n";
        if ($mean_coord_down > $mean_coord_up){
          if (exists $CDS_variants{$cds}){
            $rev_phys++;
            print LOG "Reverse physical: $chrom\t$gmap\t$locus\t$cds \[@{$CDS_variants{$cds}}\]\t$mean_coord_down\n";
          }
          else {
            $rev_phys++;
            print LOG "Reverse physical: $chrom\t$gmap\t$locus\t$cds\t$mean_coord_down\n";
          }
        }  
      }    
    }      
  }
  if ($reverse){print LOG "There are $rev_phys reverse physicals on Chromosom $chrom\n"}
  @pos_order=();
  @unique_pos_order=();
  %pos_order_to_mean_coord_locus_cds=();
  $rev_phys=0;
}

my $end=`date +%H:%M:%S`; chomp $end;

print "\nJob started at $start\n";
print "Job finished at $end\n";

#################################
# s  u  b  r  o  u  t  i  n  e  s
#################################
          
sub dataset {
    
  ##############################################################
  # use latest GFF_SPLITT folder and retrieve relevant gff files
  ##############################################################
    
  my ($dir, $query)= @_;
  opendir(DIR, $dir) || die "Can't read directory";
  my (@files, @vers);
  my @dir=readdir DIR;
    
  splice(@dir, 0,2);
  closedir (DIR);
        
  if ($query eq "folder"){
    foreach (@dir){
      if ($_ =~ /^WS(\d+)/){
        push(@vers, $1);
      }
    }
    return @vers;
  }    
  
  if ($query eq "genes"){
    foreach (@dir){
      if ($_ =~ /^CHROMOSOME_([A-Z]+).genes.gff/){
        push(@files, $&);
      }   
    }
    return @files;
  }
  
  if ($query eq "rna"){
    foreach (@dir){
      if ($_ =~ /^CHROMOSOME_([A-Z]+).rna.gff/){
        push(@files, $&);
      }
    }
    return @files;
  }
}

sub ace_output {

  #########################################################
  # write acefile for all CDSes with gmap/interpolated gmap
  # also write Map position to locus objects
  #########################################################

  my ($cds, $chrom, $gmap)=@_;

  $gmap=sprintf("%6.6f", $gmap);

  if (exists $CDS_variants{$cds}){
    foreach (@{$CDS_variants{$cds}}){
      print ACE "\nSequence : \"$_\"\n";
      print ACE "Interpolated_gMap\t$chrom\t$gmap\n";
      if (exists $predicted_gene_to_locus{$_}){
        print ACE "\nLocus : \"$predicted_gene_to_locus{$_}\"\n";
        print ACE "Map\t\"$chrom\"\tPosition\t$gmap\n";
      }  
    }
  }
  else {
    print ACE "\nSequence : \"$cds\"\n";
    print ACE "Interpolated_gMap\t$chrom\t$gmap\n";
    if (exists $predicted_gene_to_locus{$cds}){
      print ACE "\nLocus : \"$predicted_gene_to_locus{$cds}\"\n";
      print ACE "Map\t\"$chrom\"\tPosition\t$gmap\n";
    }
  }
}
