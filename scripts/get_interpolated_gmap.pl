#!/usr/local/bin/perl5.6.1 -w

# get_interpolated_gmap.pl

# by Chao-Kung Chen [030512]

# This script calculates interpolated gmap for CDS/transcripts lying between as well as outside genetic markers.
# Output ace file of such information and upload to autoace during each build
# Output also other files related. See POD

# Last updated on: $Date: 2003-07-03 14:00:02 $
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

my ($diff, $reverse, $database, $gff_location, $help, $noUpload, $map, 
    $comp, $gmap_marker_dir, $gff_dir, $curr_db, $output);

GetOptions ("diff"          => \$diff,
            "rev|reverse"   => \$reverse,
	    "db|database=s" => \$database,
	    "map"           => \$map,
	    "comp"          => \$comp, 
	    "h|help"        => \$help, 
	    "n|noUpload"    => \$noUpload
           );

$gff_dir = "/wormsrv2/autoace/GFF_SPLITS/";
$curr_db = "/nfs/disk100/wormpub/DATABASES/current_DB/"; 
$output = "/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP";

if (!defined @ARGV or $help){system ("perldoc /wormsrv2/scripts/get_interpolated_gmap.pl"); exit(0)}

my @versions=dataset($gff_dir, "folder");

my @order = sort {$a <=> $b} @versions;

if (!defined $database){
  $gff_location = "/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS";
  $gmap_marker_dir = "/wormsrv2/autoace"; 
  print "\nUsing $gmap_marker_dir as database path for genetics marker loci . . .\n";
} 

if (defined $database){
  $gff_location = "/wormsrv2/autoace/GFF_SPLITS/WS"."$order[-1]";
  $gmap_marker_dir = $database;
  print "\nUsing $gmap_marker_dir as database path for genetics marker loci . . .\n";
}

my ($revfile, $diffile);

if ($reverse){
  my $version = $order[-1]+1;
  $revfile = "$output/reverse_physicals_WS$version.$rundate.$$";
  open(REV, ">$revfile") || die $!;
  system("chmod 777 $revfile");
}

if ($diff){
  $diffile = "/wormsrv2/logs/mapping_diff.".$rundate;
  system("rm -f $diffile; chmod 777 $diffile");
  open(DIFF, ">$diffile") || die $!;
}

chdir $gff_location;

my (@data, %CDS_mapping, %CDS_isoforms_mapping, %CDS_variants, %predicted_gene_to_locus);
my (%chrom_meanCoord, %genetics_mapping);
my ($cds, $parts, @coords, $i, $mean_coords, %cds_mean, %mean_coord_cds, %chrom_mean_coord_cds);

my $version = $order[-1]+1; # autoace
my $acefile = "$output/interpolated_map_"."WS$version.$rundate.ace";
my $gacefile = "$output/interpolated_map_to_geneace_"."WS$version.$rundate.ace";
my $download = "$output/"."WS$version"."_interpolated_map.txt";

if ($map){
  open (ACE, ">$acefile") || die "Can't output file!\n";
  open (GACE, ">$gacefile") || die "Can't output file!\n";
  open (DNLD, ">$download") || die "Can't output file!\n";
  system("chmod 777 $acefile $gacefile $download");
}

######################################################
# get list of predicted CDS/Transcript linked to locus
######################################################

my $CDS_linked_to_locus=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/get_predicted_gene_to_locus.def" quit
EOF
my $transcript_linked_to_locus=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/get_transcript_to_locus.def" quit
EOF

open (FH1, "echo '$CDS_linked_to_locus' | tace $curr_db | ") || die "Couldn't access $curr_db\n";
open (FH2, "echo '$transcript_linked_to_locus' | tace $curr_db | ") || die "Couldn't access $curr_db\n";         


while (<FH1>){
  chomp($_);
  if ($_ =~ /^\".+\"$/){
    my ($predict, $locus)= split(/\s+/, $_);
    $predict =~ s/\"//g; $locus =~ s/\"//g;
    $predicted_gene_to_locus{$predict} = $locus;
  }
}
while (<FH2>){
  chomp($_);
  if ($_ =~ /^\".+\"$/){
    my ($predict, $locus)= split(/\s+/, $_);
    $predict =~ s/\"//g; $locus =~ s/\"//g;
    $predicted_gene_to_locus{$predict} = $locus;
  }
}
close FH1; close FH2;

my %locus_map;

if($comp){
  open(IN, "/wormsrv1/geneace/gMAP/interpolated_gmap_based_on_contig_map_out.all") || die $!;
  open (MAPCOMP, ">$output/compare_gmap_WS$version.$rundate") || die $!;
  system("chmod 777 $output/compare_gmap_WS*.$rundate");
  while (<IN>){
    chomp;
    my ($locus, $LG, $map) = split (/\s+/, $_);
    $locus_map{$locus} = $map;
  }
  close IN; 
}

open (INFO, ">$output/gmap_info_"."WS$version.$rundate") || die $!;
system("chmod 777 $output/gmap_info_WS*.$rundate");

#print INFO "\nAll CDS/transcript linked to locus in $curr_db ", scalar keys %predicted_gene_to_locus,"\n\n";

########################################################################################
# retrieve data from gff file: CHROMOSOME_number.genes.gff and CHROMOSOME_number.rna.gff
########################################################################################

my @gff_files_cds=dataset($gff_location, "genes");
my @gff_files_rna=dataset($gff_location, "rna");

my $cds_count=0;
my $rna_count=0;
my $variants=0;
  
foreach (@gff_files_cds){
  @data = `less $_ | cut -f 1,4,5,9`;
  foreach (@data){
    chomp;
    if ($_ =~ /^CHROM.+/){
      $cds_count++;
      my ($chrom, $left, $right, $junk1, $CDS)= split(/\s+/,$_);
      $chrom =~ s/CHROMOSOME_//;
      $CDS =~ s/\"//g;
      my $ori = $CDS;
      if ($CDS =~ /(.+\.\d+)\D+/){
        $CDS = $1;
	$variants++;
        push(@{$CDS_variants{$CDS}}, $ori); #parent of isoforms
        push(@{$CDS_isoforms_mapping{$CDS}}, $chrom, $left, $right, "DNA");
      }
      else {
        push(@{$CDS_mapping{$CDS}}, $chrom, $left, $right, "DNA");
      }
    }
  }
}

print "\nEgrepping transcripts in gff files . . .\n"; 

foreach (@gff_files_rna){
  @data = `egrep "(RNA|UNKNOWN).+(Transcript|Sequence).+" $_ | cut -f 1,2,4,5,9`;
  #@data = `egrep ".+(Transcript|Sequence).+" $_ | cut -f 1,2,4,5,9`;
  foreach (@data){
    chomp;
    if ($_ !~ /\[.+\].+/ && $_ !~ /mRNA.+/){
      $rna_count++;
      my ($chrom, $type, $left, $right, $junk, $RNA)= split(/\s+/,$_);
      $RNA =~ s/(.+):.+/$1/;
      $RNA =~ s/\"//g;
      my $ori = $RNA;  
      $chrom =~ s/CHROMOSOME_//;
      if ($RNA =~ /(.+\.\d+)\D+/){
	$RNA = $1;
	$variants++;
	push(@{$CDS_variants{$RNA}}, $ori); #parent of isoforms
        push(@{$CDS_isoforms_mapping{$RNA}}, $chrom, $left, $right, "RNA");
      }
      else {
	push(@{$CDS_mapping{$RNA}}, $chrom, $left, $right, "RNA");
      }	
    }      
  }
}

print INFO "Things fetched from gff files:\n";
print INFO "Number of CDS/transcript without isoforms: ", scalar keys %CDS_mapping, "\n";
print INFO "Number of CDS/transcript with isoforms: ", scalar keys %CDS_isoforms_mapping, "\n";
print INFO "Number of CDS/transcript variants: ", $variants, "\n";
print INFO "Total CDS/transcript: ", $cds_count + $rna_count, "\n";
print "Processing gmap positions and coordinates . . .\n";
     
###########################################################################################
# get mean value of chrom coords of each CDS/Transcript from GFF files (genes|rna|rest.gff)
# mean of isoforms is the mean of the left and right coords of the longest variant
###########################################################################################
     
my $count=0;
my $chromo;

####################
# cds has no isoform
####################

foreach (sort keys %CDS_mapping){
  
  $mean_coords = (@{$CDS_mapping{$_}}->[1] + @{$CDS_mapping{$_}}->[2]) / 2;
  push(@{$cds_mean{$_}}, $mean_coords, @{$CDS_mapping{$_}}->[3]); # key: cds(w/0 isoform) value: mean coords, "DNA" or "RNA"
  $chromo = @{$CDS_mapping{$_}}->[0];
  push (@{$chrom_mean_coord_cds{$chromo}}, $mean_coords, $_);

}

##################
# cds has isoforms
##################

foreach (sort keys %CDS_isoforms_mapping){
  $parts = scalar @{$CDS_isoforms_mapping{$_}};

  for ($i=0; $i < $parts; $i=$i+4 ){
    push(@coords, ${@{$CDS_isoforms_mapping{$_}}}[$i+1], ${@{$CDS_isoforms_mapping{$_}}}[$i+2]);  # get coords of all isoforms
  }
  @coords = sort {$a <=> $b} @coords;
  $mean_coords = ($coords[-1] + $coords[0]) / 2;
  push(@{$cds_mean{$_}}, $mean_coords, @{$CDS_isoforms_mapping{$_}}->[3]);# key: cds(w/ isoform) value: mean coords

  $chromo = @{$CDS_isoforms_mapping{$_}}->[0];
  push (@{$chrom_mean_coord_cds{$chromo}}, $mean_coords, $_);
  @coords =();
}

##################################################################################
# retrieve locus data: chromosome, map position, genomic_sequence/transcript, gmap
##################################################################################
 
my $marker_gmap_of_each_chrom=<<EOF;
Table-maker -p "/wormsrv1/geneace/wquery/marker_gmap_of_each_chrom.def" quit
EOF
  
my ($mean_coord, %chrom_pos);
my $cmp_file = "$output/cmp_gmap_with_coord_order_"."WS$version.$rundate.$$";

if ($reverse) {
  open(CMP, ">$cmp_file") || die $!;
  system ("chmod 755 $cmp_file");
  print CMP "\nGFF file version from $gff_location\n\n";
}

##############################
# including CDS and Transcript 
############################## 

open (FH, "echo '$marker_gmap_of_each_chrom' | tace $gmap_marker_dir | ") || die "Couldn't access $gmap_marker_dir\n";
while (<FH>){
  chomp($_);
  # $4=cds/transcript, $2=chrom, $1=locus, $3=gmap position
  if ($_ =~ /^\"(.+)\"\s+\"(.+)\"\s+(-\d+\.\d+|\d+\.\d+|\d+|-\d+)\s+\"(.+)\"/){   # 0 & postive & negative values
    my $locus = $1;
    my $chrom = $2;
    my $gmap = $3;
    my $cds = $4;
   
    ##################
    # CDS has isoforms
    ##################
   
    if ($cds =~ /(.+\.\d+)\D+/){ 
      my $cds = $1;
      $mean_coord=@{$cds_mean{$cds}}->[0];
      if (defined $mean_coord){push(@{$genetics_mapping{$cds}}, $chrom, $locus, $gmap, $mean_coord)}
      else {push(@{$genetics_mapping{$cds}}, $chrom, $locus, $gmap, "NA")}
      if (defined $mean_coord){push(@{$chrom_pos{$chrom}}, $gmap, $mean_coord, $locus, $cds)}
      else {push(@{$chrom_pos{$chrom}}, $gmap, "NA", $locus, $cds)}
    }
  
    #####################
    # CDS has no isoforms
    #####################

    else {
      $mean_coord=@{$cds_mean{$cds}}->[0];
      if (defined $mean_coord){push(@{$genetics_mapping{$cds}}, $chrom, $locus, $gmap, $mean_coord)}
      else {push(@{$genetics_mapping{$cds}}, $chrom, $locus, $gmap, "NA")}
      if (defined $mean_coord){push(@{$chrom_pos{$chrom}}, $gmap, $mean_coord, $locus, $cds)}
      else {push(@{$chrom_pos{$chrom}}, $gmap, "NA", $locus, $cds)}	    
    }
  }
}
close FH;

$count = 0;
if ($diff){

  #############################################
  # check CDS mapping by coords and by genetics
  #############################################  
  
  foreach (sort keys %genetics_mapping){
    
    ###############################
    # check for CDS having isoforms 
    ###############################
    
    if (exists @{$CDS_isoforms_mapping{$_}}->[0] && (@{$genetics_mapping{$_}}->[0] ne @{$CDS_isoforms_mapping{$_}}->[0])){
      $count++;
      print  DIFF "ERROR: $_ (@{$genetics_mapping{$_}}->[1]) on @{$genetics_mapping{$_}}->[0] (genetics) \/ on @{$CDS_isoforms_mapping{$_}}->[0] (coordinates)\n\n";
    }
    
    ##################################################    
    # check for transcripts and CDS having no isoforms 
    ##################################################
    
    if (exists @{$CDS_mapping{$_}}->[0] && (@{$genetics_mapping{$_}}->[0] ne @{$CDS_mapping{$_}}->[0])){
      $count++;
      print  DIFF "ERROR: $_ (@{$genetics_mapping{$_}}->[1]) on @{$genetics_mapping{$_}}->[0] (genetics) \/ on @{$CDS_mapping{$_}}->[0] (coordinates)\n\n";
    }
  }
  print DIFF "There are $count discrepancies in genetics\/chrom. coords mapping\n";
  close DIFF;
  exit(0);
}

##################################################################
# hard coded chrom length - as change in the future is very seldom
##################################################################
    
my %chrom_length = (I => 15080471, II => 15279300, III => 13783268, IV => 17493790, V => 20922239,  X => 17705013);
foreach (sort keys %chrom_length){print INFO "$_ -> $chrom_length{$_} bp\n"}

##################################################################
# sorting gmap positions of marker loci from left tip to right tip 
##################################################################

print "Generating ace files . . .\n";

my (@chroms, $chrom, $ea, @pos_order, $part, @unique_pos_order, %pos_order_to_mean_coord_locus_cds, 
    @mean_coords_order, $units, $bp_length_per_unit, @all_mean_of_each_chrom, @unique_all_mean_of_each_chrom,
    %all_mean_of_each_chrom);

@chroms=qw(I II III IV V X);
foreach $chrom (@chroms){
  $parts=scalar @{$chrom_pos{$chrom}};
  for (my $i=0; $i< $parts; $i=$i+4){
    push(@pos_order, @{$chrom_pos{$chrom}}->[$i]);  # $i=gmap, duplication due to isoforms
    my $pos = @{$chrom_pos{$chrom}}->[$i];          
    my $mean_coord = @{$chrom_pos{$chrom}}->[$i+1]; 
    my $locus = @{$chrom_pos{$chrom}}->[$i+2];      
    $cds = @{$chrom_pos{$chrom}}->[$i+3];           
    push(@{$pos_order_to_mean_coord_locus_cds{$pos}}, $mean_coord, $locus, $cds); # key is gmap
    push(@all_mean_of_each_chrom, $mean_coord); # duplication due to isoforms
  }

  @pos_order = sort {$a <=>$b} @pos_order;  # sorting gmap positions of markers from left tip to right tip      
  
  my %seen=();
  foreach (@pos_order){push(@unique_pos_order, $_) unless $seen{$_}++}
  foreach (@all_mean_of_each_chrom){push(@unique_all_mean_of_each_chrom, $_) unless $seen{$_}++}
  foreach (@unique_all_mean_of_each_chrom){$all_mean_of_each_chrom{$_}++}

  print INFO "\nUnique markers $chrom: ", scalar @unique_pos_order,"\n";

  ###########################################
  # Left/Right tip gmap
  # Total units of a chromosom
  ###########################################

  $units = ($unique_pos_order[-1] - $unique_pos_order[0]);
  print INFO "Chromosome $chrom:\tLEFT tip\t$unique_pos_order[0] - RIGHT tip\t$unique_pos_order[-1]\t-> Total: $units Units\n";

  #####################################################################################
  # divide physical length of each chromosome by gmap units of corresponding chromosome
  # $bp_length_per_unit is the baseline for interpolated gmap outside tip markers 
  #####################################################################################
   
  $bp_length_per_unit = $chrom_length{$chrom} / $units;
  print INFO "Chrom_".$chrom.": $bp_length_per_unit\tbp/U\n";

  ##############################################################
  # get interpolated gmap for CDSes lying outside of tip markers
  ##############################################################

  my $outside_L=0;
  my $outside_R=0;
  my $L_tip = $unique_pos_order[0];
  my $R_tip = $unique_pos_order[-1];

  my $L_mean_coord = @{$pos_order_to_mean_coord_locus_cds{$L_tip}}->[0];
  my $R_mean_coord = @{$pos_order_to_mean_coord_locus_cds{$R_tip}}->[0];

  ##################################################################
  ## get end coords of each chromosom -> 
  ## for calculating int. gmap of CDSes lying outside of tip markers
  ##################################################################

  my $rev_phys =0;
  my (@all_coords, @all_cds_each_chrom);
  my ($baseline, $length_diff, $position, $down_mean_coord, $up_mean_coord, $gmap_down, $gmap_up);

  my @all_values = @{$chrom_mean_coord_cds{$chrom}};

  for (my $i=0; $i< scalar @all_values; $i=$i+2){
    push(@all_coords,  $all_values[$i]);
    $mean_coord_cds{$all_values[$i]}= $all_values[$i+1];
  }
  for (my $i=0; $i< scalar @all_values; $i=$i+2){push(@all_cds_each_chrom, $all_values[$i+1])}

  @all_coords = sort {$a <=>$b} @all_coords; # all coords of CDS/transcript on each chrom
  my $L_end = $all_coords[0];  # shortes mean coord of CDS/transcript on each chrom
  my $R_end = $all_coords[-1]; # longest mean coord of CDS/transcript on each chrom

  print INFO "Chrom $chrom: $L_tip -> $R_tip | $L_mean_coord -> $R_mean_coord | $L_end to $R_end\n";

  foreach (@all_coords){
 
   my $feature = @{$cds_mean{$mean_coord_cds{$_}}}->[1];
     
    if ($_ ne "NA" && $_ < $L_mean_coord){
      $outside_L++;
      $length_diff = $L_mean_coord - $_;
      $position = $length_diff / $bp_length_per_unit;
      $position = $L_tip - $position;
      if ($map){
	ace_output($mean_coord_cds{$_}, $chrom, $position, $_, $feature, \%predicted_gene_to_locus);
	next;
      }	
    }

    if ($_ ne "NA" && $_ > $R_mean_coord){ 
      $outside_R++;
      $length_diff = $_ - $R_mean_coord;  
      $position = $length_diff / $bp_length_per_unit;
      $position = $R_tip + $position;
      if ($map){
	ace_output($mean_coord_cds{$_}, $chrom, $position, $_, $feature, \%predicted_gene_to_locus);
	next; 
      }	
    }
  
    ############################################################
    # get interpolated gmap for CDSes lying in between 2 markers
    ############################################################
  
    if (($_ ne "NA") && $_ > $L_mean_coord && $_ < $R_mean_coord){
      
      for ($i=0; $i < (scalar @unique_pos_order) - 1; $i++){
         
        $gmap_down = $unique_pos_order[$i];
        $gmap_up =   $unique_pos_order[$i+1];
        $down_mean_coord = @{$pos_order_to_mean_coord_locus_cds{$gmap_down}}->[0];
        $up_mean_coord = @{$pos_order_to_mean_coord_locus_cds{$gmap_up}}->[0];

        if (($down_mean_coord ne "NA" && $up_mean_coord ne "NA") && 
           ($down_mean_coord < $_ && $_ < $up_mean_coord)){
          $length_diff = $_ - $down_mean_coord;

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
          $position = $length_diff / $baseline;
          $position = $gmap_down + $position;
	  if (!$all_mean_of_each_chrom{$_} && $map){
            ace_output($mean_coord_cds{$_}, $chrom, $position, $_, $feature, \%predicted_gene_to_locus);       
          }
          last;        
        }
      }
    } 
  }
  print INFO "There are ", scalar @all_coords, " CDS/transcript (w/o counting isoforms) on $chrom:\n";
  print INFO "Outside: L-end ($outside_L) R-end ($outside_R). Total: ",$outside_L + $outside_R,"\n";

  #######################################################
  # output:  list of gmap position and corresponding info
  #          list of reverse physicals
  #######################################################

  my $locusf;

  if ($reverse){
    
    for ($i=0; $i < (scalar @unique_pos_order)-1; $i++){      
      
      $gmap_down = $unique_pos_order[$i];
      my $gmdn =  $gmap_down;
      $gmdn = sprintf("%8.4f", $gmdn);
      $gmap_up =   $unique_pos_order[$i+1];
      my $locus1 = @{$pos_order_to_mean_coord_locus_cds{$gmap_down}}->[1];
      $locusf = $locus1;
      $locusf = sprintf("%10s", $locusf);
      $down_mean_coord = @{$pos_order_to_mean_coord_locus_cds{$gmap_down}}->[0];
      my $dnmcoord =  $down_mean_coord;
      $dnmcoord = sprintf("%12.1f", $dnmcoord);
      $up_mean_coord = @{$pos_order_to_mean_coord_locus_cds{$gmap_up}}->[0];
      $cds = @{$pos_order_to_mean_coord_locus_cds{$gmap_down}}->[2];
      my $cdsf = $cds;
      $cdsf = sprintf("%-15s", $cdsf);

      print CMP "$chrom\t$gmdn\t$locusf\t$cdsf\t$dnmcoord\n";
      if ($down_mean_coord eq "NA") {print REV "\n** Gmap marker $cds ($locus1) on $chrom has no coordinate **\n"}

      if ($down_mean_coord ne "NA" && $up_mean_coord ne "NA" && ($down_mean_coord > $up_mean_coord)){
	if (exists $CDS_variants{$cds}){
	  $rev_phys++;
	  if ($reverse) {print REV "Reverse physical: $chrom\t$gmap_down\t$locus1\t$cds\[@{$CDS_variants{$cds}}\]\t$down_mean_coord\n"}
	}
	else {
	  $rev_phys++;
	  if ($reverse) {print REV "Reverse physical: $chrom\t$gmap_down\t$locus1\t$cds\t$down_mean_coord\n"}
	}
      }  
    }  
  }  
  if ($reverse){
    print REV "--->$rev_phys reverse physical(s) on Chromosom $chrom\n\n";
    my $gmap = $unique_pos_order[-1];
    my $gm = $gmap;
    $gm = sprintf("%8.4f", $gm);
    my $locus = @{$pos_order_to_mean_coord_locus_cds{$gmap}}->[1];
    $locusf = $locus; 
    $locusf = sprintf("%10s", $locusf);
    $mean_coord = @{$pos_order_to_mean_coord_locus_cds{$gmap}}->[0];
    my $mc = $mean_coord;
    $mc = sprintf("%12.1f", $mc);
    $cds = @{$pos_order_to_mean_coord_locus_cds{$gmap}}->[2];
    my $cdsf = $cds;
    $cdsf = sprintf("%-15s", $cdsf);

    print CMP "$chrom\t$gm\t$locusf\t$cdsf\t$mc\n"; 
  }

  @pos_order=(); @unique_pos_order=(); @all_coords=(); %pos_order_to_mean_coord_locus_cds=(); 
  $rev_phys=(); @all_mean_of_each_chrom=(); @unique_all_mean_of_each_chrom=(); %all_mean_of_each_chrom=();
}

if (!$noUpload){
  print "\n\nDeleting old interpolated map and uploading new ones to autolace . . .\n\n";

  my $tace = &tace;
  my $log = "/wormsrv2/logs/load_gmap_to_autoace_"."WS$order[-1].$rundate.$$";
 
  my $command=<<END;
find sequence * where Interpolated_map_position
edit -D Interpolated_map_position
find transcript * where Interpolated_map_position
edit -D Interpolated_map_position
find locus * where Interpolated_map_position
edit -D Interpolated_map_position

pparse $acefile
save
quit
END

  open (Load_A,"| $tace /wormsrv2/autoace/ >> $log") || die "Failed to upload to Geneace";
  print Load_A $command;
  close Load_A;
 
  system("chmod 777 $log");
  
  print "\n\nDone with uploading . . .\n\n";
}

my $end=`date +%H:%M:%S`; chomp $end;

print "\nJob started at $start\n";
print "Job finished at $end\n";

my $maintainers = "ck1\@sanger.ac.uk";
mail_maintainer($0,$maintainers,"Interpolated map position for WS$version can now be loaded to geneace.");


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
        
  if ($query eq "folder"){foreach (@dir){if ($_ =~ /^WS(\d+)/){push(@vers, $1)}} return @vers}    
  if ($query eq "genes"){foreach (@dir){if ($_ =~ /^CHROMOSOME_(I|II|III|IV|V|X).genes.gff/ ||
					    $_ =~ /^CHROMOSOME_(I|II|III|IV|V|X).clone_path.gff/){push(@files, $&)}} return @files}
  if ($query eq "rna"){foreach (@dir){if ($_ =~ /^CHROMOSOME_(I|II|III|IV|V|X).(rna|rest).gff/){push(@files, $&)}} return @files}
}

sub ace_output {

  #########################################################
  # write acefile for all CDSes with gmap/interpolated gmap
  # also write Map position to locus objects
  #########################################################

  my ($cds, $chrom, $gmap, $mean_coord, $feature, $predicted_gene_to_locus)=@_;
 
  %predicted_gene_to_locus = %$predicted_gene_to_locus;
  
  $gmap = sprintf("%8.4f", $gmap);
  my $cdsf = $cds; 
  $cdsf = sprintf ("%-15s", "$cdsf");

  if (exists $CDS_variants{$cds}){

    foreach (@{$CDS_variants{$cds}}){
      
      my $cdsf = $_;
      $cdsf = sprintf ("%-15s", "$cdsf");
   
      if ($feature eq "DNA"){print ACE "\nSequence : \"$_\"\n"; print DNLD "$cdsf\t"}
      if ($feature eq "RNA"){print ACE "\nTranscript : \"$_\"\n"; print DNLD "$cdsf\t"}
      print ACE "Interpolated_map_position\t\"$chrom\"\t$gmap\t\/\/$mean_coord (iso)\n";
      print DNLD "\t$chrom\t$gmap\t";
      if (exists $predicted_gene_to_locus{$_}){
        print ACE "\nLocus : \"$predicted_gene_to_locus{$_}\"\n";
	print GACE "\nLocus : \"$predicted_gene_to_locus{$_}\"\n";
	print DNLD "$predicted_gene_to_locus{$_}\n";
        print ACE "Interpolated_map_position\t\"$chrom\"\t$gmap\t\/\/$mean_coord (iso)\n";
	print GACE "Interpolated_map_position\t\"$chrom\"\t$gmap\t\/\/$mean_coord (iso)\n";
	if ($comp && $locus_map{$predicted_gene_to_locus{$_}}){
	  print MAPCOMP "$predicted_gene_to_locus{$_}\t\"$chrom\"\tCoords:\t$gmap\tContig:\t$locus_map{$predicted_gene_to_locus{$_}}\n";
	}
      }
      else {print DNLD "-\n"}	
    }
  }  
  else {
    if ($feature eq "DNA"){print ACE "\nSequence : \"$cds\"\n"; print DNLD "$cdsf\t"}
    if ($feature eq "RNA"){print ACE "\nTranscript : \"$cds\"\n"; print DNLD "$cdsf\t"}
    print ACE "Interpolated_map_position\t\"$chrom\"\t$gmap\t\/\/$mean_coord\n";
    print DNLD "\t$chrom\t$gmap\t";
    if (exists $predicted_gene_to_locus{$cds}){
      print ACE "\nLocus : \"$predicted_gene_to_locus{$cds}\"\n";
      print GACE "\nLocus : \"$predicted_gene_to_locus{$cds}\"\n";
      print DNLD "$predicted_gene_to_locus{$cds}\n";
      print ACE "Interpolated_map_position\t\"$chrom\"\t$gmap\t\/\/$mean_coord\n";
      print GACE "Interpolated_map_position\t\"$chrom\"\t$gmap\t\/\/$mean_coord\n";
      if ($comp && exists $locus_map{$predicted_gene_to_locus{$cds}}){
	print MAPCOMP "$predicted_gene_to_locus{$cds}\t\"$chrom\"\tCoords:\t$gmap\tContig:\t$locus_map{$predicted_gene_to_locus{$cds}}\n";
      }
    }
    else {print DNLD "-\n"}
  }
}

__END__
=head2 NAME - get_interpolated_gmap.pl  

=head3 <SCRIPT DESCRIPTION> 

            This script is run at end of build to calculate interpolated map positions for all CDS and Transcripts.
            It can generate the following files:

            (1) /wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/reverse_physicals_WSXXX.yymmdd.pid
            (2) /wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/cmp_gmap_with_coord_order_WSXXX.yymmdd.pid (list of marker loci with gmap and coords)
            (3) /wormsrv2/logs/mapping_diff.yymmdd (not run for the build, but called by geneace_check.pl)  
            (4) /wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/interpolated_gmap_WSXXX.yymmdd 
            (5) /wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/interpolated_gmap_to_geneace_WSXXX.yymmdd.ace
            (6) /wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/WSXXX_interpolated_map.txt
            (7) /wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/gmap_info_WS*yymmdd (information might be interesting in case data looks strange) 
            

=head3 <USAGE> 

            Output reverse physicals of marker loci. Ie, check if the marker loci for genetics map are aligned lineally by gff coordinates.

B<-db: / -databse:>
   
            Specifies database to look for marker map positions and gff coordinates
            Eg. 
                -db /wormsrv1/geneace 

                if this option is omitted and this script is run during the build, it points to autoace
                                
B<-map:>    
            Output interpolated map positions as ace file. Requeires -rev to check for reverse physicals
           

B<-comp:>
             
            Output a list that compares interpolated gmap of CDS/transcripts links to locus based on 
            (1) sequence coordinates 
            (2) physical contig maps, which can be view in gmap


B <-n / -noUpload):
  
            Interpolated_map_positions will not be uploaded to autoace.
