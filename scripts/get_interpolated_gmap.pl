#!/usr/local/bin/perl5.8.0 -w
#
# get_interpolated_gmap.pl
#
# by Chao-Kung Chen [030512]
#
# This script calculates interpolated genetic map positions for CDS, Transcripts 
# and Pseudogenes lying between and outside genetic markers.
#
# Last updated on: $Date: 2004-09-02 12:04:18 $
# Last updated by: $Author: dl1 $


use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Cwd 'chdir';
use Getopt::Long;
use GENEACE::Geneace;

###################################################
# variables and command-line options with aliases #
###################################################

my ($diff, $reverse, $database, $gff_location, $help, $debug, $map, $comp, $verbose);

GetOptions ("diff"          => \$diff,
            "rev|reverse"   => \$reverse,
	    "database=s" => \$database,
	    "map"           => \$map,
	    "comp"          => \$comp,
	    "h|help"        => \$help,
	    "d|debug"       => \$debug,
	    "v|verbose"     => \$verbose,
           );


my $gff_dir    = "/wormsrv2/autoace/GFF_SPLITS/";
my $output     = "/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP";

my $rundate    = &rundate;
my $start      = &runtime;
my $script_dir = "/wormsrv2/scripts/";
my $tace = glob("~wormpub/ACEDB/bin_ALPHA/tace");

if (!defined @ARGV){system ("perldoc /wormsrv2/scripts/get_interpolated_gmap.pl"); exit(0)}

# set WS version number
my $version = &get_wormbase_version;


# Use specified database for path but default to using autoace if -database not specified
if($database){
  my $prev_version = $version -1;
  $gff_location = "/wormsrv2/autoace/GFF_SPLITS/WS"."$prev_version";
}
else{
  $gff_location = "/wormsrv2/autoace/GFF_SPLITS/GFF_SPLITS";
  $database = "/wormsrv2/autoace";
}
print "\nUsing $database as database path for genetics marker loci\n";

my ($revfile, $diffile);

if ($reverse){
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

my (@data, %CDS_mapping, %CDS_isoforms_mapping, %CDS_variants, %CDS_to_gene);
my (%chrom_meanCoord, %genetics_mapping);
my ($cds, $parts, @coords, $i, $mean_coords, %cds_mean, %mean_coord_cds, %chrom_mean_coord_cds);

######################
# path of output files
######################

my $acefile = "$output/interpolated_map_"."WS$version.$rundate.ace";
my $gacefile = "$output/interpolated_map_to_geneace_"."WS$version.$rundate.ace";
my $download1 = "$output/"."WS$version"."_CDSes_interpolated_map.txt";
my $download2 = "$output/"."WS$version"."_Clones_interpolated_map.txt";
my $download3 = "$output/"."WS$version"."_Pseudogene_interpolated_map.txt";


if ($map){
  if (defined glob("$output/WS*gz")){system("rm -f $output/WS*gz")} # remove zipped map file of last build

  open (ACE, ">$acefile") || die "Can't output file!\n";
  open (GACE, ">$gacefile") || die "Can't output file!\n";
  open (CDSes, ">$download1") || die "Can't output file!\n";
  open (CLONEs, ">$download2") || die "Can't output file!\n";
  open (PSEUDO, ">$download3") || die "Can't output file!\n";
  system("chmod 777 $acefile $gacefile $download1 $download2 $download3");

  print CDSes "WS$version interpolated map positions for CDSes\n\n";
  print CLONEs "WS$version interpolated map positions for Clones\n\n";
  print PSEUDO "WS$version interpolated map positions for Pseudogenes\n\n";
}

##############################################################
# hash to retrieve Gene obj info: 
#      eg. WBGenexxxxxxxx to CGC_name/Sequence_name/Other_name
#          and vice verse
##############################################################

my $ga = init Geneace();
my %Gene_info = $ga -> gene_info($database);

#####################################################################################################
# get list of CDS/Transcript/Pseudogene linked to Gene for writing interpolated map to it
#####################################################################################################


my $CDS_linked_to_gene=<<EOF;
  Table-maker -p "$database/wquery/get_CDS_to_gene.def" quit
EOF

my $transcript_linked_to_gene=<<EOF;
  Table-maker -p "$database/wquery/get_transcript_to_gene.def" quit
EOF

my $pseudogene_linked_to_gene=<<EOF;
  Table-maker -p "$database/wquery/get_pseudogene_to_gene.def" quit
EOF

open (FH1, "echo '$CDS_linked_to_gene'        | $tace $database |") || die "Couldn't access $database\n";
open (FH2, "echo '$transcript_linked_to_gene' | $tace $database |") || die "Couldn't access $database\n";
open (FH3, "echo '$pseudogene_linked_to_gene' | $tace $database |") || die "Couldn't access $database\n";

my %sequence_linked_to_gene;

my @FHS = qw(*FH1 *FH2 *FH3);
foreach my $e (@FHS){
  while (<$e>){
    chomp($_);
    if ($_ =~ /^\"(.+)\"\s+\"(.+)\"$/){
      my $seq = $1; my $gene_id = $2;
      my $locus;
      $locus = $Gene_info{$gene_id}{'Public_name'} if exists $Gene_info{$gene_id}{'Public_name'}; # make sure no null values
      $sequence_linked_to_gene{$seq} = $locus;
    }
  }
}

close FH1; close FH2; close FH3;

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


########################################################################################
# retrieve data from gff file: CHROMOSOME_number.genes.gff and CHROMOSOME_number.rna.gff
########################################################################################

my @gff_files_cds         = &dataset($gff_location, "CDS");
my @gff_files_rna         = &dataset($gff_location, "rna");
my @gff_files_pseudogene  = &dataset($gff_location, "pseudogene");

my $cds_count=0;
my $rna_count=0;
my $pseudogene_count=0;
my $variants=0;

print "\nParsing CDS coords from gff files . . .\n";

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
        push(@{$CDS_variants{$CDS}}, $ori);  #parent of isoforms
        push(@{$CDS_isoforms_mapping{$CDS}}, $chrom, $left, $right, "CDS");
      }
      else {
        push(@{$CDS_mapping{$CDS}}, $chrom, $left, $right, "CDS"); # for clones
      }
    }
  }
}

print "\nParsing transcript coords from gff files . . .\n"; 

foreach (@gff_files_rna){

  @data = `egrep "Transcript.+" $_ | cut -f 1,3-5,9`;

  foreach (@data){
    chomp;

    $rna_count++;
    my ($chrom, $type, $left, $right, $junk, $RNA)= split(/\s+/,$_);
    if ($type =~ /Transcript/i ){ # $type is the feature of transcript obj
      #print "$chrom -> $type -> $left -> $right -> $junk -> $RNA  \n";	
      $RNA =~ s/\"//g;
      my $ori = $RNA;
      $chrom =~ s/CHROMOSOME_//;
      if ($RNA =~ /(.+\.\d+)\D+/){
	$RNA = $1;
	$variants++;
	push(@{$CDS_variants{$RNA}}, $ori) if !exists $CDS_variants{$RNA}; #get rid of duplicates of CDS listed in rna.gff as coding_transcript obj
	push(@{$CDS_isoforms_mapping{$RNA}}, $chrom, $left, $right, "RNA");
      }
      else {
	push(@{$CDS_mapping{$RNA}}, $chrom, $left, $right, "RNA");
      }	
    }
  }
}


print "\nParsing pseudogene coords from gff files . . .\n";

foreach (@gff_files_pseudogene){
  @data = `less $_ | cut -f 1,4,5,9`;
  foreach (@data){
    chomp;
    if ($_ =~ /^CHROM.+/){
      $pseudogene_count++;
      my ($chrom, $left, $right, $junk1, $pseudogene)= split(/\s+/,$_);
      $chrom =~ s/CHROMOSOME_//;
      $pseudogene =~ s/\"//g;
      my $ori = $pseudogene;
      if ($pseudogene =~ /(.+\.\d+)\D+/){
        $pseudogene = $1;
	$variants++;
        push(@{$CDS_variants{$pseudogene}}, $ori); #parent of isoforms
        push(@{$CDS_isoforms_mapping{$pseudogene}}, $chrom, $left, $right, "pseudogene");
      }
      else {
        push(@{$CDS_mapping{$pseudogene}}, $chrom, $left, $right, "pseudogene");
      }
    }
  }
}

print "\nDoing interpolation . . .\n";

###########################################################################################
# get mean value of chrom coords of each CDS/Transcript from GFF files (genes|rna|rest.gff)
# mean of isoforms is the mean of the left and right coords of the longest variant
###########################################################################################

my $chromo;

####################
# cds has no isoform
####################

foreach (sort keys %CDS_mapping){

  $mean_coords = ($CDS_mapping{$_}->[1] + $CDS_mapping{$_}->[2]) / 2;
  push(@{$cds_mean{$_}}, $mean_coords, $CDS_mapping{$_}->[3]); # key: cds(w/0 isoform) value: mean coords, "CDS" or "RNA" or "pseudogene"
  $chromo = $CDS_mapping{$_}->[0];
  push (@{$chrom_mean_coord_cds{$chromo}}, $mean_coords, $_);

}

##################
# cds has isoforms
##################

foreach (sort keys %CDS_isoforms_mapping){
  $parts = scalar @{$CDS_isoforms_mapping{$_}};

  for ($i=0; $i < $parts; $i=$i+4 ){
    push(@coords, $CDS_isoforms_mapping{$_}->[$i+1], $CDS_isoforms_mapping{$_}->[$i+2]);  # get coords of all isoforms
  }
  @coords = sort {$a <=> $b} @coords;

  $mean_coords = ($coords[-1] + $coords[0]) / 2;

  push(@{$cds_mean{$_}}, $mean_coords, $CDS_isoforms_mapping{$_}->[3]);# key: cds(w/ isoform) value: mean coords
  $chromo = $CDS_isoforms_mapping{$_}->[0];
  push (@{$chrom_mean_coord_cds{$chromo}}, $mean_coords, $_);
  @coords =();
}

if ($debug){
  foreach (sort keys %chrom_mean_coord_cds){
    print "$_ => @{$chrom_mean_coord_cds{$_}} MEAN\n";
  }
}

##################################################################################
# retrieve marker loci (linked to CDS/Transcript/Pseudogene and have map position)
##################################################################################

my $marker_gmap_of_each_chrom=<<EOF;
Table-maker -p "/wormsrv1/geneace/wquery/marker_gmap_of_each_chrom.def" quit
EOF

my ($mean_coord, %chrom_pos);
my $cmp_file = "$output/cmp_gmap_with_coord_order_"."WS$version.$rundate.$$";

if ($reverse) {
  open(CMP, ">$cmp_file") || die $!;
  system ("chmod 755 $cmp_file");
  print CMP "\nGFF file version from WS$version\n\n";
}

open (FH, "echo '$marker_gmap_of_each_chrom' | $tace $database |") || die "Couldn't access $database\n";

my %unique_locus;
while (<FH>){
  chomp($_);
  # $4=cds/transcript/pseudogene, $2=chrom, $1=gene_id, $3=gmap position

  # sequence from same class and sequence from diff. classes (eg. CDS, Pseudogene)
  if ( ($_ =~ /^\"(WBGene\d+)\"\s+\"(.+)\"\s+(-\d+\.\d+|\d+\.\d+|\d+|-\d+)\s+\"(.+\.[a-z0-9]+)\"\s+\"(.+\.[a-z0-9]+)\"/) || 
       ($_ =~ /^\"(WBGene\d+)\"\s+\"(.+)\"\s+(-\d+\.\d+|\d+\.\d+|\d+|-\d+)\s+\"(.+\.[a-z0-9]+)\"/) ){ 

    my $gene_id = $1;
    my $locus = $Gene_info{$gene_id}{'Public_name'} if exists $Gene_info{$gene_id}{'Public_name'};
    my $chrom = $2;
    my $gmap = $3;
    my $cds = $4;
    my $cds_2 ;

    if (!$5){$cds_2 = "NA"}
    else{$cds_2 = $5}

    my @CDSES = (); push(@CDSES, $cds, $cds_2);

    foreach my $cds (@CDSES){
      if (!exists $unique_locus{$locus}){ # get rid of duplication of CGC_name from isoforms
	$unique_locus{$locus} = $locus;

	##################
	# CDS has isoforms
	##################
	
	if ($cds eq "NA"){}
	elsif ($cds =~ /(.+\.\d+)\D/){
	  my $cds = $1;
	  $mean_coord=$cds_mean{$cds}->[0];
	  if (defined $mean_coord){push(@{$genetics_mapping{$cds}}, $chrom, $locus, $gmap, $mean_coord)}
	  else {push(@{$genetics_mapping{$cds}}, $chrom, $locus, $gmap, "NA")}
	  if (defined $mean_coord){push(@{$chrom_pos{$chrom}}, $gmap, $mean_coord, $locus, $cds)}
	  else {push(@{$chrom_pos{$chrom}}, $gmap, "NA", $locus, $cds)}
	}
	
	#####################
	# CDS has no isoforms
	#####################
	
	else {			# blah.t1 (tRNA eg. goes here)
	  $mean_coord=$cds_mean{$cds}->[0];
	  if (defined $mean_coord){push(@{$genetics_mapping{$cds}}, $chrom, $locus, $gmap, $mean_coord)}
	  else {push(@{$genetics_mapping{$cds}}, $chrom, $locus, $gmap, "NA")}
	  if (defined $mean_coord){push(@{$chrom_pos{$chrom}}, $gmap, $mean_coord, $locus, $cds)}
	  else {push(@{$chrom_pos{$chrom}}, $gmap, "NA", $locus, $cds)}	
	}
      }
    }
  }
}

close(FH);
%unique_locus =();

##############################################
# check gene mapping by coords and by genetics
##############################################

if ($diff){
  my $count = 0;
  foreach (sort keys %genetics_mapping){

    # check for CDS/Pseudogene/Transcript? having isoforms 
    if (exists $CDS_isoforms_mapping{$_}->[0] && $genetics_mapping{$_}->[0] ne $CDS_isoforms_mapping{$_}->[0]){
      $count++;
      print  DIFF "ERROR: $_ $genetics_mapping{$_}->[1]) on $genetics_mapping{$_}->[0] (genetics) \/ on $CDS_isoforms_mapping{$_}->[0] (coordinates)\n\n";
    }

    # check for CDS/Pseudogene/transcripts having no isoforms 
    if (exists $CDS_mapping{$_}->[0] && $genetics_mapping{$_}->[0] ne $CDS_mapping{$_}->[0]){
      $count++;
      print  DIFF "ERROR: $_ $genetics_mapping{$_}->[1]) on $genetics_mapping{$_}->[0] (genetics) \/ on $CDS_mapping{$_}->[0] (coordinates)\n\n";
    }
  }
  print DIFF "There are $count discrepancies in genetics\/chromosome coordinates mapping\n";
  close DIFF;
  exit(0);
}

###################################################
# get chrom length - as slight changes still occurs
###################################################

my @dna_file = glob("/wormsrv2/autoace/CHROMOSOMES/CHROMOSOME_*.dna");

my %chrom_length;
foreach (@dna_file){
  my $line;
  my @line = `egrep "[atcg]" $_`;
  foreach (@line){chomp; $line .= $_}
  $_ =~ /.+CHROMOSOME_(\w+)\.dna/;
  my $chrom = $1;
  $chrom_length{$chrom} = length $line if $chrom ne "MtDNA";
}

foreach (sort keys %chrom_length){print INFO "$_ -> $chrom_length{$_} bp\n"}


#########################################################################
# list of all CDS / transcripts / pseudogenes:
# to distinguish later on which of these classes a seq variant belongs to
#########################################################################

my $db = Ace->connect(-path=>$database,
                      -program =>$tace) || do { print "Connection failure: ",Ace->error; die();}; 


my @genes= $db->fetch(-query=>'find worm_genes');

my %class;

foreach my $gene (@genes){
  if (defined $gene->at('Properties.Coding.CDS')){
    $class{$gene} = "CDS";
  }
  elsif (defined $gene->at('Properties.Transcript')){
    $class{$gene} = "RNA";
  }
  elsif (defined $gene->at('Type.Coding_pseudogene') || defined $gene->at('Type.RNA_pseudogene')){
    $class{$gene} = "pseudo";
  }
}
$db->close;


# for debugging only
if ($debug){
  open (CL, ">/wormsrv2/autoace/MAPPINGS/INTERPOLATED_MAP/type");
  foreach (sort keys %class){
    print CL "$_ -> $class{$_}  (1-1)\n";
  }
  foreach (keys %CDS_isoforms_mapping){
    print CL "$_ -> $CDS_isoforms_mapping{$_}->[3] PARENT (1-2)\n"; # [3] is type (DNA, RNA or psdudo)
  }
  foreach (keys %CDS_variants){
    print CL "$_ -> @{$CDS_variants{$_}} ISO (1-3)\n";
  }
}

######################################################################################################
# hash for checking a gene_id is live: use for writing interpolated_map_position only to live gene ids
######################################################################################################

my @id_status = $ga -> gene_id_status();  # array of 2 hash refs
my %gene_id_is_live = %{$id_status[0]};      # first one is about gene_id_is_live

##################################################################
# sorting gmap positions of marker loci from left tip to right tip
##################################################################

my (@chroms, $chrom, $ea, @pos_order, $part, @unique_pos_order, %pos_order_to_mean_coord_locus_cds, 
    @mean_coords_order, $units, $bp_length_per_unit, @all_mean_of_each_chrom, @unique_all_mean_of_each_chrom,
    %all_mean_of_each_chrom);

# check for any error which will mean script needs to be rerun
my $error_check = 0;
my $rev_phys =();
my $count_rev = 0;

print "\nGenerating ace files . . .\n";

@chroms=qw(I II III IV V X);
foreach $chrom (@chroms){

  $parts=scalar @{$chrom_pos{$chrom}};
  for (my $i=0; $i< $parts; $i=$i+4){
    push(@pos_order, $chrom_pos{$chrom}->[$i]);  # $i=gmap, has duplication if isoforms
    my $pos = $chrom_pos{$chrom}->[$i];
    my $mean_coord = $chrom_pos{$chrom}->[$i+1];
    my $locus = $chrom_pos{$chrom}->[$i+2];

    $cds = $chrom_pos{$chrom}->[$i+3];

    # multiple loci may have same gmap before correction 
    push(@{$pos_order_to_mean_coord_locus_cds{$pos}}, $mean_coord, $locus, $cds);
    push(@all_mean_of_each_chrom, $mean_coord); # duplication due to isoforms
  }

  @pos_order = sort {$a <=>$b} @pos_order;  # sorting gmap positions of markers from left tip (-) to right tip (+)

  # remove duplicated mean coord or gmap positions from isoforms
  my %seen=();
  foreach (@pos_order){push(@unique_pos_order, $_) unless $seen{$_}++}
  foreach (@all_mean_of_each_chrom){push(@unique_all_mean_of_each_chrom, $_) unless $seen{$_}++}
  foreach (@unique_all_mean_of_each_chrom){$all_mean_of_each_chrom{$_}++} # hash for quick lookup of all mean coords

  print INFO "\nUnique markers $chrom: ", scalar @unique_pos_order,"\n";   # only gmap positions

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

  ################################################################################################
  # start processing interpolated gmap of CDSes/transcript/pseudogene lying outside of tip markers
  ################################################################################################

  my $outside_L=0;
  my $outside_R=0;
  my $L_tip = $unique_pos_order[0];
  my $R_tip = $unique_pos_order[-1];

  my $L_mean_coord = $pos_order_to_mean_coord_locus_cds{$L_tip}->[0];
  my $R_mean_coord = $pos_order_to_mean_coord_locus_cds{$R_tip}->[0];

  ########################################################################################
  ## get end coords of each chromosom ->
  ## for calculating int. gmap of CDSes/transcript/pseudogene lying outside of tip markers
  ########################################################################################

  $rev_phys =0;
  my (@all_coords, @all_cds_each_chrom);
  my ($baseline, $length_diff, $position, $down_mean_coord, $up_mean_coord, $gmap_down, $gmap_up);

  my @all_values = @{$chrom_mean_coord_cds{$chrom}};

  for (my $i=0; $i< scalar @all_values; $i=$i+2){
    push(@all_coords,  $all_values[$i]);
    $mean_coord_cds{$all_values[$i]}= $all_values[$i+1];
  }
  for (my $i=0; $i< scalar @all_values; $i=$i+2){push(@all_cds_each_chrom, $all_values[$i+1])}

  @all_coords = sort {$a <=>$b} @all_coords; # all mean coords of CDS/transcript/pseudogene on each chrom
  my $L_end = $all_coords[0];  # shortes mean coord of CDS/transcript/pseudogene on each chrom
  my $R_end = $all_coords[-1]; # longest mean coord of CDS/transcript/pseudogene on each chrom

  print INFO "Chrom $chrom: $L_tip -> $R_tip | $L_mean_coord -> $R_mean_coord | $L_end to $R_end\n";

  foreach (@all_coords){

    # get interpolated map if mean coord is not that of a marker loci
    if (!$all_mean_of_each_chrom{$_}){ 

      my $feature = $cds_mean{$mean_coord_cds{$_}}->[1];

      if ($_ ne "NA" && $_ < $L_mean_coord){
        $outside_L++;
        $length_diff = $L_mean_coord - $_;
        $position = $length_diff / $bp_length_per_unit;
        $position = $L_tip - $position;
        if ($map){
	  &ace_output($mean_coord_cds{$_}, $chrom, $position, $_, $feature);
	  next;
        }	
      }

      if ($_ ne "NA" && $_ > $R_mean_coord){
        $outside_R++;
        $length_diff = $_ - $R_mean_coord;
        $position = $length_diff / $bp_length_per_unit;
        $position = $R_tip + $position;
        if ($map){
	  &ace_output($mean_coord_cds{$_}, $chrom, $position, $_, $feature);
	  next;
        }	
      }

      ##############################################################################################
      # start processing interpolated gmap of CDSes/transcript/pseudogene lying in between 2 markers
      ##############################################################################################

      if (($_ ne "NA") && $_ > $L_mean_coord && $_ < $R_mean_coord){
	
        # loop thru all marker mean coords intervals	
        for ($i=0; $i < (scalar @unique_pos_order) - 1; $i++){
	
          $gmap_down = $unique_pos_order[$i];
          $gmap_up =   $unique_pos_order[$i+1];
          $down_mean_coord = $pos_order_to_mean_coord_locus_cds{$gmap_down}->[0];
          $up_mean_coord = $pos_order_to_mean_coord_locus_cds{$gmap_up}->[0];
	
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
	
	    if ($map){
	      &ace_output($mean_coord_cds{$_}, $chrom, $position, $_, $feature);
	    }
	    last;
	  }
	}
      }
    }
  }
  print INFO "There are ", scalar @all_coords, " CDS/transcript/pseudogene (w/o counting isoforms) on $chrom:\n";
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

      # multiple loci may have a same gmap position before correction
      foreach (my $i = 0; $i < scalar @{ $pos_order_to_mean_coord_locus_cds{$gmap_down}}; $i=$i+3){

	my $locus1 = $pos_order_to_mean_coord_locus_cds{$gmap_down}->[$i+1];
	$locusf = $locus1;
	$locusf = sprintf("%10s", $locusf);
	$down_mean_coord = $pos_order_to_mean_coord_locus_cds{$gmap_down}->[$i];

	my $dnmcoord = $down_mean_coord;
	$dnmcoord = sprintf("%12.1f", $dnmcoord);
	$up_mean_coord = $pos_order_to_mean_coord_locus_cds{$gmap_up}->[$i];
	$cds = $pos_order_to_mean_coord_locus_cds{$gmap_down}->[$i+2];
	my $cdsf = $cds;
	$cdsf = sprintf("%-15s", $cdsf);

	print CMP "$chrom\t$gmdn\t$locusf\t$cdsf\t$dnmcoord\n";
	if ($down_mean_coord eq "NA") {
	  print REV "\n** Gmap marker $cds ($locus1) on $chrom has no coordinate : [$down_mean_coord] **\n";
	}

	if ($down_mean_coord ne "NA" && $up_mean_coord ne "NA" && ($down_mean_coord > $up_mean_coord)){
	  if (exists $CDS_variants{$cds}){
	    $rev_phys++;
	    $error_check = 1;
	    if ($reverse) {print REV "Reverse physical: $chrom\t$gmap_down\t$locus1\t$cds\[@{$CDS_variants{$cds}}\]\t$down_mean_coord\n"}
	  }
	  else {
	    $rev_phys++;
	    $error_check = 1;
	    if ($reverse) {print REV "Reverse physical: $chrom\t$gmap_down\t$locus1\t$cds\t$down_mean_coord\n"}
	  }
	}
      }
    }
  }
  if ($reverse){
    print REV "--->$rev_phys reverse physical(s) on Chromosome $chrom\n\n";
    my $gmap = $unique_pos_order[-1];
    my $gm = $gmap;
    $gm = sprintf("%8.4f", $gm);
    my $locus = $pos_order_to_mean_coord_locus_cds{$gmap}->[1];
    $locusf = $locus; 
    $locusf = sprintf("%10s", $locusf);
    $mean_coord = $pos_order_to_mean_coord_locus_cds{$gmap}->[0];
    my $mc = $mean_coord;
    $mc = sprintf("%12.1f", $mc);
    $cds = $pos_order_to_mean_coord_locus_cds{$gmap}->[2];
    my $cdsf = $cds;
    $cdsf = sprintf("%-15s", $cdsf);

    print CMP "$chrom\t$gm\t$locusf\t$cdsf\t$mc\n"; 
  }

  @pos_order=(); @unique_pos_order=(); @all_coords=(); %pos_order_to_mean_coord_locus_cds=();

  # counting total number of rev. phys
  $count_rev += $rev_phys;
  $rev_phys=(); @all_mean_of_each_chrom=(); @unique_all_mean_of_each_chrom=(); %all_mean_of_each_chrom=();
}

# no upload of map data to DB if in debug mode
if (!$debug){

  print "Deleting old interpolated map and uploading new ones to $database\n";

  my $tace = &tace;
  my $log = "/wormsrv2/logs/load_gmap_to_autoace"."_WS$version.$rundate.$$";

  my $command=<<END;
find elegans_CDS * where Interpolated_map_position
edit -D Interpolated_map_position

find transcript * where Interpolated_map_position
edit -D Interpolated_map_position

find pseudogene * where Interpolated_map_position
edit -D Interpolated_map_position

find locus * where Interpolated_map_position
edit -D Interpolated_map_position

pparse $acefile
save
quit
END

  open (LOAD_A,"| $tace -tsuser \"genetic_map\" $database/ >> $log") || die "Failed to upload to $database";
  print LOAD_A $command;
  close LOAD_A;

  system("chmod 777 $log");
  print "\nCheck $log to double check file uploading to $database went OK\n";
}

# Finish and exit
if($error_check == 1){
  print "\nERROR: $count_rev reverse physical(s) were found and need you to do some corrections.\n";
  system("perl5.6.1 $script_dir\/update_rev_physicals.pl -panel $count_rev -rev $revfile -comp $cmp_file -v $version");
  #system("perl5.6.1 \/nfs\/team71\/worm\/ck1\/WORMBASE_CVS\/scripts\/update_rev_physicals.pl -panel $count_rev -rev $revfile -comp $cmp_file -v $version -d");
}

else {

  my ($jah, $recipients);
  $recipients = "krb\@sanger.ac.uk, jah\@bioch.ox.ac.uk, ar2\@sanger.ac.uk, dl1\@sanger.ac.uk, pad\@sanger.ac.uk";
  $recipients = "krb\@sanger.ac.uk" if $debug;

  my @rev = `cat $revfile`;

  $jah = $revfile;
  open(JAH, ">$jah") || die $!;

  print JAH "get_interpolated_map.pl started  at $start\n","                        finished at ", &runtime,"\n\n";
  print JAH "WS$version reverse physicals / genetic map modifications\n";
  print JAH "===================================================\n\n";

  if ( -e glob("$output/rev_physical_CGC_WS$version") ) {
    my $rev_update = glob("$output/rev_physical_CGC_WS$version");
    print JAH `cat $rev_update`;
  }
  else {
    print JAH "\n", @rev;
  }
  print JAH "\n\nBelow is a list of full WS$version genetics map update\n";
  print JAH "------------------------------------------------------\n\n";
  print JAH `cat $cmp_file`;
  mail_maintainer("WS$version genetics map modification & update", $recipients, $jah);
  close JAH;
  print "\n\n* * * No reverse physicals found * * *\n";
  print "Genetics Map update has been mailed to CGC (JAH)\n\n" if !$debug;
}

exit(0);


###########################################################
#
#
#      T  H  E      S  U  B  R  O  U  T  I  N  E  S
#
#
###########################################################


##############################################################
# use latest GFF_SPLITT folder and retrieve relevant gff files
##############################################################

sub dataset {

  my ($dir, $query)= @_;
  opendir(DIR, $dir) || die "Can't read directory";
  my (@files, @vers);
  my @dir=readdir DIR;

  splice(@dir, 0,2);
  closedir (DIR);

  if ($query eq "folder") {foreach (@dir){if ($_ =~ /^WS(\d+)/){push(@vers, $1)}} return @vers}

  if ($query eq "CDS") {foreach (@dir){if ($_ =~ /^CHROMOSOME_(I|II|III|IV|V|X).CDS.gff/ ||
					    $_ =~ /^CHROMOSOME_(I|II|III|IV|V|X).clone_path.gff/){push(@files, $&)}} return @files}

  if ($query eq "rna") {foreach (@dir){if ($_ =~ /^CHROMOSOME_(I|II|III|IV|V|X).(rna|rest).gff/){push(@files, $&)}} return @files}

  if ($query eq "pseudogene") {foreach (@dir){if ($_ =~ /^CHROMOSOME_(I|II|III|IV|V|X).pseudogenes.gff/){push(@files, $&)}} return @files}
}


############################################################################################################


#########################################################
# write acefile for all CDSes with gmap/interpolated gmap
# also write Map position to locus objects
#########################################################

sub ace_output {

  my ($cds, $chrom, $gmap, $mean_coord, $feature)=@_;

  $gmap = sprintf("%8.4f", $gmap);
  my $cdsf = $cds; 
  $cdsf = sprintf ("%-15s", "$cdsf");

  my $clone;
  if ($cds !~ /.+\..+/){
    $clone = $cds;
    $clone = sprintf ("%-15s", "$clone");
  }

  my $type;

  if (exists $CDS_variants{$cds}){
    foreach(@{$CDS_variants{$cds}}){
      $type = $class{$_};

      print CL "$_ -> $type (2)" if $debug;

      my $cdsf = $_;
      $cdsf = sprintf ("%-15s", "$cdsf");

      if ($type eq "CDS"){print ACE "\nCDS : \"$_\"\n";	print CDSes "$cdsf\t";}
      if ($type eq "RNA"){print ACE "\nTranscript : \"$_\"\n"; print CDSes "$cdsf\t"}
      if ($type eq "pseudo"){print ACE "\nPseudogene : \"$_\"\n"; print PSEUDO "$cdsf\t"}
      print ACE "Interpolated_map_position\t\"$chrom\"\t$gmap\t\/\/$mean_coord (iso)\n";
      print CDSes "\t$chrom\t$gmap\t" if $type ne "pseudo";
      print PSEUDO "\t$chrom\t$gmap\t" if $type eq "pseudo";

      if (exists $sequence_linked_to_gene{$_}){
	# write interpolated_map to a Gene obj
        print ACE "\nGene : \"$Gene_info{$sequence_linked_to_gene{$_}}{'Gene'}\"\n";

	# append CGC_name to a sequence
	print CDSes "$sequence_linked_to_gene{$_}\n" if $type ne "pseudo";
	print PSEUDO "$sequence_linked_to_gene{$_}\n" if $type eq "pseudo";
        print ACE "Interpolated_map_position\t\"$chrom\"\t$gmap\t\/\/$mean_coord (iso)\n";

	# writes updated interpolated_map only to live gene ids
	if ( exists $gene_id_is_live{$Gene_info{$sequence_linked_to_gene{$_}}{'Gene'}} ){
	  print GACE "\nGene : \"$Gene_info{$sequence_linked_to_gene{$_}}{'Gene'}\"\n";
	  print GACE "Interpolated_map_position\t\"$chrom\"\t$gmap\t\/\/$mean_coord (iso)\n";
	}
	if ($comp && $locus_map{$CDS_to_gene{$_}}){
	  print MAPCOMP "$CDS_to_gene{$_}\t\"$chrom\"\tCoords:\t$gmap\tContig:\t$locus_map{$CDS_to_gene{$_}}\n";
	}
      }
      else {
	print CDSes "-\n"  if $type ne "pseudo" && $type;
	print PSEUDO "-\n" if $type eq "pseudo" && $type;
      }	
    }
  }
  else {
    if ($feature eq "CDS"){
      print ACE "\nCDS : \"$cds\"\n" if !$clone;
      print ACE "\nSequence : \"$cds\"\n" if $clone;
      print CDSes "$cdsf\t" if !$clone;
      print CLONEs "$clone\t" if $clone;
    }
    if ($feature eq "RNA"){
      print ACE "\nTranscript : \"$cds\"\n";
      print CDSes "$cdsf\t";
    }
    if ($feature eq "pseudogene"){print ACE "\nPseudogene : \"$cds\"\n"; print PSEUDO "$cdsf\t"}

    print ACE "Interpolated_map_position\t\"$chrom\"\t$gmap\t\/\/$mean_coord\n";
    print CDSes"\t$chrom\t$gmap\t" if !$clone && $feature ne "pseudogene";
    print CLONEs"\t$chrom\t$gmap\t" if $clone;
    print PSEUDO"\t$chrom\t$gmap\t" if $feature eq "pseudogene";

    if (exists $sequence_linked_to_gene{$cds}){
      # write interpolated_map to a Gene obj

      print ACE "\nGene : \"$Gene_info{$sequence_linked_to_gene{$cds}}{'Gene'}\"\n";

      # append CGC_name to a sequence
      print CDSes "$sequence_linked_to_gene{$cds}\n" if $feature ne "pseudogene";
      print PSEUDO "$sequence_linked_to_gene{$cds}\n" if $feature eq "pseudogene";
      print ACE "Interpolated_map_position\t\"$chrom\"\t$gmap\t\/\/$mean_coord\n";

      # writes updated interpolated_map only to live gene ids
      if ( exists $gene_id_is_live{$Gene_info{$sequence_linked_to_gene{$cds}}{'Gene'}} ){
	print GACE "\nGene : \"$Gene_info{$sequence_linked_to_gene{$cds}}{'Gene'}\"\n";
	print GACE "Interpolated_map_position\t\"$chrom\"\t$gmap\t\/\/$mean_coord\n";
      }
      if ($comp && exists $locus_map{$CDS_to_gene{$cds}}){
	print MAPCOMP "$CDS_to_gene{$cds}\t\"$chrom\"\tCoords:\t$gmap\tContig:\t$locus_map{$CDS_to_gene{$cds}}\n";
      }
    }
    else {
      print CDSes "-\n" if !$clone && $feature ne "pseudogene";
      print CLONEs "\n" if $clone;
      print PSEUDO "-\n" if $feature eq "pseudogene";
    }
  }
}


###########################################################################################################################


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


B <-d / -debug):

            Interpolated_map_positions will not be uploaded to autoace.
