#!/usr/local/bin/perl5.6.1 -w

# get_interpolated_gmap.pl

# by Chao-Kung Chen [030113]

# This script calculates interpolated gmap for CDS/transcripts lying between as well as outside genetic markers.
# Output ace file of such information

# Last updated on: $Date: 2003-03-14 18:42:25 $
# Last updated by: $Author: ck1 $

use strict;                    
use lib "/wormsrv2/scripts/";
use Wormbase;
use Cwd 'chdir';
use Getopt::Long;

###################################################
# variables and command-line options with aliases # 
###################################################

my ($diff, $reverse);

GetOptions ("diff" => \$diff,
	    "reverse" => \$reverse,
	   );

#################################################
# look for latest WS release in GFF_SPLITS folder
#################################################

my $gff_dir = "/wormsrv2/autoace/GFF_SPLITS/";

my $ga_dir="/wormsrv1/geneace";
#my $ga_dir="/wormsrv1/chaokung/wormdata/CK1_GENEACE/";
my @versions=dataset($gff_dir, "folder");

my @order = sort {$a <=> $b} @versions;

my $gff_location = "/wormsrv2/autoace/GFF_SPLITS/WS"."$order[-1]";

my $file = glob "/wormsrv1/chaokung/DOCS/cmp_gmap_with_coord_order";
if($file){system("rm -f $file")}

open(OUT, ">>/wormsrv1/chaokung/DOCS/cmp_gmap_with_coord_order") || die $!;
print OUT "\nGFF file version from $gff_location\n\n";

chdir $gff_location;

my (@data, %CDS_mapping, %CDS_isoforms_mapping, %CDS_variants, %RNA_mapping);
my (%chrom_pos, %genetics_mapping, %I_gmap_cds_locus, %II_gmap_cds_locus, %III_gmap_cds_locus, 
    %IV_gmap_cds_locus, %V_gmap_cds_locus, %X_gmap_cds_locus); 
my ($cds, $parts, @coords, $i, $mean_coords, %cds_mean, %mean_coord_cds,
    %mean_coord_cds_I,  %mean_coord_cds_II,  %mean_coord_cds_III,  
    %mean_coord_cds_IV,  %mean_coord_cds_V,  %mean_coord_cds_X,
   );

if ($diff){&gff_gmap; exit (0)} else {&gff_gmap}

if ($reverse){
  open (LOG, ">/wormsrv2/logs/reverse_physicals") || die $!;
  print LOG "Checking linearity of gmap markers (reverse physicals) ....\n";
  print LOG "\n";
}

sub gff_gmap {

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
      if ($_ =~ /^CHR.+/){
	my ($chrom, $left, $right, $junk1, $CDS, $junk2)= split(/\s/,$_); 
     	$CDS =~ s/\"//g;
	my $ori = $CDS;
	$cds_count++;
	push(@{$CDS_isoforms_mapping{$CDS}}, $chrom, $left, $right);
	if ($CDS =~ /(.+\.\d+)\D+/){
	  $CDS = $1;
	  push(@{$CDS_variants{$CDS}}, $ori); #parent of isoforms
	}
	else {$CDS = $CDS}
	$chrom =~ s/CHROMOSOME_//;
	push(@{$CDS_mapping{$CDS}}, $chrom, $left, $right);
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
	$rna_count++;
	push(@{$CDS_isoforms_mapping{$RNA}}, $chrom, $left, $right);
	if ($RNA =~ /(.+\.\d+)\D+/){$RNA = $1}
	else {$RNA = $RNA}
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
    $parts = scalar @{$CDS_mapping{$_}};

    ####################
    # cds has no isoform
    ####################

    if ($parts == 3){  
      # print "$_ -> ${@{$CDS_mapping{$_}}}[1] --- ${@{$CDS_mapping{$_}}}[2]\n";
      $mean_coords = (${@{$CDS_mapping{$_}}}[1] + ${@{$CDS_mapping{$_}}}[2]) / 2;  
      $cds_mean{$_}= $mean_coords;
      if (${@{$CDS_mapping{$_}}}[0] eq "I"){$mean_coord_cds_I{$mean_coords}=$_}
      if (${@{$CDS_mapping{$_}}}[0] eq "II"){$mean_coord_cds_II{$mean_coords}=$_}
      if (${@{$CDS_mapping{$_}}}[0] eq "III"){$mean_coord_cds_III{$mean_coords}=$_}
      if (${@{$CDS_mapping{$_}}}[0] eq "IV"){$mean_coord_cds_IV{$mean_coords}=$_}
      if (${@{$CDS_mapping{$_}}}[0] eq "V"){$mean_coord_cds_V{$mean_coords}=$_}
      if (${@{$CDS_mapping{$_}}}[0] eq "X"){$mean_coord_cds_X{$mean_coords}=$_}
      }
    else {
    
    ##################
    # cds has isoforms
    ##################
       
      #  print $parts, "-> $_\n";
      for ($i=0; $i < $parts; $i=$i+3 ){
	 push(@coords, ${@{$CDS_mapping{$_}}}[$i+1], ${@{$CDS_mapping{$_}}}[$i+2]);  # get coords of all isoforms
       }
       #   print "$_ -> @coords\n";
       @coords = sort {$a <=> $b} @coords;
       $mean_coords = ($coords[-1] + $coords[0]) / 2;
       $cds_mean{$_}= $mean_coords;
       if (${@{$CDS_mapping{$_}}}[0] eq "I"){$mean_coord_cds_I{$mean_coords}=$_}
       if (${@{$CDS_mapping{$_}}}[0] eq "II"){$mean_coord_cds_II{$mean_coords}=$_}
       if (${@{$CDS_mapping{$_}}}[0] eq "III"){$mean_coord_cds_III{$mean_coords}=$_}
       if (${@{$CDS_mapping{$_}}}[0] eq "IV"){$mean_coord_cds_IV{$mean_coords}=$_}
       if (${@{$CDS_mapping{$_}}}[0] eq "V"){$mean_coord_cds_V{$mean_coords}=$_}
       if (${@{$CDS_mapping{$_}}}[0] eq "X"){$mean_coord_cds_X{$mean_coords}=$_}
       @coords =();
    }
  } 

  ######################################################################################
  # retrieve locus data: chromosome, map position, genomic_sequence/transcript from gmap
  ######################################################################################

  my $marker_gmap_of_each_chrom=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/marker_gmap_of_each_chrom.def" quit
EOF

  open (FH, "echo '$marker_gmap_of_each_chrom' | tace $ga_dir | ") || die "Couldn't access Geneace\n";
  while (<FH>){
    chomp($_);
    # $4=cds/transcript, $2=chrom, $1=locus, $3=gmap position
    if ($_ =~ /^\"(.+)\"\s+\"(.+)\"\s+(-\d+\.\d+|\d+\.\d+)\s+\"(.+)\"/){   # postive and negative values

      if ($2 eq "I"){push (@{$I_gmap_cds_locus{$3}}, $4, $1)}
      if ($2 eq "II"){push (@{$II_gmap_cds_locus{$3}}, $4, $1)}
      if ($2 eq "III"){push (@{$III_gmap_cds_locus{$3}}, $4, $1)}
      if ($2 eq "IV"){push (@{$IV_gmap_cds_locus{$3}}, $4, $1)}
      if ($2 eq "V"){push (@{$V_gmap_cds_locus{$3}}, $4, $1)}
      if ($2 eq "X"){push (@{$X_gmap_cds_locus{$3}}, $4, $1)}
      push(@{$genetics_mapping{$4}}, $2, $1, $3);  
      push(@{$chrom_pos{$2}}, $3);
    }
  }
  close FH;

  ##################################################################
  # output discrepancy of chrom. mapping by genetics and coordinates
  # this is called from geneace_check.pl 
  ##################################################################

  $count = 0;
  if ($diff){
    open(LOG, ">/wormsrv2/logs/mapping_diff") || die $!;
    foreach (sort keys %genetics_mapping){
      if (exists ${@{$genetics_mapping{$_}}}[0] && exists ${@{$CDS_isoforms_mapping{$_}}}[0] ){
        ${@{$CDS_isoforms_mapping{$_}}}[0] =~ s/CHROMOSOME_//;
        if (${@{$genetics_mapping{$_}}}[0] ne ${@{$CDS_isoforms_mapping{$_}}}[0]){
          $count++;
          print LOG "$_ (${@{$genetics_mapping{$_}}}[1]) is mapped to ${@{$genetics_mapping{$_}}}[0] by genetics\n";
          print LOG "$_ (${@{$genetics_mapping{$_}}}[1]) is mapped to ${@{$CDS_isoforms_mapping{$_}}}[0] by coordinates.\n\n"
        }
      }
    }
    print LOG "ERROR: There are $count discrepancies in genetics/chrom. coords mapping\n"; 
  }
  
}

#########################################################
# get each gmap position, its mean of coord on each chrom
# the variables are stored in @chrom_I/II/III/IV/V/X
#########################################################

my (@pos_order, @chrom, @chrom_I, @chrom_II, @chrom_III, @chrom_IV, @chrom_V, @chrom_X, $marker_count, $total); 
my $outside=0;

foreach (sort keys %chrom_pos){
 my $marker_count = 0;
  if ($_ eq "I"){
    @pos_order = sort {$a <=>$b} @{$chrom_pos{$_}};
    get_gmap_mean_of_each_chrom(\@pos_order, \%I_gmap_cds_locus, \%mean_coord_cds_I, "I", @chrom);
  }
  if ($_ eq "II"){
    @pos_order = sort {$a <=>$b} @{$chrom_pos{$_}};
    get_gmap_mean_of_each_chrom(\@pos_order, \%II_gmap_cds_locus, \%mean_coord_cds_II, "II", @chrom);
  }
  if ($_ eq "III"){
    @pos_order = sort {$a <=>$b} @{$chrom_pos{$_}};
    get_gmap_mean_of_each_chrom(\@pos_order, \%III_gmap_cds_locus,  \%mean_coord_cds_III, "III", @chrom);
  }
  if ($_ eq "IV"){
    @pos_order = sort {$a <=>$b} @{$chrom_pos{$_}};
    get_gmap_mean_of_each_chrom(\@pos_order, \%IV_gmap_cds_locus,  \%mean_coord_cds_IV, "IV", @chrom);
  }
  if ($_ eq "V"){
    @pos_order = sort {$a <=>$b} @{$chrom_pos{$_}};
    get_gmap_mean_of_each_chrom(\@pos_order, \%V_gmap_cds_locus,  \%mean_coord_cds_V, "V", @chrom);
  }
  if ($_ eq "X"){
    @pos_order = sort {$a <=>$b} @{$chrom_pos{$_}};
    ($outside, $total)=get_gmap_mean_of_each_chrom(\@pos_order, \%X_gmap_cds_locus,  \%mean_coord_cds_X, "X", @chrom);
  }
}

# $outside -> total cds/transcripts in genome outside of tip markers
my $percent=$outside/$total *100;
$percent = sprintf("%.2f", $percent);

print OUT "\nSum of predicted cds/transcripts (not counting isoforms) outside of tip(L) or tip(R) markers: $outside ---> $percent% in genome\n";

#############
# subroutines
#############

sub dataset {
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

sub get_gmap_mean_of_each_chrom {

  my ($pos_order_chrom, $chrom_gmap_cds_locus, $mean_coord_cds, $chromosom, @chrom) = @_;
  my %Chrom_gmap_cds_locus = %$chrom_gmap_cds_locus;
  my %Mean_coord_cds = %$mean_coord_cds;
  my ($e, @unique_order_chrom);
  my %seen=();

  foreach my $e (@$pos_order_chrom){
    push(@unique_order_chrom, $e) unless $seen{$e}++;
  }

  foreach my $each (@unique_order_chrom){ # do not sort
    #print $each, ">>\n";
    $marker_count++;
    if (exists $Chrom_gmap_cds_locus{$each}){ 
      push(@chrom, $each);   # gmap position
      print OUT "\n$chromosom\t$each\t${@{$Chrom_gmap_cds_locus{$each}}}[1]\t";  
      if (${@{$Chrom_gmap_cds_locus{$each}}}[0] =~ /(^.+\.\d+)\w/){ # cds
        #print "$cds_mean{$1}\n";   
        push(@chrom, $cds_mean{$1});  # mean coord
	print OUT "$cds_mean{$1}\t";
      }
      else {
        #print "${@{$I_gmap_cds_locus{$each}}}[0] $cds_mean{${@{$I_gmap_cds_locus{$each}}}[0]} "; 
        #print "$cds_mean{${@{$I_gmap_cds_locus{$each}}}[0]}\n"; 
        push(@chrom, $cds_mean{${@{$Chrom_gmap_cds_locus{$each}}}[0]});  # mean coord
        print OUT "$cds_mean{${@{$Chrom_gmap_cds_locus{$each}}}[0]}\t";
      }
    }
  }

  print "$marker_count markers on Chrom $_\n";         
  $marker_count = 0;
  #print "@chrom\n";
  my ($i, %gff_gmap, $mean_up, $mean_down, $gmap_up, $gmap_down, $baseline, $position, $diff_bp);
  my @num_cds=keys %Mean_coord_cds;
  my $out; # counting cds/transcripts outside of tip markers for each chromosom
  my $reverse_physical=0;

  if ($reverse){
     for($i=0; $i < (scalar @chrom)-3; $i=$i+2){           # gmap
      #print $i," ->", $i+1," -> ",$i+2," -> ",$i+3,"\n";
      $mean_down =  $chrom[$i+1]; 
      $mean_up = $chrom[$i+3];
      $gmap_up=$chrom[$i+2];
      $gmap_down=$chrom[$i];
      #print ">$mean_down -> $mean_up\n";
      if ($mean_down > $mean_up){
	$reverse_physical++;
	if(exists $CDS_variants{$Mean_coord_cds{$mean_up}}){
	  print LOG "Reverse physicals found: $Mean_coord_cds{$mean_up} (@{$CDS_variants{$Mean_coord_cds{$mean_up}}}) ($mean_up) on ${@{$CDS_mapping{$Mean_coord_cds{$mean_up}}}}[0]\n"; 
	}
	else {
	 print LOG "Reverse physicals found: $Mean_coord_cds{$mean_up} ($mean_up) on ${@{$CDS_mapping{$Mean_coord_cds{$mean_up}}}}[0]\n";
	}
      }  
    }
    print LOG "\nTotal reverse physicals on ${@{$CDS_mapping{$Mean_coord_cds{$mean_up}}}}[0]: $reverse_physical\n\n";
  }

  foreach (sort keys %Mean_coord_cds){                    # gff w/o isoforms
    if ($_ < $chrom[1] || $_ > $chrom[-1]){$out++; $outside++}    # count cds/transcripts outside of tip markers
    for($i=0; $i < (scalar @chrom)-3; $i=$i+2){           # gmap
      #print $i," ->", $i+1," -> ",$i+2," -> ",$i+3,"\n";
      $mean_down =  $chrom[$i+1]; 
      $mean_up = $chrom[$i+3];
      $gmap_up=$chrom[$i+2];
      $gmap_down=$chrom[$i];
    
      if (($_>$mean_down) && ($_<$mean_up) ){  

        ################################        
        # negative gmap VS negative gmap
        ################################

        if ($gmap_down < 0 && $gmap_up < 0){   
          $baseline=($mean_up - $mean_down) / (-$gmap_down + $gmap_up);
        }
        ################################    
        # positive gmap VS positive gmap
        ################################

        if ($gmap_down > 0 && $gmap_up > 0){   
          $baseline=($mean_up - $mean_down) / ($gmap_up - $gmap_down);
        }
        ################################
        # negative gmap VS positive gmap  
        ################################

        if ($gmap_down < 0 && $gmap_up > 0){ 
          $baseline=($mean_up - $mean_down) / ($gmap_up - $gmap_down);
        }

        $diff_bp=$_- $mean_down; 
        $position=$diff_bp/$baseline; 
        $position=$gmap_down+$position;
        #print "\n";
        #print $mean_down," ->", $mean_up,">1\n";
        #print $gmap_down," -> ",$gmap_up,">2\n";
        if(exists $CDS_variants{$Mean_coord_cds{$_}}){
	  print "### $Mean_coord_cds{$_} (@{$CDS_variants{$Mean_coord_cds{$_}}}) ($_) on ${@{$CDS_mapping{$Mean_coord_cds{$_}}}}[0] -> interpolated: $position ###\n";
	}  
	else {
	  print "### $Mean_coord_cds{$_} ($_) on ${@{$CDS_mapping{$Mean_coord_cds{$_}}}}[0] -> interpolated: $position ###\n";
	}
      }
    }
  }
  if (!$reverse){
    my @total_cds=keys %CDS_mapping;
    my @total_rna=keys %RNA_mapping;
    print OUT "\n\ntotal cds/transcripts w\/o isoforms: ", (scalar @total_cds) + (scalar @total_rna), "\n";
    print OUT "total cds/transcripts on chromosom $chromosom: ", scalar @num_cds, "\n";
    print OUT "On Chromosome $chromosom: $out cds are outside of tip(L) or tip(R) markers\n";
    if ($chromosom eq "X"){return $outside, (scalar @total_cds) + (scalar @total_rna)} 
  }
  @chrom=(); 
}

     
__END__
