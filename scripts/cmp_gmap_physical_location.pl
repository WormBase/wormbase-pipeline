#!/usr/local/bin/perl5.6.1 -w

# cmp_genetics_seq_location.pl

# by Chao-Kung Chen [030113]

# Last updated on: $Date: 2003-03-03 12:42:12 $
# Last updated by: $Author: ck1 $

use strict;                    
use lib "/wormsrv2/scripts/";
use Wormbase;
use Cwd 'chdir';

#################################################
# look for latest WS release in GFF_SPLITS folder
#################################################

my $gff_dir = "/wormsrv2/autoace/GFF_SPLITS/";

my @versions=dataset($gff_dir, "folder");

my @order = sort {$a <=> $b} @versions;
#print @order, "\n";

my $gff_location = "/wormsrv2/autoace/GFF_SPLITS/WS"."$order[-1]";
chdir $gff_location;

##########################################################
# retrieve data from gff file: CHROMOSOME_number.genes.gff
##########################################################

my @gff_files=dataset($gff_location, "files");
my (@data, %CDS_mapping);

foreach (@gff_files){
  #print $_, ">>>\n";
  @data = `less $_ | cut -f 1,4,5,9`;
  foreach (@data){
    my ($chrom, $left, $right, $junk1, $CDS, $junk2)= split(/\s/,$_); 
    $CDS =~ s/\"//g;
    $chrom =~ s/CHROMOSOME_//;
    push(@{$CDS_mapping{$CDS}}, $chrom, $left, $right);
  }
}

#print keys %CDS_mapping, "\n";
#foreach (sort keys %CDS_mapping){
#  print "$_ -> ${@{$CDS_mapping{$_}}}[0]\n";
#}

############################################################################
# retrieve locus data: chromosome, map position, genomic_sequence/transcript 
############################################################################

my $genetics_mapping=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/genetics_mapping.def" quit
EOF
my $get_gmap_positions_of_all_chroms=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/get_gmap_positions_of_all_chroms.def" quit
EOF


my $ga_dir="/wormsrv1/geneace";
my (%genetics_mapping, %chrom_pos);
  
open (FH, "echo '$genetics_mapping' | tace $ga_dir | ") || die "Couldn't access Geneace\n";
while (<FH>){
  chomp($_);
  if ($_ =~ /^\"(.+)\"\s+\"(.+)\"\s+\"(.+)\"\s+(-\d+\.\d+|\d+\.\d+)/){   # postive and negative values  
    #print "$1 -> $2 -> $3 -> $4\n";
    push(@{$genetics_mapping{$3}}, $2, $1, $4);   # $3=cds, $2=chrom, $1=locus, $4=gmap position
  }
}
close FH;

open (FH, "echo '$get_gmap_positions_of_all_chroms' | tace $ga_dir | ") || die "Couldn't access Geneace\n";
while (<FH>){
  chomp($_);
  if ($_ =~ /^\"(.+)\"\s+(-\d+\.\d+|\d+\.\d+)/){   # postive and negative values  
    push(@{$chrom_pos{$1}}, $2);
  }
}
close FH;

############################################
# calculate units of each chromosome on gmap
############################################

my (%unit_size, %chrom_left_end);


foreach (sort keys %chrom_pos){
  @order=sort {$a <=>$b} @{$chrom_pos{$_}};
  my $units = ($order[-1] - $order[0]); 
  print "Chromosome $_:\tLEFT end\t$order[0] - RIGHT end\t$order[-1]\t-> Total: $units Units\n";
  $chrom_left_end{$_}=$order[0];
  $unit_size{$_}=$units;
} 	


#foreach (keys %unit_size){             # physical length / unit
 # print "$_ -> $unit_size{$_}\n";
#}	

#####################################################################
# hard coded chrom length - 
# doesn't make much sense to run script for a "completed" genome, Or?
#####################################################################

my $chrom_I   = 15080471;   # mean ~ 5714127.372384
my $chrom_II  = 15279300;   # mean ~ 5701136.79708
my $chrom_III = 13783267;   # mean ~ 747991.53714
my $chrom_IV  = 17493790;   # mean ~ 11026176.85094
my $chrom_V   = 20922239;   # mean ~ 9242251.38708
my $chrom_X   = 17705013;   # mean ~ 7685029.219064
#my $MtDNA     = 13794;

#####################################################################################
# divide physical length of each chromosome by gmap units of corresponding chromosome
#####################################################################################

my $chr_I_leng_unit   = $chrom_I / $unit_size{'I'};
my $chr_II_leng_unit  = $chrom_II / $unit_size{'II'}; 
my $chr_III_leng_unit = $chrom_III / $unit_size{'III'};
my $chr_IV_leng_unit  = $chrom_IV / $unit_size{'IV'};
my $chr_V_leng_unit   = $chrom_V / $unit_size{'V'};
my $chr_X_leng_unit  = $chrom_X / $unit_size{'X'};

print "\n\n";
print "Chrom_I:   Length  $chrom_I\t-> $chr_I_leng_unit\tKb/U\n";
print "Chrom_II:  Length  $chrom_II\t-> $chr_II_leng_unit\tKb/U\n";
print "Chrom_III: Length  $chrom_III\t-> $chr_III_leng_unit\tKb/U\n";
print "Chrom_IV:  Length  $chrom_IV\t-> $chr_IV_leng_unit\tKb/U\n";
print "Chrom_V:   Length  $chrom_V\t-> $chr_V_leng_unit\tKb/U\n";
print "Chrom_X:   Length  $chrom_X\t-> $chr_X_leng_unit\tKb/U\n";

################################################################
# calculate mean value of left and right coordinates of each CDS
# and assign gmap position by physical mean 
################################################################

my ($mean, %mean_gmap_CH, %mean_gmap_GA, %chrom_plot_CH, %chrom_plot_GA);

system("rm -f /wormsrv1/chaokung/DOCS/gmap_ch");
system("rm -f /wormsrv1/chaokung/DOCS/gmap_ga");

open(GA, ">>/wormsrv1/chaokung/DOCS/gmap_ga") || die "Can't write to file!";
open(CH, ">>/wormsrv1/chaokung/DOCS/gmap_ch") || die "Can't write to file!";

#############################################################
# hard coded mean of each -0 "centromere"  maker gene on gmap
#############################################################

my ($chrom, $unit);
my (%mean_gmap_CH_I, %mean_gmap_CH_II, %mean_gmap_CH_III, %mean_gmap_CH_IV, %mean_gmap_CH_V, %mean_gmap_CH_X);
my (%mean_gmap_GA_I, %mean_gmap_GA_II, %mean_gmap_GA_III, %mean_gmap_GA_IV, %mean_gmap_GA_V, %mean_gmap_GA_X);

foreach (sort keys %genetics_mapping){
 # print "$_ -> @{$genetics_mapping{$_}}\n";
  if (exists ${@{$genetics_mapping{$_}}}[0] && exists ${@{$CDS_mapping{$_}}}[0] ){
   $mean = (${@{$CDS_mapping{$_}}}[1] + ${@{$CDS_mapping{$_}}}[2]) / 2;
   #create_gmap_position($mean, ${@{$CDS_mapping{$_}}}[0]);
   $chrom =  ${@{$CDS_mapping{$_}}}[0];
   if ($chrom eq "I")  {$unit = $chr_I_leng_unit; push(@{$mean_gmap_CH_I{$mean}}, ($mean / $unit), ${@{$genetics_mapping{$_}}}[1])}
   if ($chrom eq "II") {$unit = $chr_II_leng_unit; push(@{$mean_gmap_CH_II{$mean}}, ($mean / $unit), ${@{$genetics_mapping{$_}}}[1])}
   if ($chrom eq "III"){$unit = $chr_III_leng_unit; push(@{$mean_gmap_CH_III{$mean}}, ($mean / $unit), ${@{$genetics_mapping{$_}}}[1])} 
   if ($chrom eq "IV") {$unit = $chr_IV_leng_unit; push(@{$mean_gmap_CH_IV{$mean}}, ($mean / $unit), ${@{$genetics_mapping{$_}}}[1])}
   if ($chrom eq "V")  {$unit = $chr_V_leng_unit; push(@{$mean_gmap_CH_V{$mean}}, ($mean / $unit), ${@{$genetics_mapping{$_}}}[1])}
   if ($chrom eq "X")  {$unit = $chr_X_leng_unit; push(@{$mean_gmap_CH_X{$mean}}, ($mean / $unit), ${@{$genetics_mapping{$_}}}[1])}

   get_gmap_position($mean, ${@{$genetics_mapping{$_}}}[0], ${@{$genetics_mapping{$_}}}[2]);
 }
}

my @means;
@means=keys %mean_gmap_CH_I;
@means=sort {$a<=>$b} @means;
foreach (@means){
  print CH "I\t${@{$mean_gmap_CH_I{$_}}}[1]\t$_\t${@{$mean_gmap_CH_I{$_}}}[0]\t$mean_gmap_GA_I{$_}\n";
}
@means=keys %mean_gmap_CH_II;  
@means=sort {$a<=>$b} @means;
foreach (@means){
  print CH "II\t${@{$mean_gmap_CH_II{$_}}}[1]\t$_\t${@{$mean_gmap_CH_II{$_}}}[0]\t$mean_gmap_GA_II{$_}\n";
}
@means=keys %mean_gmap_CH_III;  
@means=sort {$a<=>$b} @means;
foreach (@means){
  print CH "III\t${@{$mean_gmap_CH_III{$_}}}[1]\t$_\t${@{$mean_gmap_CH_III{$_}}}[0]\t$mean_gmap_GA_III{$_}\n";
}
@means=keys %mean_gmap_CH_IV;  
@means=sort {$a<=>$b} @means;
foreach (@means){
  print CH "IV\t${@{$mean_gmap_CH_IV{$_}}}[1]\t$_\t${@{$mean_gmap_CH_IV{$_}}}[0]\t$mean_gmap_GA_IV{$_}\n";
}
@means=keys %mean_gmap_CH_V;  
@means=sort {$a<=>$b} @means;
foreach (@means){
  print CH "V\t${@{$mean_gmap_CH_V{$_}}}[1]\t$_\t${@{$mean_gmap_CH_V{$_}}}[0]\t$mean_gmap_GA_V{$_}\n";
}
@means=keys %mean_gmap_CH_X;  
@means=sort {$a<=>$b} @means;
foreach (@means){
  print CH "X\t${@{$mean_gmap_CH_X{$_}}}[1]\t$_\t${@{$mean_gmap_CH_X{$_}}}[0]\t$mean_gmap_GA_X{$_}\n";
}

#################

sub get_gmap_position {
  my ($mean, $chrom, $position) = @_;
  if ($position  <= 0){
    #my $positive_pos = -$position; 
    #print $position;
    if ($chrom eq "I") {$mean_gmap_GA_I{$mean}     = -$chrom_left_end{$chrom} + $position}
    if ($chrom eq "II") {$mean_gmap_GA_II{$mean}   = -$chrom_left_end{$chrom} + $position}   
    if ($chrom eq "III") {$mean_gmap_GA_III{$mean} = -$chrom_left_end{$chrom} + $position}   
    if ($chrom eq "IV") {$mean_gmap_GA_IV{$mean}   = -$chrom_left_end{$chrom} + $position}
    if ($chrom eq "V")  {$mean_gmap_GA_V{$mean}    = -$chrom_left_end{$chrom} + $position}
    if ($chrom eq "X")  {$mean_gmap_GA_X{$mean}    = -$chrom_left_end{$chrom} + $position}
  }
  else {
    if ($chrom eq "I")  {$position += -$chrom_left_end{$chrom}; $mean_gmap_GA_I{$mean}   = $position}
    if ($chrom eq "II") {$position += -$chrom_left_end{$chrom}; $mean_gmap_GA_II{$mean}  = $position}
    if ($chrom eq "III"){$position += -$chrom_left_end{$chrom}; $mean_gmap_GA_III{$mean} = $position} 
    if ($chrom eq "IV") {$position += -$chrom_left_end{$chrom}; $mean_gmap_GA_IV{$mean}  = $position}  
    if ($chrom eq "V")  {$position += -$chrom_left_end{$chrom}; $mean_gmap_GA_V{$mean}   = $position}
    if ($chrom eq "X")  {$position += -$chrom_left_end{$chrom}; $mean_gmap_GA_X{$mean}   = $position}
  }
}  

##################################################################
# output discrepancy of chrom. mapping by genetics and coordinates
##################################################################

foreach (sort keys %genetics_mapping){
  if (exists ${@{$genetics_mapping{$_}}}[0] && exists ${@{$CDS_mapping{$_}}}[0] ){
    if (${@{$genetics_mapping{$_}}}[0] ne ${@{$CDS_mapping{$_}}}[0]){
      print "$_ is mapped to ${@{$genetics_mapping{$_}}}[0] by genetics\n";
      print "$_ is mapped to ${@{$CDS_mapping{$_}}}[0] by coordinates.\n\n"
    }
  }
}

#############
# subroutines
#############

sub dataset {
  my ($dir, $query)= @_;
#  print $dir, ">>\n";
  opendir(DIR, $dir) || die "Can't read directory";
  my (@files, @vers);
  my @dir=readdir DIR;

  splice(@dir, 0,2);
  closedir (DIR);

  if ($query eq "folder"){
    foreach (@dir){
 #     print $_, "##\n";
      if ($_ =~ /^WS(\d+)/){
	push(@vers, $1);
      }  
    }
    return @vers;
  }
  if ($query eq "files"){
    foreach (@dir){
      if ($_ =~ /^CHROMOSOME_([A-Z]+).genes.gff/){
	#print $&, "\n";
	push(@files, $&);
      }  
    }
  }
  return @files;
}


__END__
