#!/usr/local/bin/perl5.6.1 -w

# cmp_genetics_seq_location.pl

# by Chao-Kung Chen [030113]

# Last updated on: $Date: 2003-03-07 18:21:46 $
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

my $gff_location = "/wormsrv2/autoace/GFF_SPLITS/WS"."$order[-1]";
chdir $gff_location;

##########################################################
# retrieve data from gff file: CHROMOSOME_number.genes.gff
##########################################################

my @gff_files_cds=dataset($gff_location, "genes");
my @gff_files_rna=dataset($gff_location, "rna");

my (@data, %CDS_mapping, %RNA_mapping);
my $cds_count=0;
my $rna_count=0;

foreach (@gff_files_cds){
  @data = `less $_ | cut -f 1,4,5,9`;
  foreach (@data){
    my ($chrom, $left, $right, $junk1, $CDS, $junk2)= split(/\s/,$_); 
    $CDS =~ s/\"//g;
    $cds_count++;
    if ($CDS =~ /(^.+\.\d+)\w/){$CDS = $1}
    else {$CDS = $CDS}
    $chrom =~ s/CHROMOSOME_//;
    push(@{$CDS_mapping{$CDS}}, $chrom, $left, $right);
  }
}
print "all cds(with isoforms): ",$cds_count, "\n";

foreach (@gff_files_rna){
  @data = `less $_ | cut -f 1,4,5,9,10`;
  foreach (@data){
    my ($chrom, $left, $right, $junk, $RNA)= split(/\s/,$_); 
    $RNA =~ s/\"//g;
    $rna_count++;
    $chrom =~ s/CHROMOSOME_//;
    push(@{$RNA_mapping{$RNA}}, $chrom, $left, $right);
  }
}

print "all transcripts: ",$rna_count, "\n";

##############################################################
# from GFF file
# get mean value of chrom coords of each CDS/Transcript, 
# mean of isoforms is the mean of the lowest and highest value
##############################################################

my ($cds, $parts, @coords, $i, $mean_coords, %cds_mean, %mean_coord_cds,
    %mean_coord_cds_I,  %mean_coord_cds_II,  %mean_coord_cds_III,  
    %mean_coord_cds_IV,  %mean_coord_cds_V,  %mean_coord_cds_X,
    );
my $count=0;

foreach (sort keys %CDS_mapping){
  $parts = scalar @{$CDS_mapping{$_}};
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
  #  print $parts, "-> $_\n";
    for ($i=0; $i < $parts; $i=$i+3 ){
      $count++;
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

$count=0;

foreach (sort keys %RNA_mapping){
  $parts = scalar @{$RNA_mapping{$_}};
  if ($parts == 3){	
   # print "$_ -> ${@{$RNA_mapping{$_}}}[1] --- ${@{$RNA_mapping{$_}}}[2]\n";
    $mean_coords = (${@{$RNA_mapping{$_}}}[1] + ${@{$RNA_mapping{$_}}}[2]) / 2;
    $cds_mean{$_}= $mean_coords;
    if (${@{$RNA_mapping{$_}}}[0] eq "I"){$mean_coord_cds_I{$mean_coords}=$_}
    if (${@{$RNA_mapping{$_}}}[0] eq "II"){$mean_coord_cds_II{$mean_coords}=$_}
    if (${@{$RNA_mapping{$_}}}[0] eq "III"){$mean_coord_cds_III{$mean_coords}=$_}
    if (${@{$RNA_mapping{$_}}}[0] eq "IV"){$mean_coord_cds_IV{$mean_coords}=$_}
    if (${@{$RNA_mapping{$_}}}[0] eq "V"){$mean_coord_cds_V{$mean_coords}=$_}
    if (${@{$RNA_mapping{$_}}}[0] eq "X"){$mean_coord_cds_X{$mean_coords}=$_}
  }
  else {
  #  print $parts, "-> $_\n";
    for ($i=0; $i < $parts; $i=$i+3 ){
      $count++;
      push(@coords, ${@{$RNA_mapping{$_}}}[$i+1], ${@{$RNA_mapping{$_}}}[$i+2]);  # get coords of all isoforms
    }
 #   print "$_ -> @coords\n";
    @coords = sort {$a <=> $b} @coords;
    $mean_coords = ($coords[-1] + $coords[0]) / 2;
    $cds_mean{$_}= $mean_coords;
    if (${@{$RNA_mapping{$_}}}[0] eq "I"){$mean_coord_cds_I{$mean_coords}=$_}
    if (${@{$RNA_mapping{$_}}}[0] eq "II"){$mean_coord_cds_II{$mean_coords}=$_}
    if (${@{$RNA_mapping{$_}}}[0] eq "III"){$mean_coord_cds_III{$mean_coords}=$_}
    if (${@{$RNA_mapping{$_}}}[0] eq "IV"){$mean_coord_cds_IV{$mean_coords}=$_}
    if (${@{$RNA_mapping{$_}}}[0] eq "V"){$mean_coord_cds_V{$mean_coords}=$_}
    if (${@{$RNA_mapping{$_}}}[0] eq "X"){$mean_coord_cds_X{$mean_coords}=$_}
    @coords =();
  }
} 	 

###########
# test
###########
=start
foreach (keys %cds_mean){
  if ($_ eq "F47G6.4"){
    print "$_ -> $cds_mean{$_}\n";
  }
}
foreach (keys %mean_coord_cds_X){
  print "$_ -> $mean_coord_cds_X{$_}\n";
}

=end
=cut


######################################################################################
# retrieve locus data: chromosome, map position, genomic_sequence/transcript from gmap
######################################################################################

my $marker_gmap_of_each_chrom=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/marker_gmap_of_each_chrom.def" quit
EOF

my $ga_dir="/wormsrv1/geneace";
my (%genetics_mapping, %chrom_pos);
my (%I_gmap_cds_locus, %II_gmap_cds_locus, %III_gmap_cds_locus, %IV_gmap_cds_locus, %V_gmap_cds_locus, %X_gmap_cds_locus); 

 
open (FH, "echo '$marker_gmap_of_each_chrom' | tace $ga_dir | ") || die "Couldn't access Geneace\n";
while (<FH>){
  chomp($_);
  # $4=cds/transcript, $2=chrom, $1=locus, $3=gmap position
  if ($_ =~ /^\"(.+)\"\s+\"(.+)\"\s+(-\d+\.\d+|\d+\.\d+)\s+\"(.+)\"/){   # postive and negative values
 #  print "$1 -> $2 -> $3 -> $4\n";
    $cds = $4;
    if ($2 eq "I"){push (@{$I_gmap_cds_locus{$3}}, $cds, $1)}
    if ($2 eq "II"){push (@{$II_gmap_cds_locus{$3}}, $cds, $1)}
    if ($2 eq "III"){push (@{$III_gmap_cds_locus{$3}}, $cds, $1)}
    if ($2 eq "IV"){push (@{$IV_gmap_cds_locus{$3}}, $cds, $1)}
    if ($2 eq "V"){push (@{$V_gmap_cds_locus{$3}}, $cds, $1)}
    if ($2 eq "X"){push (@{$X_gmap_cds_locus{$3}}, $cds, $1)}
#    push(@{$genetics_mapping{$4}}, $2, $1, $3);  
    push(@{$chrom_pos{$2}}, $3);
  }
}

close FH;

#########################################################
# get each gmap position, its mean of coord on each chrom
# the variables are stored in @chrom_I/II/III/IV/V/X
#########################################################

my (@pos_order, @chrom, @chrom_I, @chrom_II, @chrom_III, @chrom_IV, @chrom_V, @chrom_X, $marker_count);

foreach (sort keys %chrom_pos){
 my $marker_count = 0;
 # print "$_ -> \n";
  if ($_ eq "I"){
    @pos_order = sort {$a <=>$b} @{$chrom_pos{$_}};
    get_gmap_mean_of_each_chrom(\@pos_order, \%I_gmap_cds_locus, \%mean_coord_cds_I, "I", @chrom);
  }

  if ($_ eq "II"){
    @pos_order = sort {$a <=>$b} @{$chrom_pos{$_}};
    get_gmap_mean_of_each_chrom(\@pos_order, \%II_gmap_cds_locus, \%mean_coord_cds_II, "II", @chrom);
  }
#=start
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
    get_gmap_mean_of_each_chrom(\@pos_order, \%X_gmap_cds_locus,  \%mean_coord_cds_X, "X", @chrom);
  }
#=end
#=cut
}

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
  open(OUT, ">>/wormsrv1/chaokung/DOCS/cmp_gmap_with_coord_order") || die $!;

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
  my $inside=0;
  my $outside=0;
  my @num_cds=keys %Mean_coord_cds;
  
  foreach (sort keys %Mean_coord_cds){
    print $_, "--->\n";
    for($i=0; $i < (scalar @chrom)-3; $i=$i+2){
      #print $i," ->", $i+1," -> ",$i+2," -> ",$i+3,"\n";
      $mean_down =  $chrom[$i+1];
      $mean_up = $chrom[$i+3];
      $gmap_up=$chrom[$i+2];
      $gmap_down=$chrom[$i];
      #print $mean_down," ->", $mean_up," -> ",$gmap_down," -> ",$gmap_up,"\n";

      if (($_>$mean_down) && ($_<$mean_up) ){
        if ($gmap_down < 0 && $gmap_up < 0){  # negative gmap position
          $baseline=($mean_up - $mean_down) / (-$gmap_down + $gmap_up);
        }
        if ($gmap_down > 0 && $gmap_up > 0){  # positive gmap position
          $baseline=($mean_up - $mean_down) / ($gmap_up - $gmap_down);
        }
        if ($gmap_down < 0 && $gmap_up > 0){  # positive gmap position
          $baseline=($mean_up - $mean_down) / ($gmap_up - $gmap_down);
        }
        $diff_bp=$_- $mean_down; 
        $position=$diff_bp/$baseline; 
        $position=$gmap_down+$position;
        print $mean_down," ->", $mean_up," -> ",$gmap_down," -> ",$gmap_up,"\n";
        print "\n";
        print $mean_down," ->", $mean_up,">1\n";
        print $gmap_down," -> ",$gmap_up,">2\n";
        print "$_ -> $position#####\n";
        $inside++;
      }
    }
  }
 
  my @total_cds=keys %CDS_mapping;
  print scalar @total_cds, "\n";
  print $inside, "\n";
  print scalar @num_cds, "\n";
  print "Chromosome $chromosom: cds lying outside: ", scalar @num_cds - $inside,"\n";
  @chrom=(); 
  $inside=();
}


__END__
