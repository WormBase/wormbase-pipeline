#!/usr/local/bin/perl5.6.1 -w

# cmp_genetics_seq_location.pl

# by Chao-Kung Chen [030113]

# Last updated on: $Date: 2003-02-12 18:37:25 $
# Last updated by: $Author: ck1 $

use strict;                    
use lib "/wormsrv2/scripts/";
use Wormbase;
use Cwd 'chdir';


my $gff_dir = "/wormsrv2/autoace/GFF_SPLITS/";

my @versions=dataset($gff_dir, "folder");

my @order = sort {$a <=> $b} @versions;
#print @order, "\n";

my $gff_location = "/wormsrv2/autoace/GFF_SPLITS/WS"."$order[-1]";
chdir $gff_location;

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
}

my $genetics_mapping=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/genetics_mapping.def" quit
EOF

my $ga_dir="/wormsrv1/geneace";
my %genetics_mapping;
  
open (FH, "echo '$genetics_mapping' | tace $ga_dir | ") || die "Couldn't access Geneace\n";
while (<FH>){
  chomp($_);
  if ($_ =~ /^\"(.+)\"\s+\"(.+)\"\s\"(.+)\"/){
    push(@{$genetics_mapping{$3}}, $2, $1);   # $3=cds, $2=chrom, $1=locus
  }
}
foreach (sort keys %genetics_mapping){
  if (exists ${@{$genetics_mapping{$_}}}[0] && exists ${@{$CDS_mapping{$_}}}[0] ){
    if (${@{$genetics_mapping{$_}}}[0] ne ${@{$CDS_mapping{$_}}}[0]){
      print "$_ is mapped to ${@{$genetics_mapping{$_}}}[0] by genetics\n";
      print "$_ is mapped to ${@{$CDS_mapping{$_}}}[0] by BLAT_BEST.\n\n"
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
