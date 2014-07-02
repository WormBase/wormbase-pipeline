#!/usr/bin/env perl
# while convert a GFF3 file into BED
# usage:
# 	perl blat2bed.pl INFILE > OUTFILE
#
# note:
# 	it will convert all lines keyed by ID or Parent, so make sure 
# 	the input file contains only the lines you want converted.
#
#	does not work (yet) with unnamed features a.e. repeats, 
#	as it will merge them into one BED line

my %hits;

while(<>){
	chomp;
	if (/(ID|Parent)=([^,;=\t]+)/){ # left before right in regex
		$hits{$2}||=[];
		push (@{$hits{$2}},$_)
	}else{
	   next
	}
}

while (my($k,$v)=each %hits){
	my ($chrom,$start,$stop,$strand);
	my @blockStarts;
	my @blockSizes;
	my $blockCount = scalar @$v;
        my @blocks;
	my $name = $k;
        $name = $1 if $$v[0]=~/Target=([^,;=\s]+)/; # for blat
 
	map {my @t = split "\t";push @blocks, \@t} @$v;
	foreach my $block( sort {$$a[3] <=> $$b[3]} @blocks){
	  $chrom = $$block[0];
	  $strand = $$block[6];
	  $start ||= $$block[3];
	  $start = $$block[3] if $$block[3] < $start;
	  $stop  ||= $$block[4];
	  $stop = $$block[4] if $$block[4] > $stop;
          push @blockSizes ,($$block[4] - $$block[3] +1);
	}
       
	foreach my $block( sort {$$a[3] <=> $$b[3]} @blocks){
	   push @blockStarts, ($$block[3]-$start);
	}
        $start--;
	print "$chrom\t$start\t$stop\t$name\t1000\t$strand\t$start\t$stop\t0,0,255\t$blockCount\t".join(',',@blockSizes)."\t".join(',',@blockStarts)."\n";
}
