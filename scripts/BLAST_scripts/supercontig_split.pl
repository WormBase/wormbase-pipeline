#!/usr/local/bin/perl

use strict;

my $split_length = 100000;

open (GFF, ">split.gff") || die "cannot create gff file";
open (FASTA, ">split.fasta") || die "cannot create fasta file";

my ($id , $seq);
while (<>) {
    chomp;
    if (/^>(\S+)/) {
	my $new_id = $1;
	if ($id) {
	    my $it = int((length($seq))/$split_length);
	    for (my $i = 0; $i <= $it; $i++) {
		my $split_seq = substr ($seq, $i*$split_length, $split_length);
		print FASTA ">$id.split_$i";
		my $tmp_split_seq = reverse $split_seq;
		my $count = 0;
		while (my $base = chop $tmp_split_seq) {
		    if ($count++ % 50 == 0) {
			print FASTA "\n";
		    }
		    print FASTA $base;
		}
		print FASTA "\n";
		my $length = length($split_seq);
		print GFF "$id\tsplit\tsequence\t".(($i*$split_length)+1)."\t".(($i*$split_length)+$length)."\t.\t+\t.\t\"$id.split_$i\"\n";
	    }
	}
	$id = $new_id ; $seq = "" ;
    } 
    elsif (eof) {
	if ($id) {
	    $seq .= $_ ;
	    my $it = int((length($seq))/$split_length);
	    for (my $i = 0; $i <= $it; $i++) {
		my $split_seq = substr ($seq, $it*$split_length, $split_length);
		my $count = 0;
		print FASTA ">$id.split_$i";
		my $tmp_split_seq = reverse $split_seq;
		while (my $base = chop $tmp_split_seq) {
		    if ($count++ % 50 == 0) {
			print FASTA "\n";
		    }
		    print FASTA $base;
		}
		print FASTA "\n";
		my $length = length($split_seq);
		print GFF "$id\tsplit\tsequence\t".(($i*$split_length)+1)."\t".(($i*$split_length)+$length)."\t.\t+\t.\t\"$id.split_$i\"\n";
	    }
	}
    }
    else {
	$seq .= $_ ;
    }
}
close GFF;
close FASTA;
