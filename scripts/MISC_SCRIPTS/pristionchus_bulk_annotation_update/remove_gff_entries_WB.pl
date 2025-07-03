#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';
use Getopt::Long;
use Path::Class;

GetOptions(
    "deletions|d=s" => \my $deletions_file,
    "gff|g=s" => \my $gff_file
    );

my %deletions;
my $del_fh = file($deletions_file)->openr;
while (my $line = $del_fh->getline()) {
    chomp $line;
    $deletions{$line}++;
}

my $gff_fh = file($gff_file)->openr;
while (my $line = $gff_fh->getline()) {
    chomp $line;
    next if $line =~ /^#/;
    my @col = split("\t", $line);
    next if $col[1] eq 'history';
    next unless $col[2] eq 'CDS' || $col[2] eq 'gene' || $col[2] eq 'mRNA' || $col[2] eq 'protein_coding_primary_transcript' || $col[2] eq 'pseudogenic_transcript' || $col[2] eq 'five_prime_UTR' || $col[2] eq 'three_prime_UTR' || $col[2] eq 'exon';
    my %attr = split /[;=]/, $col[8];
    my $skip = 0;
    for my $key (qw(ID Parent sequence_name Name)) {
	next unless exists $attr{$key};
	my @value_parts = split(":", $attr{$key});
	my $id;
	if (scalar @value_parts > 1) {
	    $id = $value_parts[1];
	} else {
	    $id = $attr{$key};
	}
	if ($id =~ /^(PPA\d+)/) {
	    $id = $1;
	} elsif ($id =~ /^(WBGene\d+)/) {
	    $id = $1;
	}
	if (exists $deletions{$id}) {
	    $skip = 1;
	    if ($key eq 'Parent' && exists $attr{'ID'}) {
		my @id_parts = split(":", $attr{'ID'});
		if (scalar @id_parts > 1) {
		    $id = $id_parts[1];
		} else {
		    $id = $attr{'ID'};
		}
		if ($id =~ /^(PPA\d+)/) {
		    $id = $1;
		} elsif ($id =~ /^(WBGene\d+)/) {
		    $id = $1;
		}
		$deletions{$id}++;
	    }
	    last;
	}
    }
    next if $skip;
    print $line . "\n";
}
    
