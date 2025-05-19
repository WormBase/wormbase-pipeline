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
    if ($line =~ /^#/) {
	print $line . "\n";
    }
    my @col = split("\t", $line);
    my %attr = split /[;=]/, $col[8];
    my $skip = 0;
    if (exists $attr{'ID'}) {
	my ($id) = $attr{'ID'} =~/^(PACO\d+)\.?/;
	next if exists $deletions{$id};
    }
    if (exists $attr{'Parent'}) {
	my ($id) = $attr{'Parent'} =~/^(PACO\d+)\.?/;
	next if exists $deletions{$id};
    }
    print $line . "\n";
}
    
