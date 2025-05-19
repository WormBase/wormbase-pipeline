#!/usr/bin/perl
use strict;
use warnings;
use Path::Class;
use Getopt::Long;
GetOptions(
    "in|i=s" => \my $in_file,
    "out|o=s" => \my $out_prefix
    );

my $in_fh = file($in_file)->openr;
my $keep_fh = file($out_prefix .'to_keep.txt')->openw;
my $delete_fh = file($out_prefix . 'to_delete.txt')->openw;

while (my $line = $in_fh->getline()) {
    next if $line =~ /^#/;
    chomp $line;
    my ($delete_str, $keep_str) = split(",", $line);
    my @delete_list = split(/\|/, $delete_str);
    my @keep_list = split(/\|/, $keep_str);
    my %keep_map = map {$_ => 1} @keep_list;
    for my $delete (@delete_list) {
	next if exists $keep_map{$delete};
	$delete_fh->print($delete . "\n");
    }
    for my $keep (@keep_list) {
	$keep_fh->print($keep . "\n");
    }
}
