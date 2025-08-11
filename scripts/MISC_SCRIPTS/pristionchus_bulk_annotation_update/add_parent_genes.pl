#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Path::Class;

GetOptions(
    "in|i=s" => \my $in_file,
    "out|o=s" => \my $out_file
    );

my $in_fh = file($in_file)->openr;
my $out_fh = file($out_file)->openw;
while (my $line = $in_fh->getline()) {
    chomp $line;
    my @col = split("\t", $line);
    my %attr = split(/[;=]/, $col[8]);
    if ($col[2] eq 'mRNA' || $col[2] eq 'ncRNA') {
	if (!exists $attr{'Parent'}) {
	    my ($parent) = $attr{'ID'} =~ /^(PACO\d+)\./;
	    $out_fh->print(join("\t", $col[0], $col[1], 'gene', $col[3], $col[4], $col[5], $col[6], $col[7], 'ID=' . $parent) . "\n");
	    $out_fh->print($line . ';Parent=' . $parent . "\n");
	} else {
	    $out_fh->print($line . "\n");
	}
    } else {
	$out_fh->print($line ."\n");
    }
}
