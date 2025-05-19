#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Path::Class;

GetOptions(
    "in|i=s" => \my $in_file,
    "out|o=s"=> \my $out_file
    );

my $in_fh = file($in_file)->openr;
my %cr_owned_genes;
my @unprocessed_lines;
while (my $line = $in_fh->getline()) {
    push @unprocessed_lines, $line;
    next if $line =~ /^#/;
    chomp $line;
    my @col = split("\t", $line);
    my %attr = split(/[;=]/, $col[8]);
    if ($col[2] eq 'mRNA' || $col[2] eq 'ncRNA') {
	if ((exists $attr{'owner'} && ($attr{'owner'} eq 'christian.roedelsperger@tuebingen.mpg.de' || $attr{'owner'} eq 'christian.roedelsperger@tuebingen.mpg.de,bbraschi@ebi.ac.uk' || $attr{'owner'} eq 'bbraschi@ebi.ac.uk,christian.roedelsperger@tuebingen.mpg.de'))) {
	    $cr_owned_genes{$attr{'Parent'}}++;
	}
    }
}

my $out_fh = file($out_file)->openw;
for my $unprocessed_line (@unprocessed_lines) {
    if ($unprocessed_line =~ /^#/) {
	$out_fh->print($unprocessed_line);
	next;
    }
    chomp $unprocessed_line;
    my @col = split("\t", $unprocessed_line);
    my %attr = split(/[;=]/, $col[8]);
    if ($col[2] eq 'gene' && exists $cr_owned_genes{$attr{'ID'}} && $attr{'owner'} !~ /christian\.roedelsperger/) {
	$attr{'owner'} = 'christian.roedelsperger@tuebingen.mpg.de';
	my $new_line = join("\t", @col[0..7]) . "\t";
	my $first_attribute = 1;
	for my $attribute (keys %attr) {
	    $new_line .= ';' unless $first_attribute;
	    $first_attribute = 0;
	    $new_line .= $attribute . '=' . $attr{$attribute};
	}
	$new_line .= "\n";
	$out_fh->print($new_line);
    } else {
	$out_fh->print($unprocessed_line . "\n");
    }
}
