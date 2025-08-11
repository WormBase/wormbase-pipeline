#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Path::Class;

GetOptions(
    "gff|g=s" => \my $gff_file,
    "keep|k=s" => \my $keep_file,
    "out|o=s" => \my $out_file
    );

my $keep_fh = file($keep_file)->openr;
my (%genes_to_keep, %transcripts_to_keep);
while (my $line = $keep_fh->getline()) {
    chomp $line;
    next unless $line =~ /\-\d+$/;
    $transcripts_to_keep{$line}++;
    my ($gene_id) = $line =~ /^(.+)\-\d+$/;
    $genes_to_keep{$1}++;
}

my %transcript_ids_to_keep;
my %transcripts_found;
my $gff_fh = file($gff_file)->openr;
my $out_fh = file($out_file)->openw;
while (my $line = $gff_fh->getline()) {
    next if $line =~ /^#/;
    chomp $line;
    my @col = split("\t", $line);
    my %attr = split(/[;=]/, $col[8]);
    if ($col[2] eq 'gene') {
	if ((exists $attr{'owner'} && ($attr{'owner'} eq 'christian.roedelsperger@tuebingen.mpg.de' || $attr{'owner'} eq 'christian.roedelsperger@tuebingen.mpg.de,bbraschi@ebi.ac.uk' || $attr{'owner'} eq 'bbraschi@ebi.ac.uk,christian.roedelsperger@tuebingen.mpg.de')) || (exists $attr{'Name'} && exists $genes_to_keep{$attr{'Name'}})) {
	    $out_fh->print($line . "\n");
	}
    } elsif ($col[2] eq 'mRNA' || $col[2] eq 'ncRNA') {
	if ((exists $attr{'owner'} && ($attr{'owner'} eq 'christian.roedelsperger@tuebingen.mpg.de' || $attr{'owner'} eq 'christian.roedelsperger@tuebingen.mpg.de,bbraschi@ebi.ac.uk' || $attr{'owner'} eq 'bbraschi@ebi.ac.uk,christian.roedelsperger@tuebingen.mpg.de')) || (exists $attr{'Name'} && exists $transcripts_to_keep{$attr{'Name'}})) {
	    $out_fh->print($line . "\n");
	    $transcripts_found{$attr{'Name'}} = 1 if exists $attr{'Name'};
	    $transcript_ids_to_keep{$attr{'ID'}} = 1 if exists $attr{'ID'};
	}
    } else {
	if (exists $attr{'Parent'} && exists $transcript_ids_to_keep{$attr{'Parent'}}) {
	    $out_fh->print($line . "\n");
	}
    }
}

for my $transcript (keys %transcripts_to_keep) {
    print "Could not find $transcript\n" unless exists $transcripts_found{$transcript};
}
