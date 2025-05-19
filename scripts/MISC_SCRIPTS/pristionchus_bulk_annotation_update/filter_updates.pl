#!/usr/bin/perl
use strict;
use warnings;
use Path::Class;
use Const::Fast;
use Getopt::Long;

GetOptions(
    "wormbase|w=s" => \my $wb_gff_file,
    "new|n=s" => \my $new_gff_file,
    "updates|u=s" => \my $updates_file,
    "out|o=s"=> \my $out_file
    );

my ($new_genes, $new_transcripts) = parse_gff($new_gff_file);
my ($wb_genes, $wb_transcripts) = parse_gff($wb_gff_file);

my $updates_fh = file($updates_file)->openr;
my $out_fh = file($out_file)->openw;
while (my $line = $updates_fh->getline()) {
    chomp $line;
    if ($line =~ /^#/) {
	$out_fh->print($line . "\n");
	next;
    }

    my ($wb_gene, $new_gene) = split("\t", $line);
    my $updated = 0;
    if ($wb_genes->{$wb_gene} ne $new_genes->{$new_gene}) {
	$updated = 1;
    }
    if (scalar keys %{$wb_transcripts->{$wb_gene}} != scalar keys %{$new_transcripts->{$new_gene}}) {
	$updated = 1;
    }
    for my $wb_transcript (keys %{$wb_transcripts->{$wb_gene}}) {
	if (!exists $new_transcripts->{$new_gene}{$wb_transcript}) {
	    $updated = 1;
	}
    }

    if ($updated == 1) {
	$out_fh->print($line ."\n");
    }
}

sub parse_gff {
    my $file = shift;

    my (%genes, %transcripts);
    
    my $fh = file($file)->openr;
    while (my $line = $fh->getline()) {
	next if $line =~ /^#/;
	chomp $line;
	my @col = split("\t", $line);
	my %attr = split(/[;=]/, $col[8]);

	my $chr = $col[0];
	my $feature = $col[2];
	my $start = $col[3];
	my $end = $col[4];
	my $strand = $col[6];

	my $parent = exists $attr{'Parent'} ? $attr{'Parent'} : '';
	my $id = exists $attr{'ID'} ? $attr{'ID'} : '';

	$id =~ s/^Gene://;
	$parent =~ s/^Gene://;
	$id =~ s/^Transcript://;
	$parent =~ s/^Transcript://;
	$id =~ s/^Pseudogene://;
	$parent =~ s/^Pseudogene://;

	if ($feature eq 'gene') {
	    die "Expecting ID for line:\n$line\n" if $id eq '';
	    $genes{$id} = $chr . ':' . $start . '-' . $end . '(' . $strand. ')';
	} elsif ($feature eq 'mRNA' || $feature eq 'ncRNA' || $feature eq 'pseudogenic_transcript') {
	    die "Expecting ID and parent for line:\n$line\n" if $id eq '' || $parent eq '';
	    push @{$transcripts{$parent}{$id}}, $start . '|' . $end;
	} 
    }

    return (\%genes, \%transcripts);
}

