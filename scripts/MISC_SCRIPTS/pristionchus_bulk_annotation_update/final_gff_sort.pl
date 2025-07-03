#!/usr/bin/perl
use strict;
use warnings;

use Path::Class;
use Getopt::Long;

GetOptions(
    "in|i=s" => \my $in_file,
    "out|o=s" => \my $out_file
    );

my (%genes, %transcripts, %other, %gene_ids, %transcript_ids);
my $in_fh = file($in_file)->openr;
while (my $line = $in_fh->getline()) {
    chomp $line;
    my @col = split("\t", $line);
    my %attr = split(/[;=]/, $col[8]);

    if ($col[2] eq 'gene') {
	if (exists $genes{$col[0]}{$col[3]}{$col[4]}{$col[6]}) {
	    print "DUPLICATE POSITION FOR " . $attr{'ID'} . ' and ' . $genes{$col[0]}{$col[3]}{$col[4]}{$col[6]}{'ID'} . "\n";
	}
	$genes{$col[0]}{$col[3]}{$col[4]}{$col[6]}{'ID'} = $attr{'ID'};
	$gene_ids{$attr{'ID'}}++;
	$genes{$col[0]}{$col[3]}{$col[4]}{$col[6]}{'line'} = $line;
    } elsif ($col[2] eq 'mRNA' || $col[2] eq 'ncRNA') {
	if (exists $transcripts{$attr{'Parent'}}{$col[3]}{$col[4]}{$attr{'ID'}}) {
	    print "DUPLICATE TRANSCRIPT FOR " . $attr{'ID'} . "\n";
	}
	$transcripts{$attr{'Parent'}}{$col[3]}{$col[4]}{$attr{'ID'}} = $line;
	$transcript_ids{$attr{'ID'}}++;
    } else {
	$other{$attr{'Parent'}}{$col[3]}{$col[4]}{$col[6]}{$col[2]} = $line;
    }
}

for my $gene_id(keys %transcripts) {
    if (!exists $gene_ids{$gene_id}) {
	print "TRANSCRIPT WITHOUT PARENT " . $gene_id . "\n";
    }
}

for my $transcript_id(keys %other) {
    if (!exists $transcript_ids{$transcript_id}) {
	print "CDS/EXON WITHOUT PARENT " . $transcript_id . "\n";
    }
}
my $out_fh = file($out_file)->openw;
for my $chr (sort keys %genes) {
    for my $gene_start (sort {$a<=>$b} keys %{$genes{$chr}}) {
	for my $gene_end (sort {$a<=>$b} keys %{$genes{$chr}{$gene_start}}) {
	    for my $gene_strand (sort keys %{$genes{$chr}{$gene_start}{$gene_end}}) {
		my $gene_id = $genes{$chr}{$gene_start}{$gene_end}{$gene_strand}{'ID'};
		$out_fh->print($genes{$chr}{$gene_start}{$gene_end}{$gene_strand}{'line'} . "\n");
		for my $transcript_start (sort {$a<=>$b} keys %{$transcripts{$gene_id}}) {
		    for my $transcript_end (sort {$b<=>$a} keys %{$transcripts{$gene_id}{$transcript_start}}) {
			for my $transcript_id (sort keys %{$transcripts{$gene_id}{$transcript_start}{$transcript_end}}) {
			    $out_fh->print($transcripts{$gene_id}{$transcript_start}{$transcript_end}{$transcript_id} . "\n");
			    for my $other_start (sort {$a<=>$b} keys %{$other{$transcript_id}}) {
				for my $other_end (sort {$b<=>$a} keys %{$other{$transcript_id}{$other_start}}) {
				    for my $other_strand (sort keys %{$other{$transcript_id}{$other_start}{$other_end}}) {
					for my $other_type (sort keys %{$other{$transcript_id}{$other_start}{$other_end}{$other_strand}}) {
					    $out_fh->print($other{$transcript_id}{$other_start}{$other_end}{$other_strand}{$other_type} . "\n");
					}
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
}
