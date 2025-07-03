#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Path::Class;

GetOptions(
    "new|n=s" => \my $new_gff_file,
    "old|o=s" => \my $old_gff_file,
    "out_prefix|p=s" => \my $out_prefix
    );

my ($cr_coords, $cr_ids) = process_new_gff($new_gff_file);

my %overlapping_genes;
my $fh = file($old_gff_file)->openr;
while (my $line = $fh->getline()) {
    next if $line =~ /^#/;
    chomp $line;
    my @col = split("\t", $line);
    next unless $col[2] eq 'gene';
    next unless exists $cr_coords->{$col[0]};
    next unless exists $cr_coords->{$col[0]}{$col[6]};
    for my $cr_start (sort {$a<=>$b} keys %{$cr_coords->{$col[0]}{$col[6]}}) {
	next if $cr_start > $col[4]; # CR gene starts after gene ends 
	for my $cr_end (sort {$a<=>$b} keys %{$cr_coords->{$col[0]}{$col[6]}{$cr_start}}) {
	    next if $cr_end < $col[3]; # CR gene ends before gene starts
	    my %attr = split /[;=]/, $col[8];
	    if (exists $attr{'sequence_name'}) {
		push @{$overlapping_genes{$cr_coords->{$col[0]}{$col[6]}{$cr_start}{$cr_end}}},$attr{'sequence_name'};
	    } elsif (exists $attr{'Name'}) {
		push @{$overlapping_genes{$cr_coords->{$col[0]}{$col[6]}{$cr_start}{$cr_end}}},$attr{'Name'};
	    }
	}
    }
}

my $summary_fh = file($out_prefix . 'summary.csv')->openw;
$summary_fh->print("ID,Position,Old models,Five prime\n");
for my $cr_id (keys %{$cr_ids}) {
    $summary_fh->print($cr_id . ',' . $cr_ids->{$cr_id} . ',');
    if (exists $overlapping_genes{$cr_id} && scalar @{$overlapping_genes{$cr_id}} > 0) {
	$summary_fh->print(scalar @{$overlapping_genes{$cr_id}} . ',' . join("|", @{$overlapping_genes{$cr_id}}));
	my ($strand) = $cr_ids->{$cr_id} =~ /\(([\+\-])\)$/;
	my $five_prime_overlapping_gene;
	if ($strand eq '+') {
	    $five_prime_overlapping_gene = ${$overlapping_genes{$cr_id}}[0];
	} else {   
	    $five_prime_overlapping_gene = ${$overlapping_genes{$cr_id}}[-1];
	}
	$summary_fh->print(',' . $five_prime_overlapping_gene . "\n");
    } else {
	$summary_fh->print("0\n");
    }
}

sub process_new_gff {
    my $gff_file = shift;
    my %coords;
    my %ids;
    
    my $fh = file($gff_file)->openr;
    while (my $line = $fh->getline()) {
	next if $line =~ /^#/;
	chomp $line;
	my @col = split("\t", $line);
	next unless $col[2] eq 'gene';
	my %attr = split /[;=]/, $col[8];
	my $id;
	if ($col[1] eq 'WormBase') {
	    $id = $attr{'sequence_name'};
	    $id = $attr{'Name'} if !defined $id;
	} elsif ($col[1] eq 'iris') {
	    $id = $attr{'Name'};
	} else {
	    $id = $attr{'ID'};
	}
	if (!defined $id) {
	    print $line. "\n";
	    next;
	}
	$coords{$col[0]}{$col[6]}{$col[3]}{$col[4]} = $id;
	$ids{$id} = $col[0] . ':' . $col[3] . '-' . $col[4] . ' (' . $col[6] . ')';
    }
    return (\%coords, \%ids);
}
