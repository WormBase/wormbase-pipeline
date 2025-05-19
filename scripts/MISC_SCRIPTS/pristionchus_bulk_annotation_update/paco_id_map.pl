#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Path::Class;

GetOptions (
    "cr|c=s" => \my $cr_id_map_file,
    "gff|g=s" => \my $gff_file,
    "wbgff|w=s" => \my $wb_gff_file,
    "map|m=s" => \my $out_file,
    "unmapped|u=s" => \my $unmapped_file,
    "overlap|o=s" => \my $overlap_file
    );

my $paco_id_to_wb_name_map = parse_cr_id_map_file($cr_id_map_file);
my $wb_id_to_name_map = parse_wb_gff($wb_gff_file);
my $overlap_map = parse_overlap_file($overlap_file);

my $out_fh = file($out_file)->openw;
my @unmapped;
my $gff_fh = file($gff_file)->openr;
my %used_ids;
while (my $line = $gff_fh->getline()) {
    chomp $line;
    my @col = split("\t", $line);
    next unless ($col[2] eq 'mRNA' && $col[1] eq 'pristionchus.org') || ($col[2] eq 'gene' && $col[1] ne 'pristionchus.org');
    my %attr = split(/[;=]/, $col[8]);
    if ($col[1] eq 'WormBase') {
	if (exists $attr{'sequence_name'} && !exists $used_ids{$attr{'sequence_name'}}) {
	    $out_fh->print($attr{'sequence_name'} . ',' . $attr{'sequence_name'} . "\n");
	    $used_ids{$attr{'sequence_name'}}++;
	} else {
	    push @unmapped, $line;
	}
    } elsif ($col[1] eq 'iris') {
	my $wb_name;
	my $incoming_id;
	if (exists $attr{'Name'}) {
	    if ($attr{'Name'} =~ /^(WBGene\d+)/) {
		$incoming_id = $1;
		$wb_name = $wb_id_to_name_map->{$incoming_id};
	    } elsif ($attr{'Name'} =~ /^(PACO\d+)/) {
		$incoming_id = $1;
		if (exists $paco_id_to_wb_name_map->{$incoming_id}) {
		    $wb_name = $paco_id_to_wb_name_map->{$incoming_id};
		}
	    } elsif ($attr{'Name'} =~ /^(PPA\d+)/) {
		$wb_name = $1;
		$incoming_id = $1;
	    } else {
		$incoming_id = $attr{'Name'};
	    }
	}
	if (!defined $wb_name && exists $overlap_map->{$attr{'Name'}}) {
	    $wb_name = $overlap_map->{$attr{'Name'}};
	}
	if (defined $wb_name and !exists $used_ids{$wb_name}) {
	    $out_fh->print($incoming_id . ',' . $wb_name . "\n");
	    $used_ids{$wb_name}++;
	} else {
	    if (defined $incoming_id && !exists $used_ids{$incoming_id}) {
		$out_fh->print($incoming_id . ',' . $incoming_id . "\n");
		$used_ids{$incoming_id}++;
	    } else {
		$out_fh->print($attr{'Name'} . ',' . $attr{'Name'} . "\n");
		$used_ids{$attr{'Name'}}++;
	    }
	    push @unmapped, $line;
	}
    } else {
	my $wb_name;
	if (exists $attr{'ID'}) {
	    my ($gene_id) = $attr{'ID'} =~ /^(PACO\d+)/; 
	    if (exists $paco_id_to_wb_name_map->{$gene_id}) {
		$wb_name = $paco_id_to_wb_name_map->{$gene_id};
	    }
	    if (!defined $wb_name && exists $overlap_map->{$gene_id}) {
		$wb_name = $overlap_map->{$gene_id};
	    }
	    
	    if (defined $wb_name and !exists $used_ids{$wb_name}) {
		$out_fh->print($gene_id . ',' . $wb_name . "\n");
		$used_ids{$wb_name}++;
	    } else {
		$out_fh->print($gene_id . ',' . $gene_id . "\n");
		push @unmapped, $line;
	    }
	} else {
	    push @unmapped, $line;
	}
    }
}

my $unmapped_fh = file($unmapped_file)->openw;
for my $unmapped_line (@unmapped) {
    $unmapped_fh->print($unmapped_line . "\n");
}

sub parse_cr_id_map_file {
    my $file = shift;
    my $in_fh = file($file)->openr;
    my %map;
    while (my $line = $in_fh->getline()) {
	chomp $line;
	my ($paco_id, $wb_name) = split(" ", $line);
	$map{$paco_id} = $wb_name; 
    }

    return \%map;
}

sub parse_wb_gff {
    my $file = shift;

    my %map;
    my $in_fh = file($file)->openr;
    while (my $line = $in_fh->getline()) {
	chomp $line;
	next if $line =~ /^#/;
	my @col = split("\t", $line);
	next unless $col[2] eq 'gene';
	my %attr = split(/[;=]/, $col[8]);
	next unless exists $attr{'Name'} && exists $attr{'sequence_name'};
	$map{$attr{'Name'}} = $attr{'sequence_name'};
    }

    return \%map;
}

sub parse_overlap_file {
    my $file = shift;

    my %map;
    my $fh = file($file)->openr;
    while (my $line = $fh->getline()) {
	chomp $line;
	my @col = split(",", $line);
	next unless scalar @col == 5;
	$map{$col[0]} = $col[4];
    }

    return \%map;
}
		 
