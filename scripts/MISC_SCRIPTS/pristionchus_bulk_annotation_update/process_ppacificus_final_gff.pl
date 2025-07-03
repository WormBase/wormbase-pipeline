#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Path::Class;

GetOptions(
    "new|n=s" => \my $new_gff_file,
    "wb|w=s" => \my $wb_gff_file,
    "id|i=s" => \my $id_map_file,
    "out|o=s" => \my $out_file
    );

my ($wb_transcript_coords, $wb_max_transcript_nr) = get_existing_transcript_details($wb_gff_file);
my ($new_transcript_coords) = get_new_transcript_details($new_gff_file);
my $id_map = parse_id_map_file($id_map_file);

my $in_fh = file($new_gff_file)->openr;
my $out_fh = file($out_file)->openw;
my (%wb_gene_id_name_map, %used_gene_ids, %po_id_map);
my $tmp_id_nr = 100000;
while (my $line = $in_fh->getline()) {
    chomp $line;
    
    my @col = split("\t", $line);
    my %attr = split(/[;=]/, $col[8]);

    if ($col[1] eq 'WormBase') {
	if ($col[2] eq 'gene') {
	    my $new_attr = 'ID=' . $attr{'sequence_name'};
	    $wb_gene_id_name_map{$attr{'ID'}} = $attr{'sequence_name'};
	    $out_fh->print(join("\t", $col[0], 'WormBase/pristionchus.org', $col[2], $col[3], $col[4], $col[5], $col[6], $col[7], $new_attr) . "\n");
	    $used_gene_ids{$attr{'sequence_name'}}++;
	} elsif ($col[2] eq 'mRNA') {
	    my $new_attr = 'ID=' . $attr{'Name'} . ';Parent=' . $wb_gene_id_name_map{$attr{'Parent'}};
	    $out_fh->print(join("\t", $col[0], 'WormBase/pristionchus.org', $col[2], $col[3], $col[4], $col[5], $col[6], $col[7], $new_attr) . "\n");
	} else {
	    my $parent = $attr{'Parent'};
	    $parent =~ s/^Transcript://;
	    $parent =~ s/^Pseudogene://;
	    my $new_attr = 'Parent=' . $parent;
	    $out_fh->print(join("\t", $col[0], 'WormBase/pristionchus.org', $col[2], $col[3], $col[4], $col[5], $col[6], $col[7], $new_attr) . "\n");
	}
    } elsif ($col[1] eq 'pristionchus.org') {
	if ($col[2] eq 'gene') {
	    my $ppa_gene_id;
	    if (exists $id_map->{$attr{'ID'}} && $id_map->{$attr{'ID'}} =~ /^PPA/ && !exists $used_gene_ids{$id_map->{$attr{'ID'}}}) {
		$ppa_gene_id = $id_map->{$attr{'ID'}};
	    } else {
		$ppa_gene_id = 'tmpPPA' . $tmp_id_nr;
		$tmp_id_nr++;
	    }
	    my $new_attr = 'ID=' . $ppa_gene_id;
	    $po_id_map{$attr{'ID'}} = $ppa_gene_id;
	    $used_gene_ids{$ppa_gene_id}++;
	    $out_fh->print(join("\t", $col[0], 'WormBase/pristionchus.org', $col[2], $col[3], $col[4], $col[5], $col[6], $col[7], $new_attr) . "\n");
	} elsif ($col[2] eq 'mRNA') {
	    my $gene_id = $attr{'Parent'};
	    my $ppa_gene_id = $po_id_map{$gene_id};
	    my $new_coord_string = $new_transcript_coords->{$attr{'Parent'}}{$attr{'ID'}};
	    my $ppa_transcript_id;
	    if (exists $wb_transcript_coords->{$ppa_gene_id}) {
		for my $transcript_id (sort keys %{$wb_transcript_coords->{$ppa_gene_id}}) {
		    if ($new_coord_string eq $wb_transcript_coords->{$ppa_gene_id}{$transcript_id}) {
			$ppa_transcript_id = $transcript_id;
			last;
		    }
		}
	    }
	    if (!defined $ppa_transcript_id) {
		$wb_max_transcript_nr->{$ppa_gene_id}++;
		$ppa_transcript_id = $ppa_gene_id . '.' . $wb_max_transcript_nr->{$ppa_gene_id};
	    }
	    my $new_attr = 'ID=' . $ppa_transcript_id . ';Parent=' . $ppa_gene_id;
	    $po_id_map{$attr{'ID'}} = $ppa_transcript_id;
	    $out_fh->print(join("\t", $col[0], 'WormBase/pristionchus.org', $col[2], $col[3], $col[4], $col[5], $col[6], $col[7], $new_attr) . "\n");
	} else {
	    my $new_attr = 'Parent=' . $po_id_map{$attr{'Parent'}};
	    $out_fh->print(join("\t", $col[0], 'WormBase/pristionchus.org', $col[2], $col[3], $col[4], $col[5], $col[6], $col[7], $new_attr) . "\n");
	}
    } else {
	if ($col[2] eq 'gene') {
	    my $transcript = $attr{'Name'};
	    my $gene;
	    if ($transcript =~ /^(PPA\d+)/) {
		$gene = $1;
	    } elsif ($transcript =~ /^(PACO\d+)/) {
		$gene = $1;
	    } else {
		$gene = $transcript;
	    }

	    my $ppa_gene_id;
	    if (exists $id_map->{$gene} && !exists $used_gene_ids{$gene}) {
		$ppa_gene_id = $id_map->{$gene};
	    } else {
		$ppa_gene_id = 'tmpPPA' . $tmp_id_nr;
		$tmp_id_nr++;
	    }
	    my $new_attr = 'ID=' . $ppa_gene_id;
	    $po_id_map{$attr{'ID'}} = $ppa_gene_id;
	    $used_gene_ids{$ppa_gene_id}++;
	    $out_fh->print(join("\t", $col[0], 'WormBase/pristionchus.org', $col[2], $col[3], $col[4], $col[5], $col[6], $col[7], $new_attr) . "\n");
	} elsif ($col[2] eq 'mRNA' || $col[2] eq 'ncRNA') {
	    if ($attr{'Name'} eq '24845929-3c28-4813-a138-cbf0900a3d01') {
		$attr{'Name'} = 'PPA01653.1-000001';
	    } elsif ($attr{'Name'} eq 'new_gene_122344345') {
		$attr{'Name'} = 'PPA32540.1-000001';
	    }
	    elsif ($attr{'Name'} eq 'PPA36517.1B') {
		$attr{'Name'} = 'STRG.6469a-000001';
	    }
	    elsif ($attr{'Name'} eq 'PPA_ChrV:17948024..17957786   ') {
		$attr{'Name'} = 'PPA_ChrV-000001';
	    }
	    my ($transcript) = $attr{'Name'} =~ /^(.+)\-\d+$/;
	    print "TRANSCRIPT: " . $attr{'Name'} . "\n" unless defined $transcript;
	    my $gene;
	    if ($transcript =~ /^(PPA\d+)/) {
		$gene = $1;
	    } elsif ($transcript =~ /^(PACO\d+)/) {
		$gene = $1;
	    } else {
		$gene = $transcript;
	    }
	    my $ppa_gene_id = $po_id_map{$attr{'Parent'}};
	    my $new_coord_string = $new_transcript_coords->{$gene}{$transcript};
	    my $ppa_transcript_id;
	    if (defined $new_coord_string && exists $wb_transcript_coords->{$ppa_gene_id}) {
		for my $transcript_id (sort keys %{$wb_transcript_coords->{$ppa_gene_id}}) {
		    if ($new_coord_string eq $wb_transcript_coords->{$ppa_gene_id}{$transcript_id}) {
			$ppa_transcript_id = $transcript_id;
			last;
		    }
		}
	    }
	    if (!defined $ppa_transcript_id) {
		$wb_max_transcript_nr->{$ppa_gene_id}++;
		$ppa_transcript_id = $ppa_gene_id . '.' . $wb_max_transcript_nr->{$ppa_gene_id};
	    }
	    my $new_attr = 'ID=' . $ppa_transcript_id . ';Parent=' . $ppa_gene_id;
	    $po_id_map{$attr{'ID'}} = $ppa_transcript_id;
	    $out_fh->print(join("\t", $col[0], 'WormBase/pristionchus.org', $col[2], $col[3], $col[4], $col[5], $col[6], $col[7], $new_attr) . "\n");
	} else {
	    my $new_attr = 'Parent=' . $po_id_map{$attr{'Parent'}};
	    $out_fh->print(join("\t", $col[0], 'WormBase/pristionchus.org', $col[2], $col[3], $col[4], $col[5], $col[6], $col[7], $new_attr) . "\n");
	}
    }
}
    
sub parse_id_map_file {
    my $file = shift;

    my %map;
    my $fh = file($file)->openr;
    while (my $line = $fh->getline()) {
	chomp $line;
	my ($new_id, $ppa_id) = split(",", $line);
	$map{$new_id} = $ppa_id;
    }

    return \%map;
}

sub get_existing_transcript_details {
    my $file = shift;

    my (%coords, %max_nr);
    my $fh = file($file)->openr;
    while (my $line = $fh->getline()) {
	chomp $line;
	next if $line =~ /^#/;
	my @col = split("\t", $line);
	my %attr = split(/[;=]/, $col[8]);
	
	if ($col[2] eq 'pseudogenic_transcript') {
	    my ($gene) = $attr{'ID'} =~ /^Pseudogene:(.+)$/;
	    $coords{$gene}{$gene} = $col[3] . '-' . $col[4];
	}
	next unless $col[2] eq 'CDS';
	my ($gene) = $attr{'ID'} =~ /^CDS:(.+)$/;
	my @transcripts = split(',', $attr{'Parent'});
	for my $transcript (@transcripts) {
	    $transcript =~ s/Transcript://g;
	    my ($transcript_nr) = $transcript =~ /\.(\d+)\D?$/;
	    if (exists $coords{$gene} && exists $coords{$gene}{$transcript}) {
		$coords{$gene}{$transcript} .= '|' . $col[3] . '-' . $col[4];
	    } else {
		$coords{$gene}{$transcript} = $col[3] . '-' . $col[4];
	    }
	    if (!exists $max_nr{$gene} || $max_nr{$gene} < $transcript_nr) {
		$max_nr{$gene} = $transcript_nr;
	    }
	}
    }

    return (\%coords, \%max_nr);
}

sub get_new_transcript_details {
    my $file = shift;

    my (%coords, %transcript_gene_map, %iris_map);
    my $fh = file($file)->openr;
    while (my $line = $fh->getline()) {
	chomp $line;
	my @col = split("\t", $line);
	my %attr = split(/[;=]/, $col[8]);

	if ($col[2] eq 'ncRNA') {
	    my ($gene) = $attr{'Name'} =~ /^(.+)\-\d+$/;
	    $coords{$gene}{$gene} = $col[3] . '-' . $col[4];
	    next;
	}
	if ($col[1] eq 'pristionchus.org') {
	    if ($col[2] eq 'mRNA') {
		$transcript_gene_map{$attr{'ID'}} = $attr{'Parent'};
	    }
	    if ($col[2] eq 'CDS') {
		my $transcript = $attr{'Parent'};
		my $gene = $transcript_gene_map{$transcript};
		if (exists $coords{$gene} && exists $coords{$gene}{$transcript}) {
		    $coords{$gene}{$transcript} .= '|' . $col[3] . '-' . $col[4];
		} else {
		    $coords{$gene}{$transcript} = $col[3] . '-' . $col[4];
		}
	    }
	} elsif ($col[1] eq 'iris') {
	    if ($col[2] eq 'mRNA') {
		if ($attr{'Name'} eq '24845929-3c28-4813-a138-cbf0900a3d01') {
		    $attr{'Name'} = 'PPA01653.1-000001';
		}
		elsif ($attr{'Name'} eq 'new_gene_122344345') {
		    $attr{'Name'} = 'PPA32540.1-000001';
		}
		elsif ($attr{'Name'} eq 'PPA36517.1B') {
		    $attr{'Name'} = 'STRG.6469a-000001';
		}
		elsif ($attr{'Name'} eq 'PPA_ChrV:17948024..17957786   ') {
		    $attr{'Name'} = 'PPA_ChrV-000001';
		}
		my ($transcript) = $attr{'Name'} =~ /^(.+)\-\d+$/;
		print "TRANSCRIPT: " . $attr{'Name'} . "\n" unless defined $transcript;
		my $gene;
		
		if ($transcript =~ /^(PPA\d+)/) {
		    $gene = $1;
		} elsif ($transcript =~ /^(PACO\d+)/) {
		    $gene = $1;
		} else {
		    $gene = $transcript;
		}
		my $transcript_id = $attr{'ID'};
		$iris_map{$transcript_id}{'transcript'} = $transcript;
		$iris_map{$transcript_id}{'gene'} = $gene;
	    } elsif ($col[2] eq 'CDS') {
		my $gene = $iris_map{$attr{'Parent'}}{'gene'};
		print "GENE: " . $attr{'Parent'} ."\n"unless $gene;
		my $transcript = $iris_map{$attr{'Parent'}}{'transcript'};
		if (exists $coords{$gene} && $coords{$gene}{$transcript}) {
		    $coords{$gene}{$transcript} .= '|' . $col[3] . '-' . $col[4];
		} else {
		    $coords{$gene}{$transcript} = $col[3] . '-' . $col[4];
		}
	    }
	}
    }
	
    return (\%coords);
}
