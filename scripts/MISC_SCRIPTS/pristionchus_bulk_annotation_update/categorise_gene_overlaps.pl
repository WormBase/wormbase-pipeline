#!/usr/bin/env perl
use strict;
use warnings;
use List::Util qw(min max);

my ($new_file, $old_file) = @ARGV;

my $old_genes = parse_gff($old_file);
my $new_genes = parse_gff($new_file);

# Create hash for lookup
my (%old_by_id, %new_by_id);
$old_by_id{ $_->{id} } = $_ for @$old_genes;
$new_by_id{ $_->{id} } = $_ for @$new_genes;

# Build overlaps
my %overlap_matrix;
foreach my $old (@$old_genes) {
    foreach my $new (@$new_genes) {
        next unless $old->{chr} eq $new->{chr} && $old->{strand} eq $new->{strand};
        my $olap = overlap_length($old, $new);
        next unless $olap > 0;
        $overlap_matrix{$old->{id}}{$new->{id}} = $olap;
    }
}

# Track usage
my (%used_old, %used_new);
my (@merges, @splits, @updates, @creations, @deletions);

# Step 1: Merges (many old -> one new)
foreach my $new_id (sort keys %new_by_id) {
    next if $used_new{$new_id};
    my @old_matches = grep { !$used_old{$_} && defined $overlap_matrix{$_}{$new_id} } keys %old_by_id;
    next unless @old_matches >= 2;
    my @sorted = sort { $old_by_id{$a}{start} <=> $old_by_id{$b}{start} } @old_matches;
    next unless are_sequential(\@sorted, \%old_by_id);

    push @merges, { new => $new_id, old => [@sorted] };
    $used_new{$new_id}++;
    $used_old{$_}++ for @sorted;
}

# Step 2: Splits (one old -> many new)
foreach my $old_id (sort keys %old_by_id) {
    next if $used_old{$old_id};
    my @new_matches = grep { !$used_new{$_} && defined $overlap_matrix{$old_id}{$_} } keys %new_by_id;
    next unless @new_matches >= 2;
    my @sorted = sort { $new_by_id{$a}{start} <=> $new_by_id{$b}{start} } @new_matches;
    next unless are_sequential(\@sorted, \%new_by_id);

    push @splits, { old => $old_id, new => [@sorted] };
    $used_old{$old_id}++;
    $used_new{$_}++ for @sorted;
}

# Step 3: Updates (1-to-1 only)
foreach my $old_id (keys %old_by_id) {
    next if $used_old{$old_id};
    my @new_matches = grep { !$used_new{$_} && defined $overlap_matrix{$old_id}{$_} } keys %new_by_id;
    next unless @new_matches == 1;
    my $new_id = $new_matches[0];
    my @rev = grep { !$used_old{$_} } keys %old_by_id;
    my @rev_matches = grep { defined $overlap_matrix{$_}{$new_id} } @rev;
    next unless @rev_matches == 1;

    push @updates, { old => $old_id, new => $new_id };
    $used_old{$old_id}++;
    $used_new{$new_id}++;
}

# Step 4: Creations
foreach my $new_id (keys %new_by_id) {
    next if $used_new{$new_id};
    push @creations, $new_id;
}

# Step 5: Deletions
foreach my $old_id (keys %old_by_id) {
    next if $used_old{$old_id};
    push @deletions, $old_id;
}

# Output
write_table("calculated_merges.txt", \@merges, ["new", "old"]);
write_table("calculated_splits.txt", \@splits, ["old", "new"]);
write_table("calculated_updates.txt", \@updates, ["old", "new"]);
write_table("calculated_creations.txt", \@creations, ["new"]);
write_table("calculated_deletions.txt", \@deletions, ["old"]);

print "Classification complete.\n";

### ------------------------- Subroutines -------------------------

sub parse_gff {
    my ($file) = @_;
    open my $fh, '<', $file or die "Cannot open $file: $!";
    my (%genes, %mRNAs, %coords, %strand, %chr);
    while (<$fh>) {
        chomp;
        next if /^#/;
        my @F = split("\t");
        my ($seqid, $type, $start, $end, $strand_val, $attrs) = @F[0,2,3,4,6,8];
        my %attr = map { split /=/, $_, 2 } split /;/, $attrs;

        if ($type eq 'gene') {
            my $id = ($attr{ID} =~ s/^Gene://r);
	    $genes{$id} = 1;
            $chr{$id} = $seqid;
            $strand{$id} = $strand_val;
        } elsif ($type eq 'mRNA') {
            my $mid = ($attr{ID} =~ s/^Transcript://r);
            my $pid = ($attr{Parent} =~ s/^Gene://r);
            $mRNAs{$mid} = $pid;
	} elsif ($type eq 'CDS') {
            my $mid = $attr{Parent};
	    $mid =~ s/Transcript://g;
	    my @mids = split(/,/, $mid);
	    for my $m(@mids) {
		my $pid = $mRNAs{$m} // next;
		push @{ $coords{$pid} }, [$start, $end];
	    }
	}
    }
    close $fh;

    my @result;
    foreach my $id (keys %genes) {
        my $c = $coords{$id} || [];
        my ($min, $max) = (1e9, 0);
        for (@$c) {
            $min = $_->[0] if $_->[0] < $min;
            $max = $_->[1] if $_->[1] > $max;
        }
        next if $min > $max;
        push @result, { id => $id, chr => $chr{$id}, strand => $strand{$id}, start => $min, end => $max };
    }
    return \@result;
}

sub overlap_length {
    my ($a, $b) = @_;
    return 0 unless $a->{chr} eq $b->{chr} && $a->{strand} eq $b->{strand};
    my $s = max($a->{start}, $b->{start});
    my $e = min($a->{end}, $b->{end});
    return $e >= $s ? $e - $s + 1 : 0;
}

sub are_sequential {
    my ($ids, $lookup) = @_;
    my @starts = map { $lookup->{$_}{start} } @$ids;
    for my $i (1 .. $#starts) {
        return 0 if $starts[$i] < $starts[$i-1];
    }
    return 1;
}

sub write_table {
    my ($file, $data, $cols) = @_;
    open my $fh, '>', $file or die "Cannot write to $file: $!";
    print $fh '#' . join("\t", @$cols), "\n";
    for my $row (@$data) {
        if (ref $row eq 'HASH') {
            my @vals = map {
                ref $row->{$_} eq 'ARRAY' ? join(",", @{ $row->{$_} }) : $row->{$_}
            } @$cols;
            print $fh join("\t", @vals), "\n";
        } else {
            print $fh "$row\n";
        }
    }
    close $fh;
}
