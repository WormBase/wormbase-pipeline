#!/usr/local/ensembl/bin/perl

# Marc Sohrmann (ms2@sanger.ac.uk)

# map gff coordinates of features in file 1 to gff coordinates in file 2
# (both gff files must have the same parent sequences in the 1st field)

use strict;
use Getopt::Long;
use vars qw($file1 $file2 $condition1 $condition2 $sort);

my $usage = "map_gff.pl\n";
$usage .= "-f1        [first gff file]\n";
$usage .= "-f2        [second gff file]\n";
$usage .= "-c1        [regular expression matching 2nd field in the first gff file]\n";
$usage .= "-c2        [regular expression matching 2nd field in the second gff file]\n";
$usage .= "-sort      sort gff files by start and end coordinates\n";
$usage .= "\nIMPORTANT: both gff files must be sorted\n";
$usage .= "           by start and end coordinates,\n";
$usage .= "           or sorted using the -sort option\n\n";

die "$usage" unless (GetOptions("f1=s"    => \$file1,
                                "f2=s"    => \$file2,
                                "c1=s"    => \$condition1,
                                "c2=s"    => \$condition2,
			        "sort"    => \$sort));

unless ($file1 && $file2) {
    die "$usage";
}

open (F1, "$file1") || die "cannot open $file1";
open (F2, "$file2") || die "cannot open $file2";

#######################################
# read the first gff file
my @one_tmp = ();
if ($sort) {
    my (@a, @s, @e);
    while (<F1>) {
        chomp;
        next if /^#/;
        my @f = split /\t/;
        push @a, $_;
        push @s, $f[3];
        push @e, $f[4];
    }
    foreach my $i (sort {$s[$a] <=> $s[$b] or $e[$a] <=> $e[$b]} 0..$#a) {
        push (@one_tmp, $a[$i]);
    }
}
else {
    while (<F1>) {
        chomp;
        next if /^#/;
        push (@one_tmp, $_);
    }
}
close F1;

my %one = ();
my %one_index = ();
foreach my $line (@one_tmp) {
    my @ary = split (/\t/, $line);
    if ($condition1) {
        next unless $ary[1] =~ /$condition1/;
	    }
    $one_index{$ary[0]} = 0 unless defined $one_index{$ary[0]};
    $one{$ary[0]}->[$one_index{$ary[0]}]->{start} = $ary[3];
    $one{$ary[0]}->[$one_index{$ary[0]}]->{end} = $ary[4];
    $one{$ary[0]}->[$one_index{$ary[0]}]->{gff} = $line;
    $one_index{$ary[0]}++;
}

@one_tmp = ();

#######################################
# read and process the second gff file
# (keeping the array index of the left match with the highest end coordinate)
my @two_tmp = ();
if ($sort) {
    my (@a, @s, @e);
    while (<F2>) {
        chomp;
        next if /^#/;
        my @f = split /\t/;
        push @a, $_;
        push @s, $f[3];
        push @e, $f[4];
    }
    foreach my $i (sort {$s[$a] <=> $s[$b] or $e[$a] <=> $e[$b]} 0..$#a) {
        push (@two_tmp, $a[$i]);
    }
}
else {
    while (<F2>) {
        chomp;
        next if /^#/;
        push (@two_tmp, $_);
    }
}
close F2;

my %two = ();
my %two_index = ();
my %big_end = ();
my %big_index = ();

foreach my $line (@two_tmp) {
    my @ary = split (/\t/, $line);
    if ($condition2) {
        next unless $ary[1] =~ /$condition2/;
    }
    $two_index{$ary[0]} = 0 unless defined $two_index{$ary[0]};
    $big_index{$ary[0]} = 0 unless defined $big_index{$ary[0]};
    $two{$ary[0]}->[$two_index{$ary[0]}]->{start} = $ary[3];
    $two{$ary[0]}->[$two_index{$ary[0]}]->{end} = $ary[4];
    $two{$ary[0]}->[$two_index{$ary[0]}]->{gff} = $line;
    $two{$ary[0]}->[$two_index{$ary[0]}]->{right} = $big_index{$ary[0]};
    if ($ary[4] > $big_end{$ary[0]}) {
        $big_end{$ary[0]} = $ary[4];
        $big_index{$ary[0]} = $two_index{$ary[0]};
    }
    $two_index{$ary[0]}++;
}

@two_tmp = ();

#######################################
# do the mapping

open (FAILED,">failed_mappings") or die "cant open failure file $!\n";
my @parents_list = ();
push (@parents_list, keys %one); push (@parents_list, keys %two);
my %seen = ();
my @parents = grep { ! $seen{$_}++ } @parents_list;
foreach my $parent (sort {$a cmp $b} @parents) {
    print STDERR "processing $parent\n";
    my $starting_index = 0;
    # loop over all entries from the first gff file
    MATCH:foreach my $match (@{$one{$parent}}) {
        my $map = "";
        my $length = 1000000000000000;
        # calculate the new starting index for the search
        while ($starting_index > 0) {
            $starting_index = $two{$parent}->[$starting_index]->{right};
            last if $two{$parent}->[$starting_index]->{end} < $match->{start};
	}
        # loop over all possible entries from the second gff file
        LOOP:for (my $i = $starting_index ; $i <= $#{$two{$parent}} ; $i++) {
            # we have gone too far...
            if ($two{$parent}->[$i]->{start} > $match->{end}) {
                $starting_index = $i unless $map;
                last LOOP;
	    }
            # 1 within 2
            if (($match->{start} >= $two{$parent}->[$i]->{start} && $match->{end} <= $two{$parent}->[$i]->{end}) && (($two{$parent}->[$i]->{end})-($two{$parent}->[$i]->{start}) < $length)) {
                $starting_index = $i unless $map;
                $length = ($two{$parent}->[$i]->{end})-($two{$parent}->[$i]->{start});
                my @one = split (/\t/, $match->{gff});
                my @two = split (/\t/, $two{$parent}->[$i]->{gff});
                my $new_start = ($match->{start})-($two{$parent}->[$i]->{start})+1;
                my $new_end = ($match->{end})-($two{$parent}->[$i]->{start})+1;
		my $new_parent;
		if ($two[8] =~ /\"(\S+)\"/) {
		    $new_parent = $1;
		}
		else {
		    $new_parent = $two[8];
		}
                my $strand = "";
                if ($two[6] eq "-") {
                    $strand = ($one[6] eq "+") ? "-" : "+";
		}
                else {
                    $strand = $one[6];
		}
                $map = "$new_parent\t$one[1]\t$one[2]\t$new_start\t$new_end\t$one[5]\t$strand\t$one[7]\t$one[8]";
	    }
	}
        # sucessfull mapping
        if ($map) {
            print "$map\n";
	}
        else {
            print FAILED $match->{gff},"\n";
	}
    }
}
close FAILED;
