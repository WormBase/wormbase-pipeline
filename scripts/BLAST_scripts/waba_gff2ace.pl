#!/usr/local/bin/perl

# Marc Sohrmann (ms2@sanger.ac.uk)

use strict;
use Getopt::Std;
use vars qw ($opt_e $opt_b $opt_c);

getopts ("e:b:c");

my $usage = "waba_gff2ace.pl\n";
$usage .= "-e [elegans cosmid gff file]\n";
$usage .= "-b [briggsae supercontig fasta file]\n";

unless ($opt_e && $opt_b) {
    die "$usage";
  }

# get the length for all contigs (elegans and briggsae)
my %ce2length;
my %cb2length;

open (CG, "$opt_e") || die "cannot open elegans cosmid gff file $opt_e";
while (<CG>) {
    chomp;
    my @a = split /\t/;
    $a[8] =~ /\"(\S+)\"/;
    my $cos = $1;
    $ce2length{$cos} = $a[4]-$a[3]+1;
  }
close CG;

open (BF, "$opt_b") || die "cannot open briggsae fasta file $opt_b";
my ($id , $seq);
while (<BF>) {
    chomp;
    if (/^>(\S+)/) {
        my $new_id = $1;
	if ($id) {
	    $cb2length{$id} = length($seq);
	}
	$id = $new_id ; $seq = "" ;
    } 
    elsif (eof) {
        if ($id) {
	    $seq .= $_ ;
	    $cb2length{$id} = length($seq);
        }
    }
    else {
        $seq .= $_ ;
    }
}
close BF;

# read the gff file, and print the ace file
# (use s_children to store the waba homols)

my %stateMap = ( "coding" => "coding",
                 "high"   => "strong",
                 "low"    => "weak");

my (%ce_hash, %cb_hash);
while (<>) {
    chomp;
    my @ary = split /\t/;
    # get names, start and end of elegans and briggsae contigs
    my $ce_name = $ary[0];
    my $ce_start = $ary[3];
    my $ce_end = $ary[4];
    my $target = $ary[8];
    my $cb_name;
    my $cb_start;
    my $cb_end;
    my $cigar;
    if ($target =~ /\"\w+:(\S+)\"\s+(\-?\d+)\s+(\-?\d+)\s*(\S*)/) {
        $cb_name = $1;
        $cb_start = $2;
        $cb_end = $3;
	$cigar = $4;
    }
    elsif ($target =~ /\"(\S+)\"\s+(\-?\d+)\s+(\-?\d+)\s*(\S*)/) {
        $cb_name = $1;
        $cb_start = $2;
        $cb_end = $3; 
	$cigar = $4;
    }
    elsif ($target =~ /Sequence\s+(\S+)\s+(\-?\d+)\s+(\-?\d+)\s*(\S*)/) {
        $cb_name = $1;
        $cb_start = $2;
        $cb_end = $3;
	$cigar = $4;
    }
    else {
        die"cannot read the gff file:\n\t$_";
    }
    # check the strand, and modify the coordinates accordingly
    my ($ce_start_strand, $ce_end_strand, $cb_start_strand, $cb_end_strand);
    if ($ary[6] eq "-") {
        my $tmp;
        $tmp = $ce_start;
        $ce_start_strand = $ce_end;
        $ce_end_strand = $tmp;
        $tmp = $cb_start;
        $cb_start_strand = $cb_end;
        $cb_end_strand = $tmp;
    }
    else {
        $ce_start_strand = $ce_start;
        $ce_end_strand = $ce_end;
        $cb_start_strand = $cb_start;
        $cb_end_strand = $cb_end;
    }
    # build the gff lines (as arrays), and store them
    $ary[1] =~ /_(\S+)/;
    my $state = "WABA_".$stateMap{$1};
    unless ($state) {
        die "cannot assign state";
    }
    unless ($opt_c) {
        my @ce_ary = ("DNA_homol", "\"$cb_name\"", $state, $ary[5], $ce_start_strand, $ce_end_strand, $cb_start, $cb_end);
	my @cb_ary = ("DNA_homol", "\"$ce_name\"", $state, $ary[5], $cb_start_strand, $cb_end_strand, $ce_start, $ce_end);
	push (@{$ce_hash{$ce_name}}, \@ce_ary);
	push (@{$cb_hash{$cb_name}}, \@cb_ary);
    }
    else {
        
    }
}

# print elegans ace file
open (CE_ACE, ">ce.ace") || die "cannot create ce.ace";
foreach my $name (sort {$a cmp $b} keys %ce_hash) {
    unless (exists $ce2length{$name}) {
        print STDERR "no length info for $name -> skip it\n";
        next;
    }
    print CE_ACE "\n";
    print CE_ACE "Sequence : \"$name\"\n";
    print CE_ACE "Homol_data \"$name:waba\" 1 $ce2length{$name}\n";
    print CE_ACE "\n";
    print CE_ACE "Homol_data : \"$name:waba\"\n";
    foreach my $aref (@{$ce_hash{$name}}) {
        print CE_ACE join ("\t", @$aref), "\n";
    }
}

# print briggsae ace file
open (CB_ACE, ">cb.ace") || die "cannot create cb.ace";
foreach my $name (sort {$a cmp $b} keys %cb_hash) {
    unless (exists $cb2length{$name}) {
        print STDERR "no length info for $name -> skip it\n";
        next;
    }
    print CB_ACE "\n";
    print CB_ACE "Sequence : \"$name\"\n";
    print CB_ACE "Homol_data \"$name:waba\" 1 $cb2length{$name}\n";
    print CB_ACE "\n";
    print CB_ACE "Homol_data : \"$name:waba\"\n";
    foreach my $aref (sort {$a->[1] cmp $b->[1]} @{$cb_hash{$name}}) {
        print CB_ACE join ("\t", @$aref), "\n";
    }
}

close CE_ACE;
close CB_ACE;
