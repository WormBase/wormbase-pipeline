#!/usr/local/bin/perl


=head1 NAME 

halfwise.pl 

=head1 SYNOPSIS

perl halfwise.pl < current.versions 


=head1 DESCRIPTION

submits halfwise jobs for worm analysis. builds stuff
~wormpub/analysis/cosmids directories: firstly by
making a directory cosmid/data/wise if not there and
then by submitting a job to the nightq. It does not
submit if wise/halfwise.out is there.

=cut

use strict;

my $root = "/nfs/disk100/wormpub/analysis/cosmids/";
my $clone;
my $clonename;
my $bsubopts = "-q nightq";

while( <> ) {
    ($clone) = split;
    $clone =~ m#(\w+)/\d+# || do { print STDERR "$clone doesn't look good!\n"; next; };
    $clonename = $1;
#    system("rm -rf $root/$clone/wise");

    if( ! -d "$root$clone/wise" ) {
	print STDERR "Making wise dir! $root/$clone/wise \n";
	if( system("mkdir $root$clone/wise") != 0 ) {
	    print STDERR "Could not make dir for $clone\n";
	    next;
	}
    } 

    if( ! -e "$root/$clone/wise/seq.fasta" ) {
	if( system("ln -s $root/$clone/$clonename.fasta $root/$clone/wise/seq.fasta") != 0 ) {
	    print STDERR "Could not make link for $clone\n";
	    next;
	}
    }
    if( ! -e "$root/$clone/wise/seq.fasta.masked" ) {
	chdir("$root/$clone/wise/");
	print "No $root/$clone/wise/seq.fasta.masked\n";
	if( system("/nfs/disk100/humpub/scripts/RepeatMasker $root/$clone/wise/seq.fasta") != 0 ) {
	    print STDERR "Could not make repeat masker for $clone\n";
	    next;
	}
    }


    system("bsub $bsubopts -o $root/$clone/wise/bsub.output 'halfwise $root/$clone/wise/seq.fasta.masked -gene worm.gf -init wing -pretty -caceh -pseudo -cut 25 -aln 200 -quiet > $root/$clone/wise/halfwise.out'");
}








