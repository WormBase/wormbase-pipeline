#!/usr/local/bin/perl
#
# usage: repeatmasker_worm.pl -c [for camace output and/or] -s [for stlace output] 
#							  -m [just map existing RepeatMasker file to ace file]	
#
# runs RepeatMasker and creates ace file output
# 11.12.02 Kerstin Jekosch

use strict;
use Ace;
use Getopt::Long;

#############
# variables #
#############

my $campath = '/wormsrv2/camace';
my $stlpath = '/wormsrv2/stlace';
my $outdir  = '/nfs/disk100/wormpub/repeats';
my (%cam,%stl,$infile,$outfile,$opt_c,$opt_s,$opt_m,%acc);
&GetOptions( 'c'    => \$opt_c,
		     's'    => \$opt_s,		
             'm'    => \$opt_m,
             'f=s' => \$infile);

if ($opt_c && $opt_s) {
	print STDERR "Try s OR c, not both at a time!\n";
	exit(0);
}

unless ($infile) {
	print STDERR "Specify infile (-f [file])!";
	exit(0);
}

#########################
# get accession numbers #
#########################

my $autotace  = "/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/tace /wormsrv2/current_DB";
my $command1=<<EOF;
Table-maker -p "/wormsrv2/autoace/wquery/accession2clone.def"
quit
EOF
open (TACE, "echo '$command1' | $autotace | ");
while (<TACE>) {
	chomp;
	next if ($_ eq "");
    next if (/\/\//);
    s/acedb\>\s//g;
    s/\"//g;
	s/EMBL://g;
	my ($name,$acc) = ($_ =~ /^(\S+)\s(\S+)/);
	$acc{$name} = $acc;
}
close TACE;




#####################
# get camace clones #
#####################

if ($opt_c) {
	my $db = Ace->connect(-path=>$campath) || die "Cannot open camace",Ace->error;
	print STDERR "Checking for clones in camace\n";
	my @clones = map $_->name, $db->fetch(Genome_sequence => '*') ; 
	foreach my $clone (@clones) {
		$cam{$clone} = 1;
	}
}

#####################
# get stlace clones #
#####################

if ($opt_s) {
	my $db = Ace->connect(-path=>$stlpath) || die "Cannot open stlace",Ace->error;
	print STDERR "Checking for clones in stlace\n";
	my @clones = map $_->name, $db->fetch(Genome_sequence => '*') ; 
	foreach my $clone (@clones) {
		$stl{$clone} = 1;
	}
}

######################
# start RepeatMasker #
######################

open(F,"$infile");
print STDERR "Opening infile\n";
while (<F>) {
	my ($clone,$seq);
	next unless /\>/;
	/\>(\S+)\// and do {
		$clone = $1;
		unless ($clone =~ /\S+/) {
			die "Cannot read file $infile\n";
		}

		open(TMP,">$outdir/$clone.dna");
		open(FETCH,"pfetch $acc{$clone} |");

		print TMP ">$clone\n";
		while (<FETCH>) {
			next if /\>/;
			print TMP $_;
		}
		close FETCH;
	};
	close TMP;
	
	# run RepeatMasker
	if (($opt_c && (exists $cam{$clone})) || ($opt_s && (exists $stl{$clone}))) {
		print "Running RepeatMasker for $clone\n";
		system("RepeatMasker -nolow -lib $outdir/newelegans.lib $outdir/$clone.dna");
	}
		
	# map to ace file
	print STDERR "Mapping output for $clone\n";
	map2ace("$outdir/$clone.dna.out","$clone.repeatmasker.ace");
} 




###############
# subroutines #
###############

sub map2ace {

	my $file     = $_[0];
	my $outfile  = $_[1];
	open(IN,"$file");
	open(OUT,">$outfile");
	my ($score,$query,$qstart,$qend,$rep,$repfam,$try1,$try2,$try3);
	my ($repstart,$repend);

	while (<IN>) {
		/(\d+)\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+\(\d+\)\s+\S\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/ and do {
			($score,$query,$qstart,$qend,$rep,$repfam,$try1,$try2,$try3) = ($1,$2,$3,$4,$5,$6,$7,$8,$9);
			if ($try1 =~ /\(/) {
				$repstart = $try2;
				$repend   = $try3; 
			}
			else {
				$repstart = $try1;
				$repend   = $try2; 
			}
		};
		print OUT "Sequence \"$query\"\n";
		print OUT "Motif_homol \"$rep\" \"RepeatMasker\" $score $qstart $qend $repstart $repend\n\n";	
	}
}


 
