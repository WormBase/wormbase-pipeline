#!/usr/local/bin/perl
#
# usage: repeatmasker_worm.pl -c [for camace output] OR -s [for stlace output] -f [infile]
#						      -nolow [for NOT masking low complexity regions]
#							  -low   [for JUST masking low complexity regions]		
#
# runs RepeatMasker and creates ace file output
# 11.12.02 Kerstin Jekosch

use strict;
use Ace;
use Getopt::Long;
$|=1;

#############
# variables #
#############

my $campath = '/wormsrv2/camace';
my $stlpath = '/wormsrv2/stlace';
my $outdir  = '/nfs/disk100/wormpub/repeats';
my $tace    = &tace." /wormsrv2/current_DB";
my (%cam,%stl,$infile,$outfile,$opt_c,$opt_s,%acc,$nolow,$low);
&GetOptions( 'c'    => \$opt_c,
		     's'    => \$opt_s,		
             'f=s'  => \$infile,
			 'nolow'=> \$nolow,
			 'low'  => \$low,
			 );

if ($opt_c && $opt_s) {
	print STDERR "Try s OR c, not both at a time!\n";
	exit(0);
}

unless ($infile) {
	print STDERR "Specify infile (-f [file])!";
	exit(0);
}
unless ($low || $nolow)  {
	print STDERR "Specify -low or -nolow\n";
}

open(LOG,">$outdir/log");

#########################
# get accession numbers #
#########################

my $command=<<EOF;
Table-maker -p "/wormsrv2/autoace/wquery/accession2clone.def"
quit
EOF
open (TACE, "echo '$command' | $tace | ");
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
print LOG "Got accession numbers out of autoace\n";

#####################
# get camace clones #
#####################

if ($opt_c) {
	my $db = Ace->connect(-path=>$campath) || die "Cannot open camace",Ace->error;
	print LOG "Checking for clones in camace\n";
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
	print LOG "Checking for clones in stlace\n";
	my @clones = map $_->name, $db->fetch(Genome_sequence => '*') ; 
	foreach my $clone (@clones) {
		$stl{$clone} = 1;
	}
}

########
# MAIN #
########

open(F,"$infile");
print STDERR "Opening infile\n";
while (<F>) {
	my $clone;
	next unless /\>/;

	/\>(\S+)\// and do {
		print LOG "\nDealing with $clone:\n";
		$clone = $1;
		unless ($clone =~ /\S+/) {
			die "Cannot read file $infile\n";
		}

		# create temporary sequence file
		open(TMP,">$outdir/$clone.dna");
		open(FETCH,"pfetch $acc{$clone} |");

		print TMP ">$clone\n";
		while (<FETCH>) {
			next if /\>/;
			print TMP $_;
		}
		close FETCH;
		close TMP;
		
		# check for pfetch success
		open(TEST,"ls -l $outdir | grep $clone.dna |");
		while (<TEST>) {
			if (/\s+0\s+/) {
				print LOG "Pfetch failed for $clone\n";
			}
		} 
		close TEST;
		# run RepeatMasker
		if ($nolow) {
			if (($opt_c && (exists $cam{$clone})) || ($opt_s && (exists $stl{$clone}))) {
				print LOG "Running RepeatMasker -nolow for $clone\n";
				system("RepeatMasker -nolow -lib $outdir/lib/newelegans.lib $outdir/$clone.dna");
				system("mv $outdir/$clone.dna.out $outdir/nolow/$clone.dna.out");
				unlink("$outdir/$clone.dna") if (-e "$outdir/$clone.dna");
				unlink("$outdir/$clone.dna.stderr") if (-e "$outdir/$clone.dna.stderr");
				unlink("$outdir/$clone.dna.cat") if (-e "$outdir/$clone.dna.cat");
				unlink("$outdir/$clone.dna.masked") if (-e "$outdir/$clone.dna.masked");
			}
		}
		elsif ($low) {
			if (($opt_c && (exists $cam{$clone})) || ($opt_s && (exists $stl{$clone}))) {
				print LOG "Running RepeatMasker -int for $clone\n";
				system("RepeatMasker -int $outdir/$clone.dna");
				system("mv $outdir/$clone.dna.out $outdir/low/$clone.dna.out");
				unlink("$outdir/$clone.dna") if (-e "$outdir/$clone.dna");
				unlink("$outdir/$clone.dna.stderr") if (-e "$outdir/$clone.dna.stderr");
				unlink("$outdir/$clone.dna.cat") if (-e "$outdir/$clone.dna.cat");
				unlink("$outdir/$clone.dna.masked") if (-e "$outdir/$clone.dna.masked");
				
			}
			
		}

#		# map to ace file
#		print LOG "Mapping output for $clone\n";
#		map2ace("$outdir/$clone.dna.out","$clone.repeatmasker.ace");
		
	};
}

##########################
# produce final ace file #
##########################

if ($nolow) {
	unlink("$outdir/nolow/all.repeatmasker_nolow.ace") if (-e "$outdir/nolow/all.repeatmasker_nolow.ace");
	system("repeatmasker2ace.pl $outdir/nolow RepeatMasker > $outdir/nolow/all.repeatmasker_nolow.ace");
}
elsif ($low) {
	unlink("$outdir/low/all.repeatmasker_low.ace") if (-e "$outdir/low/all.repeatmasker_low.ace");
	system("repeatmasker2ace.pl $outdir/low low_complexity_repeat > $outdir/low/all.repeatmasker_low.ace");
}

exit(0);

__END__

# replaced internal mapping by upper system call... 

open(FINAL,">$infile.repeatmasker.ace");
print LOG "Producing final outfile\n";
foreach my $file (@files) {
	open(FILE,"$file");
	while (<FILE>) {
		unless (/\S+/) {
			print LOG "$file WAS EMPTY!\n";
		} 
		print FINAL $_;
	}
	unlink("$file");
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
		if ($query =~ /\S+/) {
			print OUT "Sequence \"$query\"\n";
			print OUT "Motif_homol \"$rep\" \"RepeatMasker\" $score $qstart $qend $repstart $repend\n\n";	
		}
	}
}


 
