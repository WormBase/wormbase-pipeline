#!/usr/local/bin/perl5.6.0
#
# evaluates confirmed introns
# 5.10.01 Kerstin Jekosch

use strict;
use Getopt::Std;
use vars qw($opt_c $opt_a $opt_h $opt_v $opt_m);
$|=1;


#######################
# variables and files #
#######################

$opt_a = ""; # produce output for autoace
$opt_c = ""; # produce output for camace
$opt_m = ""; # perform everything for mRNAs (default is EST)
$opt_h = ""; # help
$opt_v = ""; # be verbose
&getopts('achvm');
if ($opt_h) {
	exec('perldoc',$0);
	exit;	
}
my $dir = '/nfs/disk100/wormpub/blat2';
my (%seq,%intron,$temp);

if ($opt_m) {
	if ($opt_a) {
		open(CI,"$dir/autoace.ci.mRNA.ace") or die "Cannot open $dir/autoace.ci.mRNA.ace $!\n";
		open(GOOD,">$dir/autoace.good_introns.mRNA.ace") or die "Cannot open $dir/autoace.good_introns.mRNA.ace $!\n";
		open(BAD,">$dir/autoace.bad_introns.mRNA.ace") or die "Cannot open $dir/autoace.bad_introns.mRNA.ace $!\n";
	}
	elsif ($opt_c) {
		open(CI,"$dir/camace.ci.mRNA.ace") or die "Cannot open $dir/camace.ci.mRNA.ace $!\n";
		open(GOOD,">$dir/camace.good_introns.mRNA.ace") or die "Cannot open $dir/camace.good_introns.mRNA.ace $!\n";
		open(BAD,">$dir/camace.bad_introns.mRNA.ace") or die "Cannot open $dir/camace.bad_introns.mRNA.ace $!\n"; 
	}
	else {
		print STDERR "make a choice: -a for autoace, -c for camace\n";
	}
}
else{
	if ($opt_a) {
		open(CI,"$dir/autoace.ci.EST.ace") or die "Cannot open $dir/autoace.ci.EST.ace $!\n";
		open(GOOD,">$dir/autoace.good_introns.EST.ace") or die "Cannot open $dir/autoace.good_introns.EST.ace $!\n";
		open(BAD,">$dir/autoace.bad_introns.EST.ace") or die "Cannot open $dir/autoace.bad_introns.EST.ace $!\n";
	}
	elsif ($opt_c) {
		open(CI,"$dir/camace.ci.EST.ace") or die "Cannot open $dir/camace.ci.EST.ace $!\n";
		open(GOOD,">$dir/camace.good_introns.EST.ace") or die "Cannot open $dir/camace.good_introns.EST.ace $!\n";
		open(BAD,">$dir/camace.bad_introns.EST.ace") or die "Cannot open $dir/camace.bad_introns.EST.ace $!\n"; 
	}
	else {
		print STDERR "make a choice: -a for autoace, -c for camace\n";
	}
}

#########################################
# loop over file with confirmed introns #
#########################################

my $lala = $/;
$/ = "";
while (<CI>) {
	next unless /^\S/;
	if (/Sequence : \"(\S+)\"/) {
		my $link = $1;
		print "Sequence : $link\n" if ($opt_v);
		my @introns = split /\n/, $_;

		#########################
		# get the link sequence #
		#########################
		
		open(SEQ,"$dir/autoace.fa") or die "Cannot open $dir/autoace.fa $!\n";
		my $switch = 0;
		$/ = $lala;
		my @dna;
		while (<SEQ>) {
			if (/^\>$link$/) {
				$switch = 1;
				print "Started with $link\n" if ($opt_v);
    		}
			elsif (/^(\w+)$/) {
				if ($switch == 1) {
					push @dna, split(//,$1);
				}			
			}
			else { 
				$switch	= 0;		
			}
		}
		
		####################
		# evaluate introns #
		####################
		
		$/ = "";
		foreach my $test (@introns) {
			if ($test =~ /Confirmed_intron/) {
				my @f = split / /, $test;
				
				#######################################
				# get the donor and acceptor sequence #
				#######################################
				
				my ($first,$last,$pastfirst,$prelast);
				if ($f[1] < $f[2]) {
					($first,$last,$pastfirst,$prelast) = ($f[1]-1,$f[2]-1,$f[1],$f[2]-2);
				}
				else {
					($first,$last,$pastfirst,$prelast) = ($f[2]-1,$f[1]-1,$f[2],$f[1]-2);
				}		
				my $start = $dna[$first].$dna[$pastfirst];
				my $end   = $dna[$prelast].$dna[$last];
				print "Coords start $f[1] => $start, end $f[2] => $end\n" if ($opt_v);


				##################
				# map to S_Child #
				##################
				
				my $lastvirt = int((scalar @dna)/100000) + 1;
				my ($startvirt,$endvirt,$virtual);
				if ((int($first/100000) + 1 ) > $lastvirt) {
					$startvirt = $lastvirt;
				}
				else {
					$startvirt = int($first/100000) + 1;
				}
				if ((int($last/100000) + 1 ) > $lastvirt) {
					$endvirt = $lastvirt;
				}
				else {
					$endvirt = int($first/100000) + 1;
				}
				
				if ($startvirt == $endvirt) { 
					$virtual = "Confirmed_intron_EST:".$link."_".$startvirt unless ($opt_m);
					$virtual = "Confirmed_intron_mRNA:".$link."_".$startvirt    if ($opt_m);
				}
				elsif (($startvirt == ($endvirt - 1)) && (($last%100000) <= 50000)) {
					$virtual = "Confirmed_intron_EST:".$link."_".$startvirt unless ($opt_m);
					$virtual = "Confirmed_intron_mRNA:".$link."_".$startvirt    if ($opt_m);
				}


				#################
				# check introns #
				#################
				
				my $one = $f[1]%100000;
				my $two = $f[2]%100000;
				if (((($start eq 'gt') || ($start eq 'gc')) && ($end eq 'ag')) ||
		    		 (($start eq 'ct') && (($end eq 'ac') || ($end eq 'gc')))) {	 
					print GOOD "Feature_data : \"$virtual\"\n";
					print GOOD "Confirmed_intron $one $two EST\n\n";
				}  	
					
				else {
					print BAD "Feature_data : \"$virtual\"\n";
					print BAD "Confirmed_intron $one $two EST\n\n";		
				}
			}
		}
	}
}
