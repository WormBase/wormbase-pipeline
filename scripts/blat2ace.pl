#!/usr/local/bin/perl5.6.0 
#
# Exporter to map blat data to acedb and to find the best match for each EST
# 010905 by Kerstin Jekosch

use strict;
use Ace;
use Getopt::Std;
use vars qw($opt_i $opt_h $opt_a $opt_c $opt_m $opt_x);

$| = 1;

#############################
# variables and directories #
#############################
 
my $dir	      = "/wormsrv2/autoace/BLAT";
my $dbdir     = "/wormsrv2/autoace";
my $tace      = "/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/tace /wormsrv1/camace";
my $autotace  = "/nfs/disk100/acedb/RELEASE.DEVELOPMENT/bin.ALPHA_4/tace /wormsrv2/current_DB";
my (%ESTs,%camace,%hash,%best,%other,%bestclone,%match,%ci,%dir);
open (LOG, ">$dir/blat2ace.log");

$opt_i=""; # get confirmed introns
$opt_a=""; # get best matches for autoace
$opt_c=""; # get best matches for camace
$opt_m=""; # do it for mRNAs (default is ESTs)  
$opt_x=""; # do it for parasitic nematode ESTs (blatx)
$opt_h=""; # print help
getopts ('ihacmx');
&printhelp if ($opt_h);

#################
# get EST names #
#################

unless ($opt_m || $opt_x) {
	print LOG "Making table query for EST names\n\n"; 
	my $command1=<<EOF;
	Table-maker -p "$dbdir/wquery/ESTacc2names.def"
	quit
EOF
	open (TACE, "echo '$command1' | $tace | ");
	while (<TACE>) {
		chomp;
		next if ($_ eq "");
    	next if (/\/\//);
    	s/acedb\>\s//g;
    	s/\"//g;
		s/EMBL://g;
		my ($acc,$name) = ($_ =~ /^(\S+)\s(\S+)/);
		unless ($name) {
			$name = $acc;
		} 
		$ESTs{$acc} = $name;
	}
	close TACE;
}

######################################
# get orientation of ESTs (.3 or .5) #
######################################

if ($opt_i && !$opt_m && !$opt_x) {
	print LOG "Making table query for EST names\n\n"; 
	my $command2=<<EOF;
	Table-maker -p "$dbdir/wquery/ESTorient.def"
	quit
EOF
	open (TACE, "echo '$command2' | $tace | ");
	while (<TACE>) {
		chomp;
		next if ($_ eq "");
    	next if (/\/\//);
    	s/acedb\>\s//g;
		s/\"//g;
		my ($name,$orient) = ($_ =~ /^(\S+)\s+EST_(\d)/);
		if ($orient) {
			$dir{$name} = $orient;
		}
	}
	close TACE;
}

###########################
# get all links in camace #
###########################

my @camclones = qw(cTel3X cTel4X CTEL7X cTel33B CTEL54X LINK_6R55 LINK_cTel52S SUPERLINK_CB_I SUPERLINK_CB_II SUPERLINK_CB_IIIL SUPERLINK_CB_IIIR SUPERLINK_CB_IR SUPERLINK_CB_IV SUPERLINK_CB_V SUPERLINK_CB_X); 
foreach my $camclone (@camclones) {
	$camace{$camclone} = 1;
}

############################
# map the blat hits to ace #
############################

print LOG "Start mapping\n\n";
if ($opt_m) {
	open(BLAT,"$dir/mRNA_out.psl");
	open(ACE,">$dir/autoace.mRNA.ace");
	open(CCE,">$dir/camace.mRNA.ace") if ($opt_c);
}
elsif ($opt_x) {
	open(BLAT,"$dir/nematode_out.psl");
	open(ACE,">$dir/autoace.blat.nematode.ace");
	open(CCE,">$dir/camace.blat.nematode.ace") if ($opt_c);
}
else {
	open(BLAT,"$dir/est_out.psl");
	open(ACE,">$dir/autoace.EST.ace");
	open(CCE,">$dir/camace.EST.ace") if ($opt_c);
}
while (<BLAT>) {
    next unless (/^\d/);

    my @f         = split "\t";
    my $superlink = $f[13];
	my $slsize    = $f[14];
	my $lastvirt  = int($slsize/100000) + 1; 
	
    
    #############################################################
    # replace EST name (usually accession number) by yk... name #
    #############################################################
	
    my $est = $f[9];
    if ((!$opt_m) && (!$opt_x) && (exists $ESTs{$est})) {
		my $estname  = $ESTs{$est};
		if ($est ne $estname) {
	    	$est = $estname;
	    	print LOG "EST name $est was replaced by $estname\n\n";
		}
    }
    my @lengths     = split (/,/, $f[18]);
    my @eststarts   = split (/,/, $f[19]);
    my @slinkstarts = split (/,/, $f[20]);

	my $matchstart  = $f[15];
	my $matchend    = $f[16];

	###############################
	# find virtual superlink part #
	###############################
	
	my ($virtual,$startvirtual,$endvirtual);
	if ((int($matchstart/100000) +1) > $lastvirt) { $startvirtual = $lastvirt;}
	else {$startvirtual = int($matchstart/100000) +1;}  
	
	if ((int($matchend/100000) +1) > $lastvirt) { $endvirtual = $lastvirt;}
	else {$endvirtual = int($matchend/100000) +1;}  
	
	if ($startvirtual == $endvirtual) {
		$virtual = "BLAT_EST:".$superlink."_".$startvirtual unless ($opt_m || $opt_x);
		$virtual = "BLAT_mRNA:".$superlink."_".$startvirtual    if ($opt_m);
		$virtual = "BLATX_NEMATODE:".$superlink."_".$startvirtual   if ($opt_x);
	}	
	elsif (($startvirtual == ($endvirtual - 1)) && (($matchend%100000) <= 50000)) {
		$virtual = "BLAT_EST:".$superlink."_".$startvirtual unless ($opt_m || $opt_x);
		$virtual = "BLAT_mRNA:".$superlink."_".$startvirtual    if ($opt_m);
		$virtual = "BLATX_NEMATODE:".$superlink."_".$startvirtual   if ($opt_x);
	}
	else {
		print LOG "$est wasn't assigned to a virtual object as match size was too big\n";
		print LOG "Start is $matchstart, end is $matchend on $superlink\n\n";
		next;
	}
    
    ###################
    # calculate score #
    ###################
	
	my $sum 		= 0;
    foreach my $length (@lengths) {
		$sum = $sum + $length;
    }	
    my $match = $f[0];
    my $score = $match/$sum*100;
    
    my @exons = ();

    #########################
    # calculate coordinates #
    #########################
	
	# need to allow for est exons in the next virtual object, otherwise they get remapped to the start 
	# of the virtual by performing %100000
	my $calc = int(($slinkstarts[0]+1)/100000);
	
	for (my $x = 0;$x < $f[17]; $x++) {
		my $newcalc      = int(($slinkstarts[$x]+1)/100000);
		my $virtualstart;
		if ($calc == $newcalc) {	
			$virtualstart = ($slinkstarts[$x] +1)%100000;
		}
		elsif ($calc == ($newcalc-1)) {
		        $virtualstart = (($slinkstarts[$x] +1)%100000) + 100000;
		}
		my $virtualend   = $virtualstart + $lengths[$x] -1;
		my ($eststart,$estend);
		if ($opt_x) {
			my $temp;
			if (($f[8] eq '++') || ($f[8] eq '-+')) {
	    		$eststart   = $eststarts[$x] +1;
	    		$estend     = $eststart + $lengths[$x] -1;
				if ($f[8] eq '-+') {
					$temp     = $estend;
					$estend   = $eststart;
					$eststart = $temp; 
				}
			}
			elsif (($f[8] eq '--') || ($f[8] eq '+-')) {
				$temp         = $virtualstart;
				$virtualstart = $virtualend;
				$virtualend   = $temp;
	    		$eststart     = $f[10] - $eststarts[$x];
	    		$estend       = $eststart - $lengths[$x] +1;
				if ($f[8] eq '--') {
					$temp     = $estend;
					$estend   = $eststart;
					$eststart = $temp; 
				}
			}			
		}
		else {
			if ($f[8] eq '+'){
	    		$eststart   = $eststarts[$x] +1;
	    		$estend     = $eststart + $lengths[$x] -1;
			}
			elsif ($f[8] eq '-') {
	    		$eststart   = $f[10] - $eststarts[$x];
	    		$estend     = $eststart - $lengths[$x] +1;
			}		
		}		
		print LOG "$est was mapped to $virtual\n\n";
		print ACE "Homol_data : \"$virtual\"\n";
		printf ACE "DNA_homol\t\"%s\"\t\"BLAT_EST_OTHER\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$est,$score,$virtualstart,$virtualend,$eststart,$estend unless ($opt_m || $opt_x);
		printf ACE "DNA_homol\t\"%s\"\t\"BLAT_mRNA_OTHER\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$est,$score,$virtualstart,$virtualend,$eststart,$estend     if ($opt_m);
		printf ACE "DNA_homol\t\"%s\"\t\"BLATX_NEMATODE\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$est,$score,$virtualstart,$virtualend,$eststart,$estend    if ($opt_x);
		if (($opt_c) && (exists $camace{$superlink})) {
	    	print CCE "Homol_data : \"$virtual\"\n";
	    	printf CCE "DNA_homol\t\"%s\"\t\"BLAT_EST_OTHER\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$est,$score,$virtualstart,$virtualend,$eststart,$estend unless ($opt_m|| $opt_x);
	    	printf CCE "DNA_homol\t\"%s\"\t\"BLAT_mRNA_OTHER\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$est,$score,$virtualstart,$virtualend,$eststart,$estend    if ($opt_m);
	    	printf CCE "DNA_homol\t\"%s\"\t\"BLATX_NEMATODE\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$est,$score,$virtualstart,$virtualend,$eststart,$estend   if ($opt_x);
		}

		push @exons, [$virtualstart,$virtualend,$eststart,$estend];				
    }
	
   	########################
	# collect best matches #
	########################
	
    if (exists $best{$est}) {
		if ($score >= $best{$est}->{'score'}) {
	    	if ( ($score > $best{$est}->{'score'}) || ($match > $best{$est}->{'match'})) { 
	    		$best{$est}->{'score'} = $score;
				$best{$est}->{'match'} = $match;
				@{$best{$est}->{'entry'}} = ({'clone' => $virtual,'link' => $superlink,'exons' => \@exons});
	    	}
	    	elsif ($match == $best{$est}->{'match'}) {
	    		$best{$est}->{'score'} = $score;
				push @{$best{$est}->{'entry'}}, {'clone' => $virtual,'link' => $superlink,'exons' => \@exons};
	    	}
		}
    }
    else {
		$best{$est}->{'match'} = $match;
		$best{$est}->{'score'} = $score;
		@{$best{$est}->{'entry'}} = ({'clone' => $virtual,'link' => $superlink,'exons' => \@exons});
    }
}
close BLAT;
close ACE;
close CCE if ($opt_c);

    
####################################
# produce outfile for best matches #
####################################

if ($opt_m) {
	open(AUTBEST,">$dir/autoace.best.mRNA.ace");
	open(CAMBEST,">$dir/camace.best.mRNA.ace") if ($opt_c);
}
else {
	open(AUTBEST,">$dir/autoace.best.EST.ace");
	open(CAMBEST,">$dir/camace.best.EST.ace") if ($opt_c);
}
if (!$opt_x) {
	foreach my $found (sort keys %best) {
    	if (exists $best{$found}) {
			foreach my $entry (@{$best{$found}->{'entry'}}) {
				if (@{$best{$found}->{'entry'}} < 2) {
					my $virtual   = $entry->{'clone'};
					my $superlink = $entry->{'link'};
					foreach my $ex (@{$entry->{'exons'}}) {
						my $score        = $best{$found}->{'score'};
						my $virtualstart = $ex->[0];
						my $virtualend   = $ex->[1];
						my $eststart     = $ex->[2];
						my $estend       = $ex->[3];
						print  AUTBEST "Homol_data : \"$virtual\"\n";
						printf AUTBEST "DNA_homol\t\"%s\"\t\"BLAT_EST_BEST\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$found,$score,$virtualstart,$virtualend,$eststart,$estend unless ($opt_m || $opt_x);
						printf AUTBEST "DNA_homol\t\"%s\"\t\"BLAT_mRNA_BEST\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$found,$score,$virtualstart,$virtualend,$eststart,$estend    if ($opt_m);
						if (($opt_c) && (exists $camace{$superlink})) {
		    				print  CAMBEST "Homol_data : \"$virtual\"\n";
		    				printf CAMBEST "DNA_homol\t\"%s\"\t\"BLAT_EST_BEST\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$found,$score,$virtualstart,$virtualend,$eststart,$estend unless ($opt_m || $opt_x);
		    				printf CAMBEST "DNA_homol\t\"%s\"\t\"BLAT_mRNA_BEST\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$found,$score,$virtualstart,$virtualend,$eststart,$estend    if ($opt_m);
						} 
		    		}

					#############################
					# produce confirmed introns #
					#############################

					if ($opt_i) {
						my ($n) = ($virtual =~ /\S+_(\d+)$/);
						for (my $y = 1; $y < @{$entry->{'exons'}}; $y++) {
							my $last   = $y-1;
							my $first  = (${$entry->{'exons'}}[$last][1] + 1) + (($n-1)*100000);
							my $second = (${$entry->{'exons'}}[$y][0] -1) + (($n-1)*100000);
							$dir{$found} = 5 if ($opt_m);
							if (${$entry->{'exons'}}[0][2] < ${$entry->{'exons'}}[0][3]) {
								if ((${$entry->{'exons'}}[$y][2] == ${$entry->{'exons'}}[$last][3] + 1) && (($second - $first) > 2)) {
									if (exists $dir{$found} && $dir{$found} eq '3') {
										push @{$ci{$superlink}}, [$second,$first];
									}
									elsif (exists $dir{$found} && $dir{$found} eq '5') {
										push @{$ci{$superlink}}, [$first,$second];
									}
									else {
										print LOG "WARNING: Direction not found for $found\n\n";
									}
								}
							}
							elsif (${$entry->{'exons'}}[0][2] > ${$entry->{'exons'}}[0][3]) {
								if ((${$entry->{'exons'}}[$last][3] == ${$entry->{'exons'}}[$y][2] + 1) && (($second - $first) > 2)) {
									if (exists $dir{$found} && $dir{$found} eq '3') {
										push @{$ci{$superlink}}, [$first,$second];
									}
									elsif (exists $dir{$found} && $dir{$found} eq '5') {
										push @{$ci{$superlink}}, [$second,$first];
									}
									else {
										print LOG "WARNING: Direction not found for $found\n\n";
									}	
								}
							}
						}
					}
				}
			}	
    	}
	}
}
close AUTBEST;
close CAMBEST if $opt_c;

########################################################
# produce final BLAT output (including BEST and OTHER) #
########################################################

if ($opt_m) {
	open(AOTHER,"$dir/autoace.mRNA.ace");
	open(COTHER,"$dir/camace.mRNA.ace") if ($opt_c);
	open(ABEST,"$dir/autoace.best.mRNA.ace");
	open(CBEST,"$dir/camace.best.mRNA.ace") if ($opt_c);
	open(AOUT,">$dir/autoace.blat.mRNA.ace");
	open(COUT,">$dir/camace.blat.mRNA.ace") if ($opt_c);
}
elsif (!$opt_x) {	
	open(AOTHER,"$dir/autoace.EST.ace");
	open(COTHER,"$dir/camace.EST.ace") if ($opt_c);
	open(ABEST,"$dir/autoace.best.EST.ace");
	open(CBEST,"$dir/camace.best.EST.ace") if ($opt_c);
	open(AOUT,">$dir/autoace.blat.EST.ace");
	open(COUT,">$dir/camace.blat.EST.ace") if ($opt_c);
}

my (%line);
my $temp = $/;
$/ = "";
while (<ABEST>) {
#	print $_;
	if ($_ =~ /^Homol_data/) {
		$line{$_} = 1;
		print AOUT $_;
	}
}

while (<AOTHER>) {
#	print $_;
	if ($_ =~ /^Homol_data/) {
		my $line = $_;
		s/BLAT_EST_OTHER/BLAT_EST_BEST/g unless ($opt_m || $opt_x);
		s/BLAT_mRNA_OTHER/BLAT_mRNA_BEST/g   if ($opt_m);
		unless (exists $line{$_}) {
			print AOUT $line;
		}	
	}
}

while (<CBEST>) {
#	print $_;
	if ($_ =~ /^Homol_data/) {
		$line{$_} = 1;
		print COUT $_;
	}
}

while (<COTHER>) {
#	print $_;
	if ($_ =~ /^Homol_data/) {
		my $line = $_;
		s/BLAT_EST_OTHER/BLAT_EST_BEST/g unless ($opt_m || $opt_x);
		s/BLAT_mRNA_OTHER/BLAT_mRNA_BEST/g   if ($opt_m);
		unless (exists $line{$_}) {
			print COUT $line;
		}	
	}
}
$/= $temp;

###################################
# produce confirmed intron output #
###################################

if ($opt_i) {
	if ($opt_m) {
		open(ACI,">$dir/autoace.ci.mRNA.ace");
		open(CCI,">$dir/camace.ci.mRNA.ace") if ($opt_c);
	}
	elsif (!$opt_x) {
		open(ACI,">$dir/autoace.ci.EST.ace");
		open(CCI,">$dir/camace.ci.EST.ace") if ($opt_c);
	}
	foreach my $superlink (sort keys %ci) {
		my %double;
		print ACI "Sequence : \"$superlink\"\n";
		print CCI "Sequence : \"$superlink\"\n" if (($opt_c) && (exists $camace{$superlink}));
		for (my $i = 0; $i < @{$ci{$superlink}}; $i++) {
			my $merge = $ci{$superlink}->[$i][0].":".$ci{$superlink}->[$i][1];
			if (!exists $double{$merge}) {
				if ($opt_m) {
					printf ACI "Confirmed_intron %d %d mRNA\n", $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1];
					printf CCI "Confirmed_intron %d %d mRNA\n", $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1] if (($opt_c) && (exists $camace{$superlink}));
				}
				else {
					printf ACI "Confirmed_intron %d %d EST\n", $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1];
					printf CCI "Confirmed_intron %d %d EST\n", $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1] if (($opt_c) && (exists $camace{$superlink}));
				}
				$double{$merge} = 1;
			}
		}
		print ACI "\n";
		print CCI "\n" if (($opt_c) && (exists $camace{$superlink}));
	}
}


###################################

sub printhelp {
	exec('perldoc',$0);
	exit;	
}

###################################

__END__

=pod

=head1 NAME - blat2ace.pl

=head2 USAGE

blat2ace.pl maps blat output to acedb. Thereby, it produces output for autoace and camace
(autoace.ace and camace,ace). In addition, it produces files assigning the ESTs to one place 
in the genome (autoace.best.ace and camace.best.ace). ESTs that have more than one best 
match are reported in morethan1match.txt. 

blat2ace.pl  arguments:

=over 4

=item 

-a => produce output for autoace (autoace.blat.ace, /helpfiles/autoace.best.ace, /helpfiles/autoace.ace)

=item 

-c => produce output for camace (camace.blat.ace, /helpfiles/camace.best.ace, /helpfiles/camace.ace)

=item 

-i => produce output for confirmed introns (autoace.ci.ace, camace.ci.ace)

=item

-m => perform everything for mRNAs (default is EST)

=back

=cut

