#!/usr/local/bin/perl5.6.0 
#
# blat2ace.pl
# kj2
#
# Exporter to map blat data to acedb and to find the best match for each EST
#
# -i  : get confirmed introns
#
# -c  : get best matches for camace
# -s  : get best matches for stlace
#
# -e  : create output for ESTs 
# -m  : create output for mRNAs 
# -x  : create output for parasitic nematode ESTs (blatx)
# -o  : create output for other CDS
#
# -h  : print help
#
# 010905 by Kerstin Jekosch

# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2002-11-13 17:48:00 $


use strict;
use Data::Dumper;
use Ace;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Getopt::Std;
use vars qw($opt_i $opt_h $opt_s $opt_c $opt_e $opt_m $opt_o $opt_x);
$| = 1;

#############################
# variables and directories #
#############################
 
my $dir	      = "/wormsrv2/autoace/BLAT";
my $dbdir     = "/wormsrv2/autoace";
my $tace      = &tace." /wormsrv2/autoace";

my %EST_name;    # EST accession => name
my %EST_dir;     # EST accession => orientation [5|3]

my %hash;
my (%best,%other,%bestclone,%match,%ci);

our %camace;
our %stlace;

our $type = "";
our $db   = "";
our %word = (
	     EST      => 'BLAT_EST',
	     mRNA     => 'BLAT_mRNA',
	     EMBL     => 'BLAT_EMBL',
	     NEMATODE => 'BLATX_NEMATODE',
	     );

# create log file
my $rundate    = `date +%y%m%d`; 
chomp $rundate;
open (LOG, ">/wormsrv2/logs/blat2ace.log.$rundate.$$") || die "Couldn't write to the blinking log file, innit?\n";

 ########################################
 # command-line options & ramifications #
 ########################################

getopts ('csemxoih');

# Help pod documentation
&usage(0) if ($opt_h);

# Exit if no data type choosen [EST|mRNA|EMBL|NEMATODE]
&usage(1) unless ($opt_e || $opt_m || $opt_o || $opt_x); 

# Exit if multiple data types choosen [EST|mRNA|EMBL|NEMATODE]
&usage(2) if (($opt_e + $opt_m + $opt_o + $opt_x) > 1);

# assign type variable
($type = 'EST')      if ($opt_e);
($type = 'mRNA')     if ($opt_m);
($type = 'EMBL')     if ($opt_o);
($type = 'NEMATODE') if ($opt_x);

############################################
# EST data from autoace (name,orientation) #
############################################

# check to see if EST hash data exists
# make it via tablemaker queries if absent
unless (-e "/wormsrv2/autoace/BLAT/EST.dat") {
    (%EST_name,%EST_dir) = &make_EST_hash;
}
# else read it into memory
else {
    open (FH, "</wormsrv2/autoace/BLAT/EST.dat") or die "EST.dat : $!\n";
    undef $/;
    my $data = <FH>;
    eval $data;
    die if $@;
    $/ = "\n";
    close FH;
}

#########################################
# get links for database                #
#########################################

# parse links for camace
my @camclones = qw(cTel3X cTel4X cTel7X cTel33B cTel54X LINK_6R55 LINK_cTel52S SUPERLINK_CB_I SUPERLINK_CB_II SUPERLINK_CB_IIIL SUPERLINK_CB_IIIR SUPERLINK_CB_IR SUPERLINK_CB_IV SUPERLINK_CB_V SUPERLINK_CB_X); 
foreach my $camclone (@camclones) {
    $camace{$camclone} = 1;
}

# parse links for stlace
my @stlclones = qw(SUPERLINK_RW1 SUPERLINK_RW1R SUPERLINK_RW2 SUPERLINK_RW3A SUPERLINK_RW3B SUPERLINK_RW4 SUPERLINK_RW5 SUPERLINK_RWXL SUPERLINK_RWXR);
foreach my $stlclone (@stlclones) {
    $stlace{$stlclone} = 1;
}

############################
# map the blat hits to ace #
############################

print LOG "Start mapping\n\n";

# output filehandle
open(ACE,  ">$dir/autoace.$type.ace")  or die "Cannot open $dir/autoace.${type}.ace $!\n";

# input filehandle
open(BLAT, "<$dir/${type}_out.psl")  or die "Cannot open $dir/${type}_out.psl $!\n";
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
    if (($opt_e)  && (exists $EST_name{$est})) {
	my $estname  = $EST_name{$est};
	if ($est ne $estname) {
	    $est = $estname;
	    print LOG "EST name '$est' was replaced by '$estname'\n\n";
	}
    }
    my @lengths     = split (/,/, $f[18]);
    my @eststarts   = split (/,/, $f[19]);
    my @slinkstarts = split (/,/, $f[20]);


#    print "$est maps to $superlink [currentDB: $db => $camace{$superlink} | $stlace{$superlink}]\n";
#    # next if LINK is part of camace BUT we want stlace
#    next if ((defined ($camace{$superlink})) && ($opt_s));
#    
#    # next if LINK is part of stlace BUT we want camace
#    next if ((defined ($stlace{$superlink})) && ($opt_c));
#    print "$est will be processed\n";


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
	$virtual = "$word{$type}:${superlink}_${startvirtual}";
    }	
    elsif (($startvirtual == ($endvirtual - 1)) && (($matchend%100000) <= 50000)) {
	$virtual = "$word{$type}:${superlink}_${startvirtual}";
    }
    else {
	print LOG "$est wasn't assigned to a virtual object as match size was too big\n";
	print LOG "Start is $matchstart, end is $matchend on $superlink\n\n";
	next;
    }
    
    ###################
    # calculate score #
    ###################

    my $sum   = 0;
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
	
	# blatx 6-frame translation v 6-frame translation
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

	# write to output file
	print  ACE "Homol_data : \"$virtual\"\n";
	if ($type eq "NEMATODE") {
	    printf ACE "DNA_homol\t\"%s\"\t\"$word{$type}\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$est,$score,$virtualstart,$virtualend,$eststart,$estend;
	}
	else {
	    printf ACE "DNA_homol\t\"%s\"\t\"$word{$type}_OTHER\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$est,$score,$virtualstart,$virtualend,$eststart,$estend;
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

#########################################
# 
#########################################


####################################
# produce outfile for best matches #
####################################

&usage(20) if ($opt_x);

open (AUTBEST, ">$dir/autoace.best.$type.ace");
open (STLBEST, ">$dir/stlace.best.$type.ace");
open (CAMBEST, ">$dir/camace.best.$type.ace");

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
		    
		    # print line for autoace
		    print  AUTBEST "Homol_data : \"$virtual\"\n";
		    printf AUTBEST "DNA_homol\t\"%s\"\t\"$word{$type}_BEST\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$found,$score,$virtualstart,$virtualend,$eststart,$estend;
		    # camace
		    if ($camace{$superlink}) {
			print  CAMBEST "Homol_data : \"$virtual\"\n";
			printf CAMBEST "DNA_homol\t\"%s\"\t\"$word{$type}_BEST\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$found,$score,$virtualstart,$virtualend,$eststart,$estend;
		    }
		    # and stlace
		    elsif ($stlace{$superlink}) {
			print  STLBEST "Homol_data : \"$virtual\"\n";
			printf STLBEST "DNA_homol\t\"%s\"\t\"$word{$type}_BEST\"\t%.1f\t%d\t%d\t%d\t%d\n\n",$found,$score,$virtualstart,$virtualend,$eststart,$estend;
		    }
		    
		}
		    
		#############################
		# produce confirmed introns #
		#############################
		
		if ($opt_i) {
		    my ($n) = ($virtual =~ /\S+_(\d+)$/);
		    for (my $y = 1; $y < @{$entry->{'exons'}}; $y++) {
			my $last   = $y - 1;
			my $first  =  ${$entry->{"exons"}}[$last][1] + 1 + (($n-1)*100000);
		        my $second = (${$entry->{'exons'}}[$y][0]   - 1) + (($n-1)*100000);
		        $EST_dir{$found} = 5 if ($opt_m || $opt_o);
		        if (${$entry->{'exons'}}[0][2] < ${$entry->{'exons'}}[0][3]) {
	    if ((${$entry->{'exons'}}[$y][2] == ${$entry->{'exons'}}[$last][3] + 1) && (($second - $first) > 2)) {
	if (exists $EST_dir{$found} && $EST_dir{$found} eq '3') {
					push @{$ci{$superlink}}, [$second,$first];
				    }
				    elsif (exists $EST_dir{$found} && $EST_dir{$found} eq '5') {
					push @{$ci{$superlink}}, [$first,$second];
				    }
				    else {
					print LOG "WARNING: Direction not found for $found\n\n";
				    }
				}
                            }
                            elsif (${$entry->{'exons'}}[0][2] > ${$entry->{'exons'}}[0][3]) {
                                if ((${$entry->{'exons'}}[$last][3] == ${$entry->{'exons'}}[$y][2] + 1) && (($second - $first) > 2)) {
		                    if (exists $EST_dir{$found} && $EST_dir{$found} eq '3') {
                                        push @{$ci{$superlink}}, [$first,$second];
                                    }
                                    elsif (exists $EST_dir{$found} && $EST_dir{$found} eq '5') {
                                        push @{$ci{$superlink}}, [$second,$first]; 
                                    }
                                    else {
                                        print LOG "WARNING: Direction not found for $found\n\n";
                                    }
                                }
                            }
                        }
                    }
                    # this section 'produce confirmed introns'


            }	
        }
    }
}
close AUTBEST;
close CAMBEST;
close STLBEST;

########################################################
# produce final BLAT output (including BEST and OTHER) #
########################################################

&usage(20) if ($opt_x);

# autoace
open (OUT_autoace, ">$dir/autoace.blat.$type.ace") or die "$!";
# camace
open (OUT_camace,  ">$dir/camace.blat.$type.ace")  or die "$!";
# stlace
open (OUT_stlace,  ">$dir/stlace.blat.$type.ace")  or die "$!";


#open(AOUT,  ">$dir/autoace.blat.$type.ace");

my (%line);
my $temp = $/;
$/ = "";

my $superlink = "";

# assign 
open(ABEST,  "<$dir/autoace.best.$type.ace");
while (<ABEST>) {
#   print $_;
    if ($_ =~ /^Homol_data/) {
	$line{$_} = 1;
	($superlink) = (/\"BLAT\_$type\:(\S+)\_\d+\"/);

#	Homol_data : "BLAT_EST:SUPERLINK_RW5_45"

	print OUT_autoace "// Source $superlink\n\n";
	print OUT_autoace $_;
	    
	# camace
	if ($camace{$superlink}) {
	    print OUT_camace $_;
	}
	# and stlace
	elsif ($stlace{$superlink}) {
	    print OUT_stlace $_;
	}
    }
}
close ABEST;


open(AOTHER, "<$dir/autoace.$type.ace");
while (<AOTHER>) {
#	print $_;
    if ($_ =~ /^Homol_data/) {
	my $line = $_;
	s/BLAT_EST_OTHER/BLAT_EST_BEST/g unless ($opt_m || $opt_o || $opt_x);
	s/BLAT_mRNA_OTHER/BLAT_mRNA_BEST/g   if ($opt_m);
	s/BLAT_EMBL_OTHER/BLAT_EMBL_BEST/g   if ($opt_o);

	unless (exists $line{$_}) {
	    print OUT_autoace $line;

	    # camace
	    if ($camace{$superlink}) {
		print OUT_camace $line;
	    }
	    # and stlace
	    elsif ($stlace{$superlink}) {
		print OUT_stlace $line;
	    }
	    
	}	
    }
}
close AOTHER;

$/= $temp;

###################################
# produce confirmed intron output #
###################################

if ($opt_i) {

    open(CI_auto, ">$dir/autoace.ci.${type}.ace");
    open(CI_cam,  ">$dir/camace.ci.${type}.ace");
    open(CI_stl,  ">$dir/stlace.ci.${type}.ace");
   
    foreach my $superlink (sort keys %ci) {
	my %double;
	
	print CI_auto "\nSequence : \"$superlink\"\n";
	print CI_stl  "\nSequence : \"$superlink\"\n" if ($stlace{$superlink});
	print CI_cam  "\nSequence : \"$superlink\"\n" if ($camace{$superlink});
	
	for (my $i = 0; $i < @{$ci{$superlink}}; $i++) {
	    my $merge = $ci{$superlink}->[$i][0].":".$ci{$superlink}->[$i][1];
	    if (!exists $double{$merge}) {
		if ($opt_m) {
		    printf CI_auto "Confirmed_intron %d %d mRNA\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1];
		    (printf CI_cam "Confirmed_intron %d %d mRNA\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1]) if ($camace{$superlink});
		    (printf CI_stl "Confirmed_intron %d %d mRNA\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1]) if ($stlace{$superlink});
		}
		if ($opt_o) {
		    printf CI_auto "Confirmed_intron %d %d Homol\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1];
		    (printf CI_cam "Confirmed_intron %d %d Homol\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1]) if ($camace{$superlink});
		    (printf CI_stl "Confirmed_intron %d %d Homol\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1]) if ($stlace{$superlink});
		}
		if ($opt_e) {
		    printf CI_auto "Confirmed_intron %d %d EST\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1];
		    (printf CI_cam "Confirmed_intron %d %d EST\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1]) if ($camace{$superlink});
		    (printf CI_stl "Confirmed_intron %d %d EST\n",  $ci{$superlink}->[$i][0], $ci{$superlink}->[$i][1]) if ($stlace{$superlink});
		}
		$double{$merge} = 1;
	    }
	}
    }
    
    close CI_auto;
    close CI_cam;
    close CI_stl;

}

##############################
# hasta luego                #
##############################

exit(0);

#################################################################################
### Subroutines                                                               ###
#################################################################################

#########################################
# get EST names  (-e option only)       #
#########################################

sub commands {
    
my $command1=<<EOF;
Table-maker -p "/wormsrv2/autoace/wquery/ESTacc2names.def"
quit
EOF

my $command2=<<EOF;
Table-maker -p "/wormsrv2/autoace/wquery/ESTorient.def"
quit
EOF

    return($command1,$command2);

}

sub make_EST_hash {
    
    my ($command1,$command2) = &commands;
    my ($acc,$name,$orient);

    my %EST_name = ();
    my %EST_dir  = ();

    # get EST names  (-e option only)       #
    open (TACE, "echo '$command1' | $tace | ");
    while (<TACE>) {
	chomp;
	next if ($_ eq "");
	next if (/\/\//);
	s/acedb\>\s//g;
    	s/\"//g;
	s/EMBL://g;
	($acc,$name) = ($_ =~ /^(\S+)\s(\S+)/);
	$name = $acc unless ($name);
	$EST_name{$acc} = $name;
    }
    close TACE;

    # get EST orientation (5' or 3')    #
    open (TACE, "echo '$command2' | $tace | ");
    while (<TACE>) {
	chomp;
	next if ($_ eq "");
	next if (/\/\//);
	s/acedb\>\s//g;
	s/\"//g;
	($name,$orient) = ($_ =~ /^(\S+)\s+EST_(\d)/);
	$EST_dir{$name} = $orient if ($orient);
    }
    close TACE;

    # Data::Dumper write hash to /wormsrv2/autoace/BLAT/EST.dat
    open (OUT, ">/wormsrv2/autoace/BLAT/EST.dat") or die "EST.dat : $!";
    print OUT Data::Dumper->Dump([\%EST_name],['*EST_name']);
    print OUT Data::Dumper->Dump([\%EST_dir],['*EST_dir']);
    close OUT;

    return (%EST_name,%EST_dir);


}





###################################

sub usage {
    my $error = shift;
    
    if ($error == 1) {
	# No data-type choosen
	print "\nNo data option choosen [-e|m|o|x]\n";
	print "Run with one of the above options\n\n";
	exit(0);
    }
    if ($error == 2) {
	# 'Multiple data-types choosen
	print "\nMultiple data option choosen [-e|m|o|x]\n";
	print "Run with one of the above options\n\n";
	exit(0);
    }
    if ($error == 3) {
	# 'chromosome.ace' file is not there or unreadable
	print "\nThe WormBase 'chromosome.ace' file you specified does not exist or is non-readable.\n";
	print "Check File: ''\n\n";
	exit(0);
    }
    if ($error == 20) {
	# 
	print "\Don't want to do this for the -x option.\n";
	print "hasta luego\n\n";
	exit(0);
    }
    elsif ($error == 0) {
	# Normal help menu
	exec ('perldoc',$0);
    }
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

=head1 AUTHOR

Kerstin Jekosch (kj2@sanger.ac.uk)

=cut

