#!/usr/local/bin/perl5.6.0
#
# confirm_introns.pl
# kj2
#
# evaluates confirmed introns
#
# confirm_introns.pl [-c -a -v -o -m -h] 
#
# -a = ""; # produce output for autoace
# -c = ""; # produce output for camace
# -m = ""; # perform everything for mRNAs (default is EST)
# -h = ""; # help
# -o = ""; # perform everthing for EMBL other CDS  (default is EST)
# -v = ""; # be verbose
# 5.10.01 Kerstin Jekosch

use strict;
use Getopt::Std;
use vars qw($opt_c $opt_s $opt_a $opt_h $opt_v $opt_o $opt_m $opt_e);
$|=1;

#######################
# variables and files #
#######################

my $dir = '/wormsrv2/autoace/BLAT';
my (%seq,%intron,$temp);
 
our $data_type = "";
our $db = "";
our %word = (
	     EMBL => 'BLAT_EST',
	     mRNA => 'BLAT_mRNA',
	     EMBL => 'BLAT_EMBL',
	     );

########################################
# command-line options & ramifications #
########################################

getopts('csemohv');

# Help pod documentation
&usage(0) if ($opt_h);

# Exit if no data type choosen [EST|mRNA|EMBL]
&usage(1) unless ($opt_e || $opt_m || $opt_o); 

# Exit if multiple data types choosen [EST|mRNA|EMBL]
&usage(2) if (($opt_e + $opt_m + $opt_o) > 1);

# Exit if the autoace.fa file is absent
&usage(3) unless ("-e $dir/autoace.fa");

# assign data type variable
($data_type = 'EST')  if ($opt_e);
($data_type = 'mRNA') if ($opt_m);
($data_type = 'EMBL') if ($opt_o);

# if no database option given then do both.
# (i know but oranges are not the only fruit)
($db = "camace")  if ($opt_c);
($db = "stlace")  if ($opt_s);
($db = "autoace") if (!$opt_c && !$opt_s);

########################################
# opening output files                 #
########################################

open(GOOD,">$dir/$db.good_introns.$data_type.ace") or die "Cannot open $dir/$db.good_introns.$data_type.ace $!\n";
open(BAD, ">$dir/$db.bad_introns.$data_type.ace")  or die "Cannot open $dir/$db.bad_introns.$data_type.ace $!\n"; 


#########################################
# loop over file with confirmed introns #
#########################################

my $lala = $/;
$/ = "";
open(CI,  "<$dir/$db.ci.$data_type.ace")           or die "Cannot open $dir/$db.ci.$data_type.ace $!\n";
while (<CI>) {
    next unless /^\S/;
    if (/Sequence : \"(\S+)\"/) {
	my $link = $1;
	print "Sequence : $link\n" if ($opt_v);
	my @introns = split /\n/, $_;
	
	#########################
	# get the link sequence #
	#########################
	
	open(SEQ,"$dir/autoace.fa") || &usage(3);
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
		# map to S_child #
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
		    $virtual = "Confirmed_intron_EST:" .$link."_".$startvirt unless ($opt_m);
		    $virtual = "Confirmed_intron_mRNA:".$link."_".$startvirt     if ($opt_m);
		}
		elsif (($startvirt == ($endvirt - 1)) && (($last%100000) <= 50000)) {
		    $virtual = "Confirmed_intron_EST:" .$link."_".$startvirt unless ($opt_m);
		    $virtual = "Confirmed_intron_mRNA:".$link."_".$startvirt     if ($opt_m);
		}
		
		#################
		# check introns #
		#################
		
		my $firstcalc = int($f[1]/100000);
		my $seccalc   = int($f[2]/100000);
		print STDERR "Problem with $test\n" unless (defined $firstcalc && defined $seccalc); 
		my ($one,$two);
		if ($firstcalc == $seccalc) {
		    $one = $f[1]%100000;
		    $two = $f[2]%100000;
		}
		elsif ($firstcalc == ($seccalc-1)) {
		    $one = $f[1]%100000;
		    $two = $f[2]%100000 + 100000;
		    print STDERR "$virtual: $one $two\n";
		}
		elsif (($firstcalc-1) == $seccalc) {
		    $one = $f[1]%100000 + 100000;
		    $two = $f[2]%100000;
		    print STDERR "$virtual: $one $two\n";
		} 
		print STDERR "Problem with $test\n" unless (defined $one && defined $two); 
		
		if ( ( (($start eq 'gt') || ($start eq 'gc')) && ($end eq 'ag')) ||
		     (  ($start eq 'ct') && (($end eq 'ac') || ($end eq 'gc')) ) ) {	 
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
close CI;

##############################
# hasta luego                #
##############################

exit(0);

#################################################################################
### Subroutines                                                               ###
#################################################################################

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
	# 'autoace.fa' file is not there or unreadable
	print "\nThe WormBase 'autoace.fa' file you does not exist or is non-readable.\n";
	print "Check File: '/wormsrv2/autoace/BLAT/autoace.fa'\n\n";
	exit(0);
    }
    elsif ($error == 0) {
	# Normal help menu
	exec ('perldoc',$0);
    }
}

__END__

=pod

=head2 NAME - confirm_introns.pl

=head1 USAGE 

=over 4

=item confirm_introns.pl [-options]

=back

confirm_introns.pl checks the introns from the BLAT mappings for gt...ag consistency.

confirm_introns.pl mandatory arguments:

=over 4

=item -e, create virtual objects for EST mappings

=item -m, create virtual objects for mRNA mappings

=item -o, create virtual objects for other EMBL CDS mappings

=item -x, create virtual objects for Parasitic nematode consensus cluster mappings

=item -f <file>, use alternative chromosome.ace file 
=back

confirm_introns.pl  OPTIONAL arguments:

=over 4

=item -c, create virtual objects for camace only

=item -s, create virtual objects for stlace only

=item -h, Help pages

=back

Examples:

confirm_introns.pl -e
Uses the chromosome.ace dump file to generate virtual objects for EST mappings.
This covers all Subsequence children of each CHROMOSOME in the current WS release.

confirm_introns.pl -cm
Uses the chromosome.ace dump file to generate virtual objects for mRNA mappings.
This covers all Subsequence children of each CHROMOSOME in camace.

=head1 AUTHOR

Kerstin Jekosch (kj2@sanger.ac.uk)

=cut
