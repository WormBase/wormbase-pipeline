#!/usr/local/bin/perl5.6.0 -w
#
# atcg.pl
#
# Nucleotide and dinucleotide composition script
#
# Usage: atcg.pl - options
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2002-08-01 10:36:36 $
#
##############################

use strict;
use lib "/wormsrv2/scripts/";

use Getopt::Std;
our $opt_i = "";
our $opt_c = "";
getopts ('i:c');

print "Specify FASTA format input file using -i option\n" if (!defined($opt_i));

print "\n";
print "atcg.pl - a program to calculate nucleotide and dinucleotide composition\n";
print "Can analyse multiple sequences, but will concatenate the results\n";
print "together. Keith Bradnam, 1998\n\n";


# Process sequence

my $seq = &read_input_file($opt_i);
my $total_length=length($seq);       # this part to calculate GC content
print "Length of sequence is $total_length bp:\n\n";

my($A,$T,$C,$G,$AA,$AT,$AC,$AG,$TA,$TT,$TC,$TG,$CA,$CT,$CC,$CG,$GA,$GT,$GC,$GG);
$A=$T=$C=$G=$AA=$AT=$AC=$AG=$TA=$TT=$TC=$TG=$CA=$CT=$CC=$CG=$GA=$GT=$GC=$GG=0;


# mononucleotides
&calculate_mononucleotides($seq);

# dinucleotides
&calculate_dinucleotides($seq);

exit(0);


#######################
# Read input file
#######################

sub read_input_file{
	my $input = shift;
	open(IN,"$input") || die "Couldn't open input file.\n\n";
	my $seq;
	my $temp=<IN>;     	
	while($temp=<IN>){
	    if($temp !~ />/){
        	chomp($temp);
	        $temp =~  tr/a-z/A-Z/;
	        $seq.=$temp;
	    }
	}
	close(IN);
	return($seq);
}

####################################
# Calculate basic stats
###################################
sub calculate_mononucleotides{
	my $seq =shift;
	my $A = $seq =~ tr/A/A/; 
	my $T = $seq =~ tr/T/T/;
	my $C = $seq =~ tr/C/C/;
	my $G = $seq =~ tr/G/G/;

	my $p_A = sprintf "%.2f",(($A/$total_length)*100);
	my $p_T = sprintf "%.2f",(($T/$total_length)*100);
	my $p_C = sprintf "%.2f",(($C/$total_length)*100);
	my $p_G = sprintf "%.2f",(($G/$total_length)*100);

	$A = sprintf "%10d",$A;
	$T = sprintf "%10d",$T;
	$C = sprintf "%10d",$C;
	$G = sprintf "%10d",$G;

	print "Nucleotide count:\n";
	print "A = $A ($p_A%)\n";
	print "T = $T ($p_T%)\n";
	print "C = $C ($p_C%)\n";
	print "G = $G ($p_G%)\n\n";

}


sub calculate_dinucleotides{
	my $seq =shift;
	my $counter =0;
	foreach my $i(0..length($seq)){
	    $counter++;
	    my $tmp = substr($seq,$i,2);   
	  SWITCH: {                    
	      if ($tmp =~ /AA/){$AA++; last SWITCH;}
	      if ($tmp =~ /AT/){$AT++; last SWITCH;}
	      if ($tmp =~ /AC/){$AC++; last SWITCH;}
	      if ($tmp =~ /AG/){$AG++; last SWITCH;}
	      if ($tmp =~ /TA/){$TA++; last SWITCH;}
	      if ($tmp =~ /TT/){$TT++; last SWITCH;}
	      if ($tmp =~ /TC/){$TC++; last SWITCH;}
	      if ($tmp =~ /TG/){$TG++; last SWITCH;}
	      if ($tmp =~ /CA/){$CA++; last SWITCH;}
	      if ($tmp =~ /CT/){$CT++; last SWITCH;}
	      if ($tmp =~ /CC/){$CC++; last SWITCH;}
	      if ($tmp =~ /CG/){$CG++; last SWITCH;}
	      if ($tmp =~ /GA/){$GA++; last SWITCH;}
	      if ($tmp =~ /GT/){$GT++; last SWITCH;}
	      if ($tmp =~ /GC/){$GC++; last SWITCH;}
	      if ($tmp =~ /GG/){$GG++; last SWITCH;}
	  }
	}

	my $p_AA = sprintf "%.2f",(($AA/$counter)*100);
	my $p_AT = sprintf "%.2f",(($AT/$counter)*100);
	my $p_AC = sprintf "%.2f",(($AC/$counter)*100);
	my $p_AG = sprintf "%.2f",(($AG/$counter)*100);
	my $p_TA = sprintf "%.2f",(($TA/$counter)*100);
	my $p_TT = sprintf "%.2f",(($TT/$counter)*100);
	my $p_TC = sprintf "%.2f",(($TC/$counter)*100);
	my $p_TG = sprintf "%.2f",(($TG/$counter)*100);
	my $p_CA = sprintf "%.2f",(($CA/$counter)*100);
	my $p_CT = sprintf "%.2f",(($CT/$counter)*100);
	my $p_CC = sprintf "%.2f",(($CC/$counter)*100);
	my $p_CG = sprintf "%.2f",(($CG/$counter)*100);
	my $p_GA = sprintf "%.2f",(($GA/$counter)*100);
	my $p_GT = sprintf "%.2f",(($GT/$counter)*100);
	my $p_GC = sprintf "%.2f",(($GC/$counter)*100);
	my $p_GG = sprintf "%.2f",(($GG/$counter)*100);

	$AA = sprintf "%10d",$AA;
	$AT = sprintf "%10d",$AT;
	$AC = sprintf "%10d",$AC;
	$AG = sprintf "%10d",$AG;
	$TA = sprintf "%10d",$TA;
	$TT = sprintf "%10d",$TT;
	$TC = sprintf "%10d",$TC;
	$TG = sprintf "%10d",$TG;
	$CA = sprintf "%10d",$CA;
	$CT = sprintf "%10d",$CT;
	$CC = sprintf "%10d",$CC;
	$CG = sprintf "%10d",$CG;
	$GA = sprintf "%10d",$GA;
	$GT = sprintf "%10d",$GT;
	$GC = sprintf "%10d",$GC;
	$GG = sprintf "%10d",$GG;



	print "Dinucleotide count:\n";
	print "AA = $AA ($p_AA%)\n";
	print "AT = $AT ($p_AT%)\n";
	print "AC = $AC ($p_AC%)\n";
	print "AG = $AG ($p_AG%)\n";
	print "TA = $TA ($p_TA%)\n";
	print "TT = $TT ($p_TT%)\n";
	print "TC = $TC ($p_TC%)\n";
	print "TG = $TG ($p_TG%)\n";
	print "CA = $CA ($p_CA%)\n";
	print "CT = $CT ($p_CT%)\n";
	print "CC = $CC ($p_CC%)\n";
	print "CG = $CG ($p_CG%)\n";
	print "GA = $GA ($p_GA%)\n";
	print "GT = $GT ($p_GT%)\n";
	print "GC = $GC ($p_GC%)\n";
	print "GG = $GG ($p_GG%)\n\n";



	if ($opt_c ne ""){

		my $tot_obs = $counter;

		my $AA_exp = (($AA+$AT+$AC+$AG)*($AA+$TA+$CA+$GA))/$tot_obs;
		my $AT_exp = (($AA+$AT+$AC+$AG)*($TA+$TT+$TC+$TG))/$tot_obs;
		my $AC_exp = (($AA+$AT+$AC+$AG)*($CA+$CT+$CC+$CG))/$tot_obs;
		my $AG_exp = (($AA+$AT+$AC+$AG)*($GA+$GT+$GC+$GG))/$tot_obs;

		my $TA_exp = (($TA+$TT+$TC+$TG)*($AA+$AT+$AC+$AG))/$tot_obs;
		my $TT_exp = (($TA+$TT+$TC+$TG)*($TA+$TT+$TC+$TG))/$tot_obs;
		my $TC_exp = (($TA+$TT+$TC+$TG)*($CA+$CT+$CC+$CG))/$tot_obs;
		my $TG_exp = (($TA+$TT+$TC+$TG)*($GA+$GT+$GC+$GG))/$tot_obs;

		my $CA_exp = (($CA+$CT+$CC+$CG)*($AA+$AT+$AC+$AG))/$tot_obs;
		my $CT_exp = (($CA+$CT+$CC+$CG)*($TA+$TT+$TC+$TG))/$tot_obs;
		my $CC_exp = (($CA+$CT+$CC+$CG)*($CA+$CT+$CC+$CG))/$tot_obs;
		my $CG_exp = (($CA+$CT+$CC+$CG)*($GA+$GT+$GC+$GG))/$tot_obs;

		my $GA_exp = (($GA+$GT+$GC+$GG)*($AA+$AT+$AC+$AG))/$tot_obs;
		my $GT_exp = (($GA+$GT+$GC+$GG)*($TA+$TT+$TC+$TG))/$tot_obs;
		my $GC_exp = (($GA+$GT+$GC+$GG)*($CA+$CT+$CC+$CG))/$tot_obs;
		my $GG_exp = (($GA+$GT+$GC+$GG)*($GA+$GT+$GC+$GG))/$tot_obs;


		my $chi=((($AA-$AA_exp)**2)/$AA_exp)+((($AT-$AT_exp)**2)/$AT_exp)+
		    ((($AC-$AC_exp)**2)/$AC_exp)+((($AG-$AG_exp)**2)/$AG_exp)+
		    ((($TA-$TA_exp)**2)/$TA_exp)+((($TT-$TT_exp)**2)/$TT_exp)+
		    ((($TC-$TC_exp)**2)/$TC_exp)+((($TG-$TG_exp)**2)/$TG_exp)+
		    ((($CA-$CA_exp)**2)/$CA_exp)+((($CT-$CT_exp)**2)/$CT_exp)+
		    ((($CC-$CC_exp)**2)/$CC_exp)+((($CG-$CG_exp)**2)/$CG_exp)+
		    ((($GA-$GA_exp)**2)/$GA_exp)+((($GT-$GT_exp)**2)/$GT_exp)+
		    ((($GC-$GC_exp)**2)/$GC_exp)+((($GG-$GG_exp)**2)/$GG_exp);

		printf  "Chi squared value = %6.2f\n", $chi;
	                                
		print "Significance level at 5% = 16.92\n";
		print "Significance level at 1% = 21.67\n";
	}

}

