#!/usr/local/bin/perl -w

print "\natcg - a program to calculate nucleotide and dinucleotide composition\natcg can analyse multiple sequence
s, but will concatenate the results\ntogether.\n\nKeith Bradnam, 1998\n\n";

open(IN,"$ARGV[0]") || die "Couldn't open file, please specify the FASTA sequence file you wish atcg\nto analyse o
n the command line, e.g. atcg filename.\n\n";


$temp=<IN>;     

while($temp=<IN>){
    if($temp !~ />/){
        chomp($temp);
        $temp =~  tr/a-z/A-Z/;
        $seq.=$temp;
    }
}
close(IN);
    
$total_length=length($seq);       # this part to calculate GC content
  
print "Length of sequence is $total_length bp:\n\n";

$A=$T=$C=$G=$AA=$AT=$AC=$AG=$TA=$TT=$TC=$TG=$CA=$CT=$CC=$CG=$GA=$GT=$GC=$GG=0;
$A = $seq =~ tr/A/A/; 
$T = $seq =~ tr/T/T/;
$C = $seq =~ tr/C/C/;
$G = $seq =~ tr/G/G/;

print "Nucleotide count:\n";
print "A = $A\tT = $T\tC = $C\tG = $G\n\n";


foreach $i(0..length($seq)){
    $tmp = substr($seq,$i,2);   
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
      $nothing=1;
  }
}
print "Dinucleotide count: (rows show 5' nucleotide, columns show 3' nucleotide\n";
print "\tA\tT\tC\tG\t\n";
print "\t--------------------------\n";
print "A\t$AA\t$AT\t$AC\t$AG\n";
print "\n";

print "T\t$TA\t$TT\t$TC\t$TG\n";
print "\n";

print "C\t$CA\t$CT\t$CC\t$CG\n";                                
print "\n";

print "G\t$GA\t$GT\t$GC\t$GG\n";                                #
print "\n";



$tot_obs= $AA+$AT+$AC+$AG+$TA+$TT+$TC+$TG+$CA+$CT+$CC+$CG+$GA+$GT+$GC+$GG;
$AA_exp = (($AA+$AT+$AC+$AG)*($AA+$TA+$CA+$GA))/$tot_obs;
$AT_exp = (($AA+$AT+$AC+$AG)*($TA+$TT+$TC+$TG))/$tot_obs;
$AC_exp = (($AA+$AT+$AC+$AG)*($CA+$CT+$CC+$CG))/$tot_obs;
$AG_exp = (($AA+$AT+$AC+$AG)*($GA+$GT+$GC+$GG))/$tot_obs;

$TA_exp = (($TA+$TT+$TC+$TG)*($AA+$AT+$AC+$AG))/$tot_obs;
$TT_exp = (($TA+$TT+$TC+$TG)*($TA+$TT+$TC+$TG))/$tot_obs;
$TC_exp = (($TA+$TT+$TC+$TG)*($CA+$CT+$CC+$CG))/$tot_obs;
$TG_exp = (($TA+$TT+$TC+$TG)*($GA+$GT+$GC+$GG))/$tot_obs;

$CA_exp = (($CA+$CT+$CC+$CG)*($AA+$AT+$AC+$AG))/$tot_obs;
$CT_exp = (($CA+$CT+$CC+$CG)*($TA+$TT+$TC+$TG))/$tot_obs;
$CC_exp = (($CA+$CT+$CC+$CG)*($CA+$CT+$CC+$CG))/$tot_obs;
$CG_exp = (($CA+$CT+$CC+$CG)*($GA+$GT+$GC+$GG))/$tot_obs;

$GA_exp = (($GA+$GT+$GC+$GG)*($AA+$AT+$AC+$AG))/$tot_obs;
$GT_exp = (($GA+$GT+$GC+$GG)*($TA+$TT+$TC+$TG))/$tot_obs;
$GC_exp = (($GA+$GT+$GC+$GG)*($CA+$CT+$CC+$CG))/$tot_obs;
$GG_exp = (($GA+$GT+$GC+$GG)*($GA+$GT+$GC+$GG))/$tot_obs;


$chi=((($AA-$AA_exp)**2)/$AA_exp)+((($AT-$AT_exp)**2)/$AT_exp)+
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
