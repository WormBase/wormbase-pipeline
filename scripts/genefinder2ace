#!/usr/local/bin/perl 

# parser for genefinder 
# which has a very unfriendy format. 
#
# Steven Jones, Sanger Centre. 

if ($#ARGV !=1) {print "$0 <genfinder output> <sequence name>\n";exit;}

open(genefinder,"$ARGV[0]");

@genename=split(/\s+/,"a b c d e f g h i j k l m n o p q r s t u v w x y z aa bb cc dd ee ff gg hh jj kk ll mm nn oo pp qq rr ss tt uu vv ww xx yy zz");
$name=0;

while (<genefinder>) {

    if (/\.\./) {
	#this is the difficult bit having to work out the frames of incomplete genes
	$frame=0;$nostart="";
	if (/^\*.+U(\d)$/) {$nostart="yes";$frame=$1; "gene has no start\n";}
	if (/^U(\d)$/) {$nostart="yes";$frame=$1;print "gene has no start\n";}
	if (/\[\s+(\S+)\s+\]/) {print STDERR "$score\n";$score=$1;}
	s/\[\s+\]//g;
	#if we run out of the alphabet :-)
	if ($genename[$name] eq "") 
	{print STDERR "running  out of names!!! adding zzN names \n";$extra++;$genename[$name]="zz".$extra;}

	@line=split(/\s+/,$_);
	#print @line;
	$strand="+";
	undef @exons;$count="";
	$first="";$last="";
	foreach $exon (@line) {
		if ($exon=~/\*/) {$strand="-";}
		
		#sort of problems with the + strand for genes with no ATG
		if ($nostart eq "yes" && $strand eq "+" && $first eq "") {$exon=~/(\d+)..(\d+)/;$exon=$1+(3-$frame)."..".$2;
								      }
		if ($exon=~/(\d+\.\.\d+)$/) {$count++;$exon{$count}=$1;}
		if ($exon=~/\d+\.\.(\d+)/) {$last=$1;} 
		if ($exon=~/(\d+)\.\.\d+$/ && $first eq "") {$first=$1;}
	}

	#sort of the problem of - strand genes with no ATG
	if ($nostart eq "yes" && $strand eq "-") {$exon{$count}=~/(\d+)..(\d+)/;$exon{$count}=$1."..".($2-(3-$frame));}

	for ($i=1;$i<=$count;$i++) {push(@exons,$exon{$i});}

	print "\nSequence $ARGV[1]\n";
	if ($strand eq "+") {print "Subsequence $ARGV[1].$genename[$name] $first $last\n\n";}
	if ($strand eq "-") {print "Subsequence $ARGV[1].$genename[$name] $last $first\n\n";}
	
	print "\nSequence $ARGV[1].$genename[$name]\n"; 
	print "CDS\nMethod Genefinder\n";
	print "CDS_predicted_by Genefinder $score\n";

	$name++;

	sub numberically {$a <=> $b;}

	$i=0;undef @revexons;
	if ($strand eq "-") {foreach $exon (@exons) {$exon{$i}=$exon;$i++;}
			     for ($i--;$i>=0;$i--) {push (@revexons,$exon{$i});}
			     undef @exons;
			     @exons=@revexons;		
		}
	

	foreach $exon (@exons) {$exon=~/(\d+)\.\.(\d+)/;$begin=$1;$end=$2;
	if ($strand eq "+")  {print "Source_exons ",$begin-$first+1,"\t",$end-$first+1,"\n";}

	if ($strand eq "-")  {
	                      $exonstart=$end-$last-1;$exonend=$begin-$last-1;
			      if ($exonstart < 0) {$exonstart=$exonstart*-1;}
			      if ($exonend < 0) {$exonend=$exonend*-1;}
				print "Source_exons $exonstart\t$exonend\n";}
		}
	}
}

print "\n";


