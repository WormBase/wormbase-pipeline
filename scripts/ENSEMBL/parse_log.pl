#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  parse_log.pl
#
#        USAGE:  ./parse_log.pl 
#
#  DESCRIPTION:  
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:   (), <>
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  03/07/06 16:56:43 BST
#     REVISION:  ---
#===============================================================================

use strict;
use YAML;

my $chromosome;
my $gene;
my $transcript;
my %exon;
while (<>){
	if (/^Slice = chromosome:CEL160:(\w+):/){$chromosome= $1}
	elsif (/^Gene\s\d+(\S+)\stype\swormbase/){$gene = $1}
	elsif (/^Transcript\s\d+\s(\S+)/){$transcript = $1}
	elsif (/\s+Exon\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)\s+\d+\s+(\S+)/){
		my $orientation=$4==-1?'-':'+';
		$exon{$1}=[$2,$3,$orientation,$5+1];
	}
	elsif (/Short exon\s\((.*)\).*for exon\s(\d+)/){
		printf "CHROMOSOME_$chromosome\tshort_exon\texon\t%i\t%i\t.\t%s\t%i\tTranscript \"$transcript\" ($1)\n",@{$exon{$2}}
	}
	elsif (/Long exon\s\((.*)\).*for exon\s(\d+)/){
		printf "CHROMOSOME_$chromosome\tlong_exon\texon\t%i\t%i\t.\t%s\t%i\tTranscript \"$transcript\" ($1)\n",@{$exon{$2}}
	}
	elsif (/Non consensus 5\' intron splice site sequence \((\w+)\) after exon (\d+)/){
		printf "CHROMOSOME_$chromosome\tnon_canonical_3_prime_splice_site\texon\t%i\t%i\t.\t%s\t%i\tTranscript \"$transcript\" ($1)\n",@{$exon{$2}}
	}
	elsif (/Non consensus 3\' intron splice site sequence \((\w+)\) before exon (\d+)/){
		printf "CHROMOSOME_$chromosome\tnon_canonical_5_prime_splice_site\texon\t%i\t%i\t.\t%s\t%i\tTranscript \"$transcript\" ($1)\n",@{$exon{$2}}
	}
	elsif (/Short intron\s\((.*)\).*for exon\s(\d+)/){
		printf "CHROMOSOME_$chromosome\tshort_intron\texon\t%i\t%i\t.\t%s\t%i\tTranscript \"$transcript\" ($1)\n",@{$exon{$2}}
	}


}

