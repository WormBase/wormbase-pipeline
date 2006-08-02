#!/usr/bin/perl -w
#===============================================================================
#
#         FILE:  gff2ace.pl
#
#        USAGE:  ./gff2ace.pl 
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
#      CREATED:  18/07/06 15:01:03 BST
#     REVISION:  ---
#===============================================================================

use strict;
use YAML;
use Class::Struct;

struct( Exon => [ start => '$', stop => '$', type => '$', id => '$' ] );
struct( Gene => [ start => '$', stop => '$', exons => '@','orientation' => '$','parent'=>'$' ] );

use Map_Helper;

foreach my $file(@ARGV){

my %genes;
Map_Helper::get_from_gff( $file, 'gene', 'coding_exon', \%genes );

while (my ($k,$v)=each(%genes)){
	my $parent=$v->parent;
	$parent=~s/Contig/Supercontig/;
	print 'Sequence : "',$parent,'"',"\n";
	if ($v->orientation eq '+'){ print 'CDS_child "',$k,'" ',$v->start,' ',$v->stop,"\n\n"}
	else {print 'CDS_child "',$k,'" ',$v->stop,' ',$v->start,"\n\n"}
	print 'CDS : "',$k,'" ',"\n";
	print 'Method "WashU_merged_genes"',"\n";
	my ($start,$stop)=sort{$a <=> $b}($v->start,$v->stop);
	foreach my $exon(sort {$a->start <=> $b->start} @{$v->exons}){
		if ($v->orientation eq '+'){print 'Source_exons ',$exon->start-$start+1,' ',$exon->stop-$start+1,"\n" }
		else{print 'Source_exons ',$stop-$exon->stop+1,' ',$stop-$exon->start+1,"\n" }
	}
	print "\n";
}

}
