#! /usr/local/bin/perl

#  Usage: conf_intr_genes.pl <-c EST -c homology ...> <GFF file>
#
#  This script takes a sequence file in GFF format 
#  and produces a table with three columns:
#   1. Name of the gene
#   2. Number of introns in the gene
#   3. Number of those introns that are confimed in other ways 
#      (given at command line, e.g. EST, cDNA, homology. If no
#       args are given, the script looks for introns confirmed by 
#       anything)
#   Note: Introns confirmed in UTR are ignored in all totals unless
#         the -utr flag is set.


use strict;
use Getopt::Long;

my ($include_utr_introns, %genes, @confirmed_by);

&GetOptions("c=s@" => \@confirmed_by,
	    'utr' => \$include_utr_introns);

while(<>) {
    /^\#/ and next;
    
    /^\S+\s+\S+\s+intron\s+(\d+)\s+\d+\s+\S+\s+\S+\s+\S+\s+(.+)$/ and do {
	my ($start, $rest) = ($1, $2);
	my ($sub_seq_id, $confirmed) = $rest =~ /Sequence\s+\"(\S+)\"(?:\s+;\s+Confirmed_(\S+))?$/;

	next if not $include_utr_introns and $confirmed eq "in_UTR";

	$genes{ $sub_seq_id }->{'start'} = $start;
	$genes{ $sub_seq_id }->{'introns'}++;

	if ($confirmed) {
	    if (! @confirmed_by or grep { $confirmed eq $_ } @confirmed_by) {
		$genes{ $sub_seq_id }->{'confirmed'}++;
	    }
	}
    };
}

foreach my $gene_id (sort { $genes{$a}->{'start'} <=> $genes{$b}->{'start'} } keys %genes) {
    printf("%-20s%4d %4d\n", 
	   $gene_id, 
	   $genes{ $gene_id }->{'introns'}, 
	   $genes{ $gene_id }->{'confirmed'});
}
