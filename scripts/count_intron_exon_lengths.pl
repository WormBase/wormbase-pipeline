#!/usr/local/bin/perl


open (GFF, "</wormsrv2/autoace/CHROMOSOMES/CHROMOSOME_X.gff");
while (<GFF>) {

    if (/\#\#sequence-region (\S+) (\d+) (\d+)/) {
	($sequence,$start,$end) = ($1,$2,$3);
	next;
    }
    

    if (/intron\s+(\d+)\s+(\d+)/) {
	$intron++;
	$intron_length = $intron_length + ($2-$1) + 1;
	print "# intron\t$intron\t$intron_length\n";
    }

    if (/exon\s+(\d+)\s+(\d+)/) {
	$exon++;
	$exon_length = $exon_length + ($2-$1) + 1;
	print "# exon\t$exon\t$exon_length\n";
    }


}
close (GFF);

$remainder = $end - $intron_length  - $exon_length;

print "\n\n";
print "# sequence\t\t$end\n";
print "# intron\t$intron\t$intron_length\n";
print "# exon   \t$exon\t$exon_length\n";
print "# intergenic\t\t$remainder\n\n";

exit 0;

