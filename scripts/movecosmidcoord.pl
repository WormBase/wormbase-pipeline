#!/usr/local/bin/perl


$file = shift;
$offset = shift;

open (FILE, $file);
while (<FILE>) {

    if (/DNA\s+(\S+)\s+(\d+)/) {
	print "DNA      $1 " . ($2+$offset). "\n";
    }
    elsif (/Subsequence\s+(\S+)\s+(\d+)\s+(\d+)/) {
	print "Subsequence\t$1 " .  ($2+$offset) . " " . ($3+$offset) . "\n";
    }
    elsif (/Overlap_right\s+(\S+)\s+(\d+)/) {
	print "Overlap_right    $1 ".  ($2+$offset) . "\n";
    }
    elsif (/Clone_left_end\s+(\S+)\s+(\d+)/) {
	print "Clone_left_end    $1 ".  ($2+$offset) . "\n";
    }
    elsif (/Clone_right_end\s+(\S+)\s+(\d+)/) {
	print "Clone_right_end    $1 ".  ($2+$offset) . "\n";
    }
    elsif (/Confirmed_intron\s+(\d+)\s+(\d+)\s+(\S+)/) {
	print "Confirmed_intron         " .  ($1+$offset) . " " . ($2+$offset) . " $3\n";
    }
    elsif (/Assembly_tags\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+.+XS)/) {
	print "Assembly_tags    $1 " . ($2+$offset) . " " . ($3+$offset) . " $4\n";
    }
    elsif (/DNA_homol\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
	print "DNA_homol        $1 $2 $3 ".  ($4+$offset) . " " . ($5+$offset) . " $6 $7\n";
    }
    elsif (/Pep_homol\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
	print "Pep_homol        $1 $2 $3 ".  ($4+$offset) . " " . ($5+$offset) . " $6 $7\n";
    }
    elsif (/Motif_homol\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
	print "Motif_homol        $1 $2 $3 ".  ($4+$offset) . " " . ($5+$offset) . " $6 $7\n";
    }
    elsif (/Feature\s+(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\S+.+)/) {
	print "Feature  $1 ".  ($2+$offset) . " " . ($3+$offset) . " $4 $5\n";
    }
    else {
	print "$_";
    }


}
