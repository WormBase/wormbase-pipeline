#!/usr/local/bin/perl


$file = shift;
$offset = shift;


open (FILE, $file);
while (<FILE>) {

    if (/Subsequence\s+(\S+)\s+(\d+)\s+(\d+)/) {
	print "Subsequence\t$1 " .  ($2+$offset) . " " . ($3+$offset) . "\n";
    }
    else {
	print "$_";
    }


}
