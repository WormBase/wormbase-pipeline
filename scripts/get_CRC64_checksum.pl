#!/usr/local/bin/perl5.6.0 

# Chao-kung Chen 
# July 29, 2002 @ recorded temperature 30 degree Celcius
# What it does: get CRC64 checksum appeared in SQ line of SwissProt flatflie
# This script revises Dan's get_checksum.pl, which generates GCG checksum


use lib "/wormsrv2/scripts";
use CRC64;

my $file = shift;
my $seq = "";
my $name = "";


open (SEQ, "<$file") or die "Can't open sequence file $!\n";
while (<SEQ>) {
    chomp;
    if (/^>(\S+)/) {
	$name = $1;
	next;
    }
    $seq .= $_;
}
close SEQ;

# lower to upper-case
$seq =~ tr/a-z/A-Z/;

$crc = SWISS::CRC64::crc64("$seq");  #returns a string like "E3DCADD69B01ADD1"

print "CRC64_checksum: $crc\n";



