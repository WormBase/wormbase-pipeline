#!/usr/local/bin/perl5.6.0 
#
# get_checksum
#
# dl



$|=1;
BEGIN {
  unshift (@INC,"/nfs/disk92/PerlSource/Bioperl/Releases/bioperl-0.05");
}

use Bio::Seq;
use Getopt::Std;

getopts ('f');

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

# bioperl stuff
my $bioseq = Bio::Seq->new(-seq=>$seq,-ffmt=>'Fasta',-type=>'Dna',);
my $chksum = $bioseq->GCG_checksum;
   

if ($opt_f) {
    print ">$name $chksum\n";
    my ($newline,$linestart);
    my $size     = length ($seq);
    my $no_lines = int ($size / 60) +1; 
    for ($i = 0; $i < $no_lines; $i++) {
	$linestart = $i * 60;
	$newline = substr($seq,$linestart,60);
	print "$newline\n";
    }
}
else {
    print "Chksum for sequence $name is $chksum\n";
}

exit(0); 
