#!/software/worm/perl_510/bin/perl
# takes a .wig file and spltis it by every thousand sequences
# appends a .#number

use strict;

my $c=1;
my $f=$ARGV[0];
open(my $outf,">","$f.$c");

while (<>){
	if ( /^track/ && ($c++ % 1000)== 0){
		close $outf;
		open($outf,">","$f.$c") || die @!;
	}
	print $outf $_;
}
