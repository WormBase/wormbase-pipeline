#!/usr/local/bin/perl

while( <> ) {
	($dir) = split;
	if( -s "$dir/wise/halfwise.out" ) {
		next;
	} else {
		print "$dir\n";
	}
}
