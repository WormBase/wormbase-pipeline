#!/usr/local/bin/perl

while (<>){
    if (/^LongText : "(.+)"/) {
	if ($inside) { print "non-termination line $.: $1 inside $inside\n" ; }
	else { $inside = $1 ; }
    }
    elsif (/^LongText/) { 
	chop ; 
	print "wierd line $.: $_" ; 
	if ($inside) { print " inside $inside" ; }
	print "\n" ;
    }
    if (/^\*\*\*LongTextEnd\*\*\*\n$/) { undef $inside ; }
}
