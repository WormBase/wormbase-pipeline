#!/usr/local/bin/perl



# get data from camace

use Ace ;
$campath = "/wormsrv2/camace" ;
$camdb = Ace->connect(-path=>"$campath") || 
    die "failed to connect to database\n" ;

# then some info about current links

$it = $camdb->fetch_many(Sequence => 'SUPERLINK*') ;
while ($obj = $it->next) {

    print "\nSUPERLINK object : \'$obj\'\n";
    foreach $a ($obj->at('Structure.Subsequence')) {
	($seq, $start, $end) = $a->row ;
	print " $seq \t$start -> $end\n";

	if ($seq =~ /\./) {
	    print "CDS object in SUPERLINK. Abort \'makelinks\' \n";
	   # exit;
	}

    }
}


exit;


