#!/software/bin/perl -w
#
# get_flanking_sequences.pl
#
# Script to grab two unique flanking sequences after specifying a sequence 
# and two coordinates
#
# Last updated by: $Author: ar2 $     
# Last updated on: $Date: 2009-09-25 08:53:21 $      

use strict;
use lib $ENV{'CVS_DIR'};                  
use Wormbase;
use Feature_mapper;
use Storable;
use Getopt::Long;


my ($store, $wormbase, $species);
my ($seq,$x,$y,$db);
GetOptions ("store:s" => \$store,
	    "x:i"     => \$x,
	    "y:i"     => \$y,
	    "db:s"    => \$db,
	    "seq:s"   => \$seq,
	    "species:s"=>\$species
	    );

my $flanking_seq_length    = 30;

# if you only provide one coordinate, use the x again (i.e. 1 bp feature)
if (!defined $y) {$y = $x+1;$x--}

print "Looking for flanking sequences to $x - $y in $seq\n";

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new(-organism => $species
			     );
}

$db |= $wormbase->orgdb;
my $mapper = Feature_mapper->new($db,undef,$wormbase);
my($l, $r) = $mapper->get_flanking_sequence($seq, $x, $y);
if($l and $r) {print "Left  flank: $l\nRight flank: $r\n";}
else { print "ERROR: cant find flanks for $seq:$x-$y\n";}

exit(0);
