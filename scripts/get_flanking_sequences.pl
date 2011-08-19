#!/software/bin/perl -w
#
# get_flanking_sequences.pl
#
# Script to grab two unique flanking sequences after specifying a sequence 
# and two coordinates
#
# Last updated by: $Author: mt3 $     
# Last updated on: $Date: 2011-08-19 13:49:14 $      

use strict;
use lib $ENV{'CVS_DIR'};                  
use Wormbase;
use Feature_mapper;
use Storable;
use Getopt::Long;


my ($store, $wormbase, $species, $min_length, $test, $is_zero);

my ($seq,$x,$y,$db);

&GetOptions ("store:s"   => \$store,
             "test"      => \$test,
             "species:s" => \$species,
             "start:i"   => \$x,
             "end:i"     => \$y,
             "db:s"      => \$db,
             "zero"      => \$is_zero,
             "seq:s"     => \$seq,
             "length:s"  => \$min_length, 
	    );

$min_length = 30 if not defined $min_length;

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new(
    -organism => $species,
    -test => $test);
}

$y = $x if not defined $y;

print STDERR "Looking for flanking sequences to $x - $y in $seq\n";

$db ||= $wormbase->orgdb;
my $mapper = Feature_mapper->new($db,0,$wormbase);
my($l, $r) = $mapper->get_flanking_sequence_for_feature($seq, $x, $y, $is_zero, $min_length);
if($l and $r) {
  print "Left  flank: $l\nRight flank: $r\n";
}
else { 
  print "ERROR: cant find flanks for $seq:$x-$y\n";
}

exit(0);


__END__

=pod

=head2 NAME - get_flanking_sequence.pl

=head1 USAGE

=over 4

=item get_flanking_sequence.pl  [-options]

=back

This script finds the flanking sequences of a region on the genome.

=over 4

=item -seq - the sequence that you are specifying the coordinates of - this can be a chomosome, a superlink or a clone.

=item -start - starting coordinate in the sequence of your object

=item -end - ending coordinate in the sequence of your object

=item -zero - raise this if your extent is a 2bp one defining a 0-bp object (i.e. between bases)

=item -length - the length of the flank requried - the default is 30bp

=item -species - the species you are working on - the default is elegans

=item -db, the database you wisk to look in - the default is autoace

=back

=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item mh6 (mh6@sanger.ac.uk)

=back

=cut
