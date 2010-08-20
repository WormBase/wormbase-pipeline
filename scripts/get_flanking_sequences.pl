#!/software/bin/perl -w
#
# get_flanking_sequences.pl
#
# Script to grab two unique flanking sequences after specifying a sequence 
# and two coordinates
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2010-08-20 09:01:34 $      

use strict;
use lib $ENV{'CVS_DIR'};                  
use Wormbase;
use Feature_mapper;
use Storable;
use Getopt::Long;


my ($help, $store, $wormbase, $species);
my ($seq,$x,$y,$db);
GetOptions ("store:s" => \$store,
	    "x:i"     => \$x,
	    "y:i"     => \$y,
	    "db:s"    => \$db,
	    "seq:s"   => \$seq,
	    "species:s"=>\$species,
	    "help"    => \$help,
	    );

# Display help if required
&usage("Help") if ($help);


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

$db ||= $wormbase->orgdb;
my $mapper = Feature_mapper->new($db,0,$wormbase);
my($l, $r) = $mapper->get_flanking_sequence($seq, $x, $y);
if($l and $r) {print "Left  flank: $l\nRight flank: $r\n";}
else { print "ERROR: cant find flanks for $seq:$x-$y\n";}

exit(0);


##############################################################
#
# Subroutines
#
##############################################################



##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - get_flanking_sequence.pl

=head1 USAGE

=over 4

=item get_flanking_sequence.pl  [-options]

=back

This script finds the flanking sequences of a region on the genome.

script_template.pl MANDATORY arguments:

=over 4

=item -x, starting coordinate in the sequence of your object

=back

=over 4

=item -y, ending coordinate in the sequence of your object

=back

=over 4

=item -seq, the sequence that you are specifying the coordinates of - this can be a chomosome, a superlink or a clone.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -species, the species you are working on
 
=back

=over 4

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
