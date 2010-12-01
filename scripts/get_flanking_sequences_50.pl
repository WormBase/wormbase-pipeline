#!/software/bin/perl -w
#
# get_flanking_sequences.pl
#
# Script to grab two unique flanking sequences after specifying a sequence 
# and two coordinates
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2010-12-01 11:59:18 $      

use strict;
use lib $ENV{'CVS_DIR'};                  
use Wormbase;
use Feature_mapper;
use Storable;
use Getopt::Long;


my ($help, $store, $wormbase, $species);
my ($seq,$x,$y,$db, $input, $output);
GetOptions ("store:s" => \$store,
            "db:s"    => \$db,
            "species:s"=>\$species,
	    "input:s" => \$input,
	    "output:s" =>\$output,
            "help"    => \$help,
            );

# Display help if required
&usage("Help") if ($help);


my $flanking_seq_length    = 50;

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new(-organism => $species
                             );
}

$db ||= $wormbase->orgdb;
my $mapper = Feature_mapper->new($db,0,$wormbase);

# get the input file
open(IN, "<$input") || die "Can't open $input\n";
open(OUT, ">$output") || die "Can't open $output\n";

while (my $line = <IN>) {

  next if ($line =~ /^\s*$/);
  my ($seq, $x, $y, $id) = split /\s+/, $line;

  my ($l, $r) = $mapper->get_flanking_sequence($seq, $x, $y, $flanking_seq_length);
  if ($l and $r) {
    print OUT "SNP: $id\n5'_FLANK: $l\n3'_FLANK: $r\nCOMMENT:  $seq\n||\n\n"
  } else { 
    print "ERROR: cant find flanks for $id $seq:$x-$y\n";
  }
}

close(OUT);
close(IN);


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

=head2 NAME - get_flanking_sequence_50.pl

=head1 USAGE

=over 4

=item get_flanking_sequence.pl  [-options]

=back

This script finds 50bp flanking sequences of a region on the genome. This length of flank is required for submissions to dbSNP.

script_template.pl MANDATORY arguments:

=over 4

=item -input, file containing the following columns: sequence_name, start, end, ID

=back

=over 4

=item -output, output filename

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
