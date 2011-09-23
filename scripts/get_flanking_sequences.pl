#!/software/bin/perl -w
#
# get_flanking_sequences.pl
#
# Script to grab two unique flanking sequences after specifying a sequence 
# and two coordinates
#
# Last updated by: $Author: klh $     
# Last updated on: $Date: 2011-09-23 10:38:26 $      

use strict;
use lib $ENV{'CVS_DIR'};                  
use Wormbase;
use Feature_mapper;
use Storable;
use Getopt::Long;


my ($store, $wormbase, $species, $min_length, $test, $is_zero, $infile, $id, $outfile);

my ($seq,$x,$y,$db);

&GetOptions ("store:s"   => \$store,
             "test"      => \$test,
             "species:s" => \$species,
             "start:i"   => \$x,
             "end:i"     => \$y,
             "id:s"      => \$id,
             "db:s"      => \$db,
             "zero"      => \$is_zero,
             "seq:s"     => \$seq,
             "infile:s"  => \$infile,
             "outfile:s" => \$outfile,
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

my ($outfh, @feats);

if (defined $outfile) {
  open $outfh, ">$outfile" or die "Could not open $outfile for writing\n";
} else {
  $outfh = \*STDOUT;
}

if (defined $infile) {
  print STDERR "Reading input from file...\n";
  open my $in_fh, $infile or die "Could not open $infile\n";
  while(<$in_fh>) {
    if (/^(\S+)\s+(\d+)\s+(\d+)\s+(\S+)/) {
      my $feat = {
        seq   => $1,
        start => $2, 
        end   => $3, 
        id    => $4,
      };
      push @feats, $feat;
    }
  }
} elsif (defined $seq and defined $x) {
  $y = $x if not defined $y;

  print STDERR "Looking for flanking sequences to $x - $y in $seq\n";

  push @feats, {
    seq   => $seq,
    start => $x,
    end   => $y,
    id    => defined($id) ? $id : "NO_ID",
  };
} else {
  die "You must supply either -infile <file>, or -seq <seq> -start <start> -end <end>\n";
}


$db ||= $wormbase->orgdb;
my $mapper = Feature_mapper->new($db,0,$wormbase);

foreach my $feat (@feats) {
  my($l, $r) = $mapper->get_flanking_sequence_for_feature($feat->{seq}, $feat->{start}, $feat->{end}, $is_zero, $min_length);

  if($l and $r) {
    print $outfh "SNP: ", $feat->{id}, "\n";
    print $outfh "5'FLANK: $l\n";
    print $outfh "3'FLANK: $r\n";
    print $outfh "COMMENT: $feat->{seq}\n";
    print $outfh "||\n\n"
  }
  else { 
    printf "ERROR: cant find flanks for %s/%d-%d\n\n", $feat->{seq}, $feat->{start}, $feat->{end};
  }
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
