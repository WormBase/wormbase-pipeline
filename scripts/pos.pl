#!/software/bin/perl -w
#
# script_template.pl                           
# 
# by Gary Williams
#
# This is a example of a good script template
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2012-01-12 16:30:47 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Storable;
use Coords_converter;

######################################
# variables and command-line options # 
######################################

my ($test, $store, $wormbase);

$test = 1;

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $ENV{"USER"},
                             -test    => $test,
			     );
}

# parse out the parameters
my ($position1, $position2, $chromosome, $clone);

foreach my $arg (@ARGV) {
  #print "arg = $arg\n";
  if ($arg =~ /^CHROMOSOME_([IVX]+|MtDNA)$/) {
    $chromosome = $arg;
  } elsif ($arg =~ /^\d+$/) {
    if (!defined $position1) {
      $position1 = $arg;
    } else {
      $position2 = $arg;
    }
  } elsif ($arg =~ /^([IVX]+|MtDNA)$/) {
    $chromosome = $arg;
  } elsif ($arg =~ /(\S+)\:(\d+)(\.\.|\-)(\d+)/) {
    $chromosome = $1;
    $position1 = $2;
    $position2 = $4;
  } elsif ($arg =~ /(\d+)(\.\.|\-)(\d+)/) {
    $position1 = $1;
    $position2 = $3;
  } elsif ($arg =~ /\w+/) {
    $clone = $arg;
  }
}



#################################
# Set up some useful paths      #
#################################


my   $db = $wormbase->database('current');


##########################
# MAIN BODY OF SCRIPT
##########################

my $coords = Coords_converter->invoke($db, undef, $wormbase);
my @clone_coords;

# sanity check

if (defined $clone) {
  $clone = uc $clone;
  if ($coords->isa_chromosome($clone)) {
    $chromosome = $clone;
    $clone = undef;
  }
} elsif (defined $chromosome) {
  #$chromosome = uc $chromosome;
  if ($coords->isa_clone($chromosome)) {
    $clone = uc $chromosome;
    $chromosome = undef;
  }
}

print "chromosome=$chromosome\t" if (defined $chromosome);
print "clone=$clone\t" if (defined $clone);
print "pos1=$position1\t" if (defined $position1);
print "pos2=$position2\t" if (defined $position2);
print "\n\n";


if (defined $chromosome && defined $position1) {
  
  if (! defined $position2) {
    $position2 = $position1;
  }


  @clone_coords = $coords->LocateSpan($chromosome, $position1, $position2);


  print "Clone: $clone_coords[0] $clone_coords[1]..$clone_coords[2]\n";

} elsif (defined $clone) {

  my ($chrom, $pos1, $pos2);

  if (! defined $position1) {
    my ($chrom, $offset)       = $coords->CloneOffset($clone);
    print "Chromosome: $chrom $offset\n";

  } else {
    ($chrom, $pos1) = $coords->Coords_2chrom_coords( $clone, $position1);

    
    if (defined $position2) {
      ($chrom, $pos2) = $coords->Coords_2chrom_coords( $clone, $position2);
      print "Chromosome: $chrom $pos1..$pos2\n";
    } else {
      print "Chromosome: $chrom $pos1\n";
    }
  }
}
exit(0);


