#!/usr/local/bin/perl5.6.1 -w

use strict;
use Carp;

##################
# variables etc. #
##################

#my $dir      = "/wormsrv2/wormbase/camace/";
#my $file    = glob("$dir/camace_Sequence.ace");

# reset input line separator
$/ = "";

# pattern to match timestamps
my $ts = "-O \"(\\d{4}\-\\d{2}\-\\d{2}_\\d{2}:\\d{2}:\\d{2}\.?\\d?_\\S+|original)\"";

open(IN,"<$ARGV[0]") || carp "Can't open input file\n";
while(<IN>){
  s/Structure\s+$ts\s+Subsequence\s+$ts/CDS_child/g;
  s/Visible\s+$ts\s+Matching_Genomic/Matching_CDS/g;
  if (m/Properties\s+$ts\s+Coding\s+$ts\s+CDS/){
    # Convert to new CDS class
    s/^Sequence :/CDS :/;

    # Need to add tags for SMap
    s/Structure\s+$ts\s+From\s+$ts\s+Source/Sequence/;

    # Get rid of this line (now removed in camace)
    s/Properties\s+$ts\s+Status\s+$ts\s+Annotated\s+$ts\s+\d{4}-\d{2}-\d{2}\s+$ts\s//g;

    # Get rid of Interpolated_gmap tag (will be succeeded by Chao-Kung's script
    s/Properties\s+$ts\s+Status\s+$ts\s+Annotated\s+$ts\s+\d{4}-\d{2}-\d{2}\s+$ts\s//g;
    s/Interpolated_gMap/Interpolated_map_position/g;
  }
  # Misc tidy up of bits of models
  s/From\s+$ts\s+Source_Exons/Source_exons/g;
  s/From_Laboratory/From_laboratory/g;

  print;

}
close(IN);

exit(0);

