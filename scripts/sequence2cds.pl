#!/usr/local/bin/perl5.8.0 -w
#
# sequence2cds.pl
# 
# by Keith Bradnam                         
#
# A script to take ?Sequence objects and make ?CDS (and ?Pseudogene) objects
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2003-08-01 15:56:45 $     

use strict;
use Carp;

##################
# variables etc. #
##################

# Needs ace file dump of sequence class as input
#my $dir      = "/wormsrv2/wormbase/camace/";
#my $file    = glob("$dir/camace_Sequence.ace");

# reset input line separator
$/ = "";

# pattern to match timestamps
my $ts = "-O \"(\\d{4}\-\\d{2}\-\\d{2}_\\d{2}:\\d{2}:\\d{2}\.?\\d?_\\S+|original)\"";

open(IN,"<$ARGV[0]") || carp "Can't open input file\n";
while(<IN>){

  # Misc tidy up of bits of models
  s/From\s+$ts\s+Source_Exons/Source_exons/g;
  s/From_Laboratory/From_laboratory/g;


  # convert things which have CDS tag to ?CDS objects
  if (m/Properties\s+$ts\s+Coding\s+$ts\s+CDS/){
    # Convert to new CDS class
    s/^Sequence :/CDS :/;

    # Need to add tags for SMap
    s/Structure\s+$ts\s+From\s+$ts\s+Source/Sequence/;

    # Get rid of this line (now removed in camace)
    s/Properties\s+$ts\s+Status\s+$ts\s+Annotated\s+$ts\s+\d{4}-\d{2}-\d{2}\s+$ts\s//g;
  }

  # Convert Sequences with Pseudogene tag into ?Pseudogene objects
  elsif(m/Properties\s+$ts\s+Pseudogene/){
    # Convert to new Pseudogene class
    s/^Sequence :/Pseudogene :/;

    # No need for Pseudogene tag now, but these are all Sequence objects so can
    # set 'Type' tag to be 'Coding_pseudogene'
    s/Properties\s+$ts\s+Pseudogene\s+$ts/Coding_pseudogene/;

    # Need to add tags for SMap
    s/Structure\s+$ts\s+From\s+$ts\s+Source/Sequence/;

    
    print;
  }

  # Make changes in Parent sequence objects that might link to CDS objects
  else{
    s/Structure\s+$ts\s+Subsequence\s+$ts/CDS_child/g;
    s/Visible\s+$ts\s+Matching_Genomic/Matching_CDS/g;
    
  }

}
close(IN);

exit(0);

