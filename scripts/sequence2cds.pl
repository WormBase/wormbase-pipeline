#!/usr/local/bin/perl5.8.0 -w
#
# sequence2cds.pl
# 
# by Keith Bradnam                         
#
# A script to take ?Sequence objects and make ?CDS (and ?Pseudogene) objects
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2003-08-05 14:26:50 $     

use strict;
use Carp;
use Getopt::Long;


##################
# variables etc. #
##################

my ($camace,$stlace, $db, $dbdir);

# pattern to match timestamps
my $ts = "-O \"(\\d{4}\-\\d{2}\-\\d{2}_\\d{2}:\\d{2}:\\d{2}\.?\\d?_\\S+|original)\"";

# store pseudogene names in hash
my %pseudogenes;

GetOptions ("camace"      => \$camace,
            "stlace"      => \$stlace);


if($camace){
  $db = "camace";
  $dbdir = "/wormsrv2/wormbase/camace/";
}

if($stlace){
  $db = "stlace";
  $dbdir = "/wormsrv2/wormbase/stlace/";
}


# reset input line separator
$/ = "";



# Two output files, one for pseudogenes, one for everything else.
open(OUT,">${dbdir}krb_${db}_Sequence.ace") || croak "Couldn't open output file\n";
open(PSEUDOGENES,">${dbdir}krb_${db}_Pseudogenes.ace") || croak "Couldn't open pseudogenes file\n";


open(IN,"<${dbdir}${db}_Sequence.ace") || carp "Can't open input file\n";

while(<IN>){

  # Misc tidy up of bits of models
  s/From\s+$ts\s+Source_Exons/Source_exons/g;
  s/From_Laboratory/From_laboratory/g;


  # convert things which have CDS tag to ?CDS objects
#  if (m/Properties\s+$ts\s+Coding\s+$ts\s+CDS/){
    # Convert to new CDS class
#    s/^Sequence :/CDS :/;

    # Need to add tags for SMap
#    s/Structure\s+$ts\s+From\s+$ts\s+Source/Sequence/;

    # Get rid of this line (now removed in camace)
#    s/Properties\s+$ts\s+Status\s+$ts\s+Annotated\s+$ts\s+\d{4}-\d{2}-\d{2}\s+$ts\s//g;
#  }

  # Convert Sequences with Pseudogene tag into ?Pseudogene objects
  if(m/Properties\s+$ts\s+Pseudogene/){
 
    # Convert to new Pseudogene class
    s/^Sequence :/Pseudogene :/;

    # Track pseudogene names for changing parent ?Sequence objects later
    m/^Pseudogene : (\".*?\")/;
    $pseudogenes{$1} = 1;

    # No need for Pseudogene tag now, but these are all Sequence objects so can
    # set 'Type' tag to be 'Coding_pseudogene'
    s/Properties\s+$ts\s+Pseudogene\s+$ts/Coding_pseudogene/;

    # Need to add tags for SMap
    s/Structure\s+$ts\s+From\s+$ts\s+Source/Sequence/;

    
    print PSEUDOGENES;
  }


  # Make changes in Parent sequence objects that might link to CDS objects
  else{
#    s/Structure\s+$ts\s+Subsequence\s+$ts/CDS_child/g;
#    s/Visible\s+$ts\s+Matching_Genomic/Matching_CDS/g;
    print OUT;
  }

}

close(IN);
close(OUT);
close(PSEUDOGENES);


# Now need to re-read the Sequence output file to replace 'Subsequence' with
# 'Pseudogene' for those things which will be pseudogenes

# Need new output file
open(IN,"<${dbdir}krb_${db}_Sequence.ace") || croak "Can't open input file\n";

open(OUT,">${dbdir}krb2_${db}_Sequence.ace") || croak "Couldn't open output file\n";

# reset input line separator
$/ = "\n";


while(<IN>){

  if(m/^Structure\s+.*\s+Subsequence\s+$ts\s+(\".*\")\s+$ts\s+\d+/){
    my $subsequence = $2;
    s/Subsequence/Pseudogene/ if ($pseudogenes{$subsequence});           
  }

  print OUT;
}


close(IN);
close(OUT);

system("mv ${dbdir}krb2_${db}_Sequence.ace ${dbdir}krb_${db}_Sequence.ace") && croak "Couldn't overwrite file\n";

exit(0);

