#!/usr/local/bin/perl5.8.0 -w
#
# sequence2cds.pl
# 
# by Keith Bradnam                         
#
# A script to take ?Sequence objects and make ?CDS objects
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2003-11-06 18:50:22 $     

use strict;
use Carp;
use Getopt::Long;


##################
# variables etc. #
##################

my $input; # where input file is located
my $outdir;   # path to output directory

# pattern to match timestamps
my $ts = "-O \"(\\d{4}\-\\d{2}\-\\d{2}_\\d{2}:\\d{2}:\\d{2}\.?\\d?_\\S+|original)\"";


GetOptions ("input=s"   => \$input,
	    "outdir=s"  => \$outdir);


if(!$outdir){
  $outdir = "/nfs/disk100/wormpub/DATABASES/TEST_DBs/cds_ace/acefiles";
}

# Open output/input streams
open(IN,"<$input") || carp "Can't open input file\n";
open(SEQ,">$outdir/camace_new_Sequence.ace") || croak "Couldn't open output file\n";
open(CDS,">$outdir/camace_CDS.ace") || croak "Couldn't open output file\n";


# reset input line separator
$/ = "";


while(<IN>){

  # Misc tidy up of bits of models
  s/From\s+$ts\s+Source_Exons/Source_exons/g;
  s/From_Laboratory/From_laboratory/g;
  s/From_Database/From_database/g;
  s/From_Author/From_author/g;


  # convert things which have CDS tag to ?CDS objects
  if (m/Properties\s+$ts\s+Coding\s+$ts\s+CDS/){
    # Convert to new CDS class
    s/^Sequence :/CDS :/;

    # Need to add tags for SMap
    s/Structure\s+$ts\s+From\s+$ts\s+Source/Sequence/;

    # Get rid of this line (now removed in camace)
    s/Properties\s+$ts\s+Status\s+$ts\s+Annotated\s+$ts\s+\d{4}-\d{2}-\d{2}\s+$ts\s//g;

    # Change Has_allele tag
    s/Has_allele/Allele/;
    
    print CDS;
  }


  # Make changes in Parent sequence objects that might link to CDS objects
  else{
    s/Structure\s+$ts\s+Subsequence\s+$ts/CDS/g;
    s/Visible\s+$ts\s+Matching_Genomic/Matching_CDS/g;
    print SEQ;
  }

}

close(IN);
close(SEQ);
close(CDS);


# reset input line separator
$/ = "\n";


exit(0);

