#!/usr/local/bin/perl5.8.0 -w
#
# locus2gene.pl
# 
# by Keith Bradnam, aged 12 and half
#
# A script to convert ?Locus objects into the new ?Gene object
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2004-02-06 17:44:57 $   

use strict;
use lib -e "/wormsrv2/scripts" ? "/wormsrv2/scripts" : $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;


################################################################
# set up input/output streams
################################################################

my $file; # where input file is located
GetOptions ("file=s"   => \$file);

open(IN,"<$file") || die "Can't open input file\n";
open(OUT,">$file.gene") || die "Couldn't open output file\n";



################################################################
# misc. variables
################################################################

# pattern to match timestamps
my $ts = "-O \"(\\d{4}\-\\d{2}\-\\d{2}_\\d{2}:\\d{2}:\\d{2}\.?\\d?_\\S+|original)\"";
my $id = 1; # the numerical part of the gene ID
my $date = "2004-02-04";
$/ = ""; # reset input line separator



# process input file, changing lines as appropriate

while(<IN>){

  my $id_padded = sprintf "%08d" , $id;
  my $name = "WBGene$id_padded";  
  my ($public_name, $cgc_name, $non_cgc_name, $sequence_name, $other_name);



  # Change Non_CGC_name field to Other_name

#  s/Non_CGC_name\s+$ts\s+(\"[\w\d]+\.[\d\w:]+\")\s+$ts\s+(\d+)\s+$ts\s+(\d+)\s+$ts/CDS $3 $5 $7/g;
  if(m/Non_CGC_name/){
    s/Non_CGC_name\s+$ts\s+(\S+)\s+$ts/Other_name $1 $2 $3/;
    $non_cgc_name = $2;
  }
  # Change Transcript_name or Pseudogene_name to Sequence_name
  s/Transcript_name/Sequence_name/g if (m/Transcript_name/);
  s/Pseudogene_name/Sequence_name/g if (m/Pseudogene_name/);



  # grab other names when present: e.g. cgc_name, sequence_name or other_name
  if(m/\s+CGC_name\s+$ts\s+(\S+)\s+$ts/){
    $cgc_name = $2;
  }
  if(m/\s+Sequence_name\s+$ts\s+(\S+)\s+$ts/){
    $sequence_name = $2;
  }
  if(m/\s+Other_name\s+$ts\s+(\S+)\s+$ts/){
    $other_name = $2;
  }


  # Assign a 'public_name' based on a priority CGC_name -> Sequence_name -> Other_name
  # For non-elegans genes, maybe not use sequence name (where it exists) for public_name
  if($cgc_name){
    $public_name = $cgc_name;
  }
  elsif($sequence_name){
    $public_name = $sequence_name;
  }
  elsif($non_cgc_name){
    $public_name = $non_cgc_name;
  }
  else{
    $public_name = $other_name;
  }


  # Main conversion to Gene object with basic history/identity information
  s/^Locus :.*/Gene : \"$name\"\nLive\nVersion 1\nVersion_change 1 now \"WBPerson1971\" Imported \"Initial conversion from geneace\"\nPublic_name $public_name/;

  # Remove CGC_approved, now assumed if CGC_name is present
  s/Type\s+$ts\s+Gene\s+$ts\s+CGC_approved\s+$ts//;

  # get rid of Type.Gene part of tree as this is no longer applicable
  s/Type\s+$ts\s+Gene\s+$ts\s+//g;  

  # change tag names to new shorter formats
  s/Gene_information\s+/Gene_info /g;
  s/Molecular_information\s+/Molecular_info /g;


  print;

  $id++;
  last if $id>205;
}


# tidy up
close(IN);
close(OUT);
exit(0);


 

###############################################################################

