#!/bin/env perl
# script to create lists for the CGC
# usage: 
#    cgc_dump.pl -database DATABASE_DIRECTORY -labs 
#    cgc_dump.pl -database DATABASE_DIRECTORY -geneclass 

use Ace;
use Getopt::Long;

my ($db,$geneclass,$labs);

GetOptions(
  '-database=s'  => \$db,
  '-geneclass'   => \$geneclass,
  '-labs'        => \$labs,
)||die(@!);

my $db = Ace->connect(-path => $db)||die(@!);


if($geneclass){
  my @gene_classes = $db->fetch(Gene_class => '*');

  foreach my $gc (sort @gene_classes){
   printf("\"%s\"\t\"%s\"\t\"%s\"\n",$gc,$gc->Designating_laboratory,join(",",$gc->Description));
 }
}

if ($labs){
  my @laboratories = $db->fetch(-query => 'find Laboratory *;CGC');
  foreach my $l (sort @laboratories){
     my $rep = $l->Representative?$l->Representative->Standard_name : '';
     printf("\"%s\"\t\"%s\"\t\"%s\"\t\"%s\"\n",$l,$l->Allele_designation,$rep,$l->Mail);
  }
}
