#!/usr/local/bin/perl5.8.0 -w
#
# locus2gene.pl
# 
# by Keith Bradnam, aged 12 and half
#
# A script to convert ?Locus objects into the new ?Gene object
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2004-02-04 16:43:50 $   

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

  s/^Locus :.*/Gene : \"$name\"\nLive\nVersion 1\nVersion_change 1 now \"Bradnam KR\" Imported \"Initial conversion from geneace\"\nRemark \"Gene created on $date\"/;

  $id++;
}


# tidy up
close(IN);
close(OUT);
exit(0);


  # Calculate 'Sequence_name' if Genomic_sequence tag is present
  if(defined($locus->at('Molecular_information.Genomic_sequence'))){
    @sequence_names = $locus->at('Molecular_information.Genomic_sequence');
    foreach my $sequence (@sequence_names){
      $sequence =~ s/[a-z]$// if ($sequence =~ m/\d[a-z]$/);
      $object .= "Sequence_name \"$sequence\"\n";
      $sequence_name = "$sequence";
    }
  }
  
  # Assign a 'public_name' based on a priority CGC_name -> Sequence_name -> Other_name
  # For non-elegans genes, maybe not use sequence name (where it exists) for public_name
  if($cgc_name){
    $public_name = $cgc_name;
  }
  elsif($sequence_name){
    $public_name = $sequence_name;
    ($public_name = $other_name) if($species ne "Caenorhabditis elegans");
  }
  else{
    $public_name = $other_name;
  }
  $object .= "Public_name \"$public_name\"\n\n";
 
  return($object,$warnings);

}

###############################################################################
# Tidy up ace file to remove leading tags and convert tags for new gene model #
# fetch timestamps for *some* terminal parts of tree object                   #
###############################################################################


sub tidy_object{
  my $object = shift;
  my $locus = shift;
  my $ts;


  # modify Name subtree
  $ts = &ts($locus,"Name.Other_name");
  $object =~ s/Name\s+Other_name\s+(\S+)/Other_name $1 -O \"$ts\" /g;
  $ts = &ts($locus,"Name.Gene_class");
  $object =~ s/Name\s+Gene_class\s+(\"\w+\")/Gene_class $1 -O \"$ts\" /g;
  $ts = &ts($locus,"Name.Old_name");
  $object =~ s/Name\s+Old_name\s+(\S+)/Other_name $1 -O \"$ts\" /g;

  # modify Gene_information subtree
  $object =~ s/Gene_information\s+//g;


  # Modify Type subtree
  $object =~ s/Type\s+Gene\s+Reference_Allele/Reference_allele /g;
  $object =~ s/Type\s+Gene\s+Phenotype/Phenotype /g;
  $object =~ s/Type\s+Gene\s+RNAi_result/RNAi_result /g;
  $object =~ s/Type\s+Gene\s+Complementation_data/Complementation_data /g;
  $object =~ s/Type\s+Gene\s+Controlled_phenotype/Controlled_phenotype /g;
  $object =~ s/Type\s+Gene\s+CGC_approved/CGC_approved /g;
  $object =~ s/Type\s+Gene\s+CGC_unresolved/CGC_unresolved /g;
  $object =~ s/Type\s+Gene\s+Pseudogene/Pseudogene /g;

  # Replace any remaining single Gene tags
  $object =~ s/Type\s+Gene\s+//g;
  $object =~ s/Type\s+Polymorphism\s+SNP/SNP /g;
  $object =~ s/Type\s+Polymorphism\s+Status/Status /g;
  $object =~ s/Type\s+Polymorphism\s+SNP_assay/SNP_assay /g;
  $object =~ s/Type\s+Polymorphism\s+RFLP/RFLP /g;
  $object =~ s/Type\s+Polymorphism\s+Transposon_insertion/Transposon_insertion /g;
  $object =~ s/Type\s+Polymorphism\s+Detection_method/Detection_method /g;
  $object =~ s/Type\s+Polymorphism/Polymorphism /g;
  $object =~ s/Type\s+PCR_product/PCR_product /g;

  # modify Molecular_information subtree, 
  # Genomic_sequence becomes CDS but only for non RNA genes.  Otherwise, it becomes a transcript

  if ($object =~ m/\nRNA_gene/){
    $object =~ s/Molecular_information\s+Genomic_sequence\s+/Transcript /g;
    $object =~ s/RNA_gene//;
  }
  else{
    $object =~ s/Molecular_information\s+Genomic_sequence\s+/CDS /g;
  }
  $object =~ s/Molecular_information\s+Other_sequence\s+/Other_sequence /g;
  $object =~ s/Molecular_information\s+Product\s+/Product /g;
  
  # modify Positive/Negative/Mapping_data subtrees
  $object =~ s/Positive\s+Positive_clone/Positive_clone /g;
  $object =~ s/Positive\s+Inside_rearr/Inside_rearr /g;
  $object =~ s/Negative\s+Negative_clone/Negative_clone /g;
  $object =~ s/Negative\s+Outside_rearr/Outside_rearr /g;
  $object =~ s/Mapping_data\s+Well_ordered/Well_ordered /g;

  #new format for 2_point
  $object =~ s/Mapping_data\s+2_point/Two_pt /g;
  $object =~ s/Mapping_data\s+Multi_point/Multi_point /g;
  $object =~ s/Mapping_data\s+Pos_neg_data/Pos_neg_data /g;

  # modify Contains/Contained_in tags to agree with new Gene_model
  $object =~ s/Contained_in\s+/In_cluster /g;

  # get rid of tags off end of model
  $object =~ s/(Positive_clone\s+\"[\w\#\.]+\"\s+\"[\w \-\,\'\?]+\")\s+(\"[\w \(\),\?\-\;\[\]\.]+\")/$1 -C $2/g;
  $object =~ s/(Negative_clone\s+\"[\w\#\.]+\"\s+\"[\w \-\,\'\?]+\")\s+(\"[\w \(\),\?\-\;\[\]\.]+\")/$1 -C $2/g;
  $object =~ s/(Allele\s+\"[\w\#]+\")\s+(\"[\w \-\,\(\)\@\:\.]+\")/$1 -C $2/g;
  $object =~ s/(Description\s+\"[\w\#\(\)\[\]\.\, \-\:]+\")\s+(\"[\w \-\,]+\")/$1 -C $2/g;
  $object =~ s/(Genomic_sequence\s+\"[\w\#\(\)\[\]\.\, \-]+\")\s+(\"[\w \-\,]+\")/$1 -C $2/g;
  $object =~ s/(Complementation_data\s+\"[\w\#\(\)\[\]\.\, \-]+\")\s+(\"[\w \-\,]+\")/$1 -C $2/g;
  $object =~ s/(Multi_point\s+\"[\w\#\(\)\[\]\.\, \-]+\")\s+(\"[\w \-\,]+\")/$1 -C $2/g;


  # get rid of CGC_approved tag (will become CGC_name in new Gene model)
  $object =~ s/CGC_approved\s+//;

  # get rid of Representative_for associated value (will be generated by XREF from Hide_under)
  $object =~ s/Representative_for\s+\".*\"/Representative_for/g;

  # check for other locus names in object using Hide_under tag
  if ($object =~ m/Hide_under\s+\"(.*)\"/){
    my $match = $1;
    if($genes{$match}){
      $object =~ s/Hide_under\s+\".*\"/Hide_under \"$genes{$match}\"/;
    }
    else{
      # Assign new wormbase ID to linked locus name
      my $id_padded = sprintf "%08d" , $id;
      $id++;
      $object =~ s/Hide_under\s+\".*\"/Hide_under \"WBgn:$id_padded\"/;
      $genes{$match} = "WBgn:".$id_padded;

    }
    
  }


