#!/usr/local/bin/perl5.8.0 -w
#
# test2test
# 
# by Keith Bradnam, aged 12 and half
#
# A script to convert ?Locus objects into the new ?Gene object
#
# Last updated by: $Author: krb $     
# Last updated on: $Date: 2004-01-30 13:39:01 $   

use strict;
use lib "/wormsrv2/scripts/";
use Wormbase;
use Ace;


# open a local database connection
my $db = Ace->connect(-path  => '/wormsrv1/geneace/');

# two output streams for good (i.e. data converted ok) and bad data
system("cp good_data.ace good_data.backup");
system("cp bad_data.ace bad_data.backup");

open(GOOD,">good_data.ace") || die "Couldn't write to 'good_data.ace'\n";
open(BAD,">bad_data.ace") || die "Couldn't write to 'bad_data.ace'\n";


###################################################################
# Main loop, grab loci objects and loop through each locus object #
###################################################################

$db->timestamps(1);

my @loci = $db->fetch(-class => 'Locus',
		      -name  => '*');

my $errors = 0;
my $id = 1;

# need hash to keep locus name -> WBgn name connections
my %genes;

foreach my $locus (@loci){
  my $cgc_name;
  my $public_name;
  my $object;
  my $species = $locus->Species;

  print "\nLOCUS: $locus\n";

  # check for various potential errors in the Locus object
  my $warnings = &test_locus_for_errors($locus);

  # grab Locus as an object
  $object = $locus->asAce;
#  my $ts = $locus->timestamp if $locus;
#  print "$ts\n";
  ##########################################
  # NO warnings? Ok to process object
  ##########################################
  if(!defined($warnings)){
    
    ###############################################################################
    # Process if Gene, and try to catch anything that falls through if not a Gene
    ###############################################################################
    if(defined($locus->at('Type.Gene'))){  
      # assign new Gene id and public_name
      ($object,$warnings) = &assign_gene_details($locus,$object,$species,$warnings);
    }    
    elsif(!defined($locus->at('Type.Polymorphism')) && !defined($locus->at('Type.Gene_cluster'))){
      $warnings .= "WARNING: $locus - Unknown locus object\n";
      $errors++;
    }    
  }

  ######################################################
  # Warnings? Don't process, print to separate out file
  ########################################################

  # tidy up object
  $object = &tidy_object($object,$locus);

  if($warnings){    
    print BAD "$warnings\n";
    print BAD "$object\n\n" if ($object);
  }
  else{
    print GOOD "$object\n\n" if ($object);
  }   
}

print "\nNumber of errors to fix = $errors\n";


################
# C'est finis! #
################
close(GOOD);
close(BAD);
exit(0);


#####################################################################################
# The Subroutines!!!
#####################################################################################


sub test_locus_for_errors{
  my $locus = shift;
  my $warnings;

  # test for Map AND !NEXT
  if($locus->at('Map')){
    my $map = $locus->Map;
    if (!defined($map)){
        $warnings .= "ERROR 1: $locus has a 'Map' tag but that tag has no map value!\n";
	$errors++;
      }
  }

  # test for more than one gene class
  if(defined($locus->at('Name.Gene_class'))){
    my @gene_classes = $locus->at('Name.Gene_Class');
    if(scalar(@gene_classes) > 1){
      $warnings .= "ERROR 2: $locus has more than one Gene_class\n";
      $errors++;
    }    
  }

  # test for no Type tag
  if(!defined($locus->at('Type'))){  
    $warnings .= "ERROR 3: $locus has no Type tag present\n";
    $errors++;
  }

  # test for Gene AND !Species 
  if(defined($locus->at('Type.Gene')) && !defined($locus->at('Species'))){
    $warnings .= "ERROR 4: $locus has a 'Gene' tag but not a 'Species' tag\n";;
    $errors++;
  }

  # test for !Gene AND Gene_class 
  if(!defined($locus->at('Type.Gene')) && defined($locus->at('Name.Gene_class'))){
    $warnings .= "ERROR 5: $locus has a 'Gene_class' tag but not a 'Gene' tag\n";
    $errors++;
  }

  # test for polymorphisms without a defined type
#  if(defined($locus->at('Type.Polymorphism')) && !defined($locus->Polymorphism)){
#    $warnings .= "ERROR 6: $locus has a 'Polymorphism' tag but no associated info\n";
#    $errors++;
#  }

  # test for no Gene tag AND Genomic_sequence tag
  if(!defined($locus->at('Type.Gene')) && defined($locus->at('Molecular_information.Genomic_sequence'))){
    $warnings .= "ERROR 7: $locus has 'Genomic_sequence' tag but no 'Gene' tag\n";
    $errors++;
  }

  # test for Genomic_sequence tag but no value   
  if(defined($locus->at('Molecular_information.CDS')) && !defined($locus->CDS)){
    $warnings .= "ERROR 8: $locus has 'CDS' tag but no associated value\n";
    $errors++;
  }

  # test for more than one Genomic_sequence tag, but need to allow for splice variants. A bit tricky this
  # and I think my RE (which just looks for word.number.letter) might allow some errors to creep through.  
  if(defined($locus->at('Molecular_information.CDS'))){
    my @genomic_sequences = $locus->CDS;
    if(scalar(@genomic_sequences)>1){
      my @problems = $locus->at('Molecular_information.CDS');
      foreach my $problem (@problems){
	if ($problem !~ m/[\w\d]+\.\d+[a-z]/){
	  $warnings .= "ERROR 9: $locus has multiple 'CDS' tags (see $problem)\n";
	  $errors++;
	}
      }
    }
  }

  # test for Polymorphism AND Gene tags both present
  if(defined($locus->at('Type.Gene')) && defined($locus->at('Type.Polymorphism'))){  
    $warnings .= "ERROR 10: $locus has a 'Polymorphism' tag AND 'Gene' tag\n";
    $errors++;
  }


  # test for Canonical_gene tag present
  if(defined($locus->at('Molecular_information.Canonical_gene'))){  
    $warnings .= "ERROR 12: $locus has a 'Canonical_gene' tag\n";
    $errors++;
  }
  

  # Look for Gene_class tag in non-gene objects 
  if(!defined($locus->at('Type.Gene'))){
    if(defined($locus->at('Name.Gene_class'))){
      $warnings .= "ERROR 16: Gene_class tag in non-gene!\n";
      $errors++;
    }
  }

  # test for Other_name tag but no value   
  if(defined($locus->at('Name.Other_name')) && !defined($locus->Other_name)){
    $warnings .= "ERROR 17: $locus has 'Other_name' tag but no associated value\n";
    $errors++;
  }
  return($warnings);
}

################################################################
# Routine to assign a new gene id and CGC_name, Public_name etc.
################################################################

sub assign_gene_details{
  my $locus = shift;
  my $object = shift;
  my $species = shift;
  my $warnings = shift;
  my $cgc_name;
  my $sequence_name;
  my @sequence_names;
  my $other_name;
  my $public_name;
  my $id_padded = sprintf "%08d" , $id;
  $id++;

  # get approximate age of gene from timestamp
  my $date = &ts($locus,"Type.Gene");
  $date =~ s/_.*//;

  # Add wormbase ID, Live tag, version number, and version change line, 
  # but first check that gene ID hasn't already been assigned by an earlier cross link.

  my $assigned_gene_name;
  if(defined($genes{$locus})){
    print "$locus = $genes{$locus} - already existed\n";
    $assigned_gene_name = $genes{$locus};
  }
  else{
    $genes{$locus} = "WBgn:".$id_padded;
    print "$locus = WBgn:$id_padded\n";
    $assigned_gene_name = "WBgn:$id_padded";  
  }
  $object =~ s/^Locus :.*/Gene : \"$assigned_gene_name\"\nLive\nVersion 1\nVersion_change 1 now \"Bradnam KR\" Imported \"Initial conversion from geneace\"\nRemark \"Gene created on $date\"/;
  chomp($object); # necessary if we are to append to object, have to remove the trailing "\n"
  
  # create CGC_name if CGC_approved tag is present
  if(defined($locus->at('Type.Gene.CGC_approved'))){
    $cgc_name = $locus;
    $object .= "CGC_name \"$cgc_name\"\n";
  }
  else{
    # Have to be careful, if it is a Gene but not CGC_approved, you need to add the Locus name
    # as an 'Other_name' in $object as well as setting $other_name
    $other_name = $locus;
    $object .= "Other_name \"$other_name\"\n";	
  }
  
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
  # Gene clusters, do similar thing as above
#  if ($object =~ m/Gene_cluster\s+Contains\s+\"(\w{3}-[\.\d]+)\"/){
    while(($object =~ m/Gene_cluster\s+Contains\s+\"(\w{3}-[\d\.]+)\"/)){
      my $match = $1;
      if($genes{$match}){
	$object =~ s/Type\s+Gene_cluster\s+Contains\s+\".*\"/Contains \"$genes{$match}\"/;
      }
      else{
	# Assign new wormbase ID to linked locus name
	my $id_padded = sprintf "%08d" , $id;
	$id++;
	$object =~ s/Type\s+Gene_cluster\s+Contains\s+\".*\"/Contains \"WBgn:$id_padded\"/;
	$genes{$match} = "WBgn:".$id_padded;	
      }
    }
 # }


  return($object);

}

######################################################################################################
# quick and dirty subroutine to grab timestamps for the value associated with any tag (if tag present)
######################################################################################################
sub ts{
  my $locus = shift;
  my $tag = shift;
  my $ts = "";
  $ts = $locus->at($tag)->timestamp if $locus->at($tag);
  return($ts);
}
