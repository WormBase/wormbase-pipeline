#!/usr/local/bin/perl5.6.0

use Ace;
use Wormbase;
use strict;

my $maintainer = "All";
my $rundate    = `date +%y%m%d`; chomp $rundate;

my $log = "/wormsrv2/logs/geneace_check.log.$rundate";
open(LOG,">$log")|| die "cant open $log";
print LOG "$0 started at ",`date`,"\n";
print LOG "=============================================\n";


# open a local database connection
my $db = Ace->connect(-path  => '/wormsrv1/geneace/');

my @loci = $db->fetch(-class => 'Locus',
		      -name  => '*');
my $locus;
my $logstring;
my $checklocus = 0;
if ($checklocus == 1)
  {
    foreach $locus (@loci)
      {
	#print "\nLOCUS: $locus\n";
	my $warnings;
	# check for various potential errors in the Locus object
	$warnings = &test_locus_for_errors($locus);
	if(defined($warnings))
	  {
	    print $warnings;
	  }
	
	$logstring .= OtherNameCheck($db,$locus,@loci);
      }
  }
      
my $CDBstring .= FindCurrent_DB_loci(@loci);

print LOG "$CDBstring";

#print the output of the OtherNameCheck subroutine
print LOG "\n\n-------------------------------------\n
These are loci that have Name.Other_name tags, but the locus object being\n
referered to does not have a Name.New_name tag\n\n";
print LOG "$logstring\n";
print LOG "===========================================\n";

#######################################################################################
#This part finds all sequence objs that have a Prov, Detailed or Concise descptn AND a 
#Locus_genomic_seq.

#Erich Schwarz wants this info - emsch@its.caltech.edu
my $findLGS = 0;

if ($findLGS == 1)
{

my $LSGlog;
  my $query = <<END;
  Find Sequence Locus_genomic_seq && (Provisional_description || Detailed_description || Concise_description)
END
      
  my @oddies = $db->fetch(-query=>"$query");
$LSGlog .= "\n\n\nThe folllowing sequences have one of the description fields\n
(Provisional,Detailed or Concise) and also refer to a Locus_genomic_seq\n
  ------------------------------------------------------\n\n";

foreach my $seq(@oddies)
  {
    ($locus) = $seq->at('Visible.Locus_genomic_seq');
    $LSGlog .= "$seq - $locus\n";
  }

$LSGlog .= "==============================================END\n";

#mail Erich
my $interested ="$maintainer"; # "emsch\@its.caltech.edu";
my $mailname = "Annotated Sequences with loci_genomic_seq tags";
#open (OUTLOG,  "|/bin/mailx -s \"$name    \" $maintainer ");
open (OUTLOG,  "|/bin/mailx -s \"$mailname\" $interested");
print OUTLOG "This mail is generated automatically.\nIf you have any queries please email ar2\@sanger.ac.uk\n\n";
print OUTLOG "$LSGlog";
close OUTLOG;

#print to main log file
print LOG "$LSGlog";
}

$db->close;
close LOG;

#$maintainer = "ar2\@sanger.ac.uk";
&mail_maintainer($0,$maintainer,$log);

################################################
sub FindCurrent_DB_loci
{
  my $returnstring;
  my @GEN_loci = @_;
  #make hash of loci (from wormsrv1/geneace)
  my %GEN_loci_check;
  
  foreach my $loci(@GEN_loci)
    {
      $GEN_loci_check{$loci} = 0;
    }
  
  # open a local database connection
  my $db = Ace->connect(-path  => '/wormsrv2/current_DB');

  #check for non-polymorphism loci that are in current_DB but not geneace
  my $query = <<END;
  Find Locus !Polymorphism
END

  my @CDB_loci = $db->fetch(-query=>"$query");
  $db->close;

  $returnstring .= "\n\n\n----------------------------------------------
These loci are in current_DB but not geneace -\n\n";
  foreach my $loci(@CDB_loci)
    {
      unless(defined ($GEN_loci_check{$loci}))
	{
	  $returnstring .= "\t$loci is in current_DB but not geneace!\n"
	}
    }
  $returnstring .= "==========================================\n\n\n\n";
}



###############################################
sub OtherNameCheck
  {
    my($db,$locus,@loci) = @_;
    my ($othername) = $locus->at('Name.Other_name');
    #print "@loci\n";
    if(defined($othername))
      {
	if ("@loci" =~ m/$othername/)
	  {
	    #this locus has $othername as Other_name
	    #but if $othername has this as its Name.New_name thats OK!

	    my ($newloci) = $db->fetch(-class=>'Locus',-name=>"$othername");
	    if (defined($newloci))
	      {
		unless($newloci->at('Name.New_name'))
		  {
		    return  "$locus has Name.Other_name $othername - this exists as separate Locus\n";
		  }
	      }
	  }
      }
    return "";
  }

#################################################
sub test_locus_for_errors{
  my $locus = shift;
  my $warnings;
  my $errors;

  #test for Map AND !NEXT
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
  if(defined($locus->at('Molecular_information.Genomic_sequence')) && !defined($locus->Genomic_sequence)){
    $warnings .= "ERROR 8: $locus has 'Genomic_sequence' tag but no associated value\n";
    $errors++;
  }

  # test for more than one Genomic_sequence tag, but need to allow for splice variants. A bit tricky this
  # and I think my RE (which just looks for word.number.letter) might allow some errors to creep through.  
  if(defined($locus->at('Molecular_information.Genomic_sequence'))){
    my @genomic_sequences = $locus->Genomic_sequence;
    if(scalar(@genomic_sequences)>1){
      my @problems = $locus->at('Molecular_information.Genomic_sequence');
      foreach my $problem (@problems){
	if ($problem !~ m/[\w\d]+\.\d+[a-z]/){
	  $warnings .= "ERROR 9: $locus has multiple 'Genomic_sequence' tags (see $problem)\n";
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

  # test for Enzyme tag (doesn't fit in with new Gene model)
  if(defined($locus->at('Molecular_information.Enzyme'))){  
    $warnings .= "ERROR 11: $locus has an 'Enzyme' tag\n";
    $errors++;
  }

  # test for Canonical_gene tag present
  if(defined($locus->at('Molecular_information.Canonical_gene'))){  
    $warnings .= "ERROR 12: $locus has a 'Canonical_gene' tag\n";
    $errors++;
  }
  
  # test for New_name tag and other information apart from Gene and Species
  if(defined($locus->at('Name.New_name'))){
    if(defined($locus->at('Species')) && defined($locus->at('Type.Gene'))){
    # check for any other info in object
      my @tags = $locus->tags();
      if(scalar(@tags)>3){
	$warnings .= "ERROR 13: New_name tag + extra info. Transfer into new gene?\n";
	$errors++;
      }
    }
    else{
      $warnings .= "ERROR 14: No species and/or Gene tag present\n";
      $errors++;
    }
  }
  
  # Test for Polymorphisms with no P in their title
  if(defined($locus->at('Type.Polymorphism'))){
    if($locus !~ /P/){
      $warnings .= "ERROR 15: No 'P' in title\n";
      $errors++;
    }
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
