#!/usr/local/bin/perl5.6.0 -w
# 
# geneace_check.pl
#
# by Anthony Rogers based on a script by Keith Bradnam
#
# Script to run consistency checks on the geneace database
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2002-07-10 10:57:29 $

use Ace;
use lib "/wormsrv2/scripts/"; 
use Wormbase;
use strict;

our $log;
our $erichlog;
our $errors;
&create_log_files;

# open a connection to geneace and grab list of loci
our $tace = glob("~wormpub/ACEDB/bin.ALPHA_4/tace");   # tace executable path
my $db = Ace->connect(-path  => '/wormsrv1/geneace/',
		      -program =>$tace) || do { print LOG "Connection failure: ",Ace->error; die();};

my @loci = $db->fetch(-class => 'Locus',
		      -name  => '*');

my $size =scalar(@loci);

# Loop through loci checkingfor various potential errors in the Locus object
print "\nChecking loci for errors:\n";
foreach my $locus (@loci){
  print "$locus\n";
  my $warnings;
  my $erich_warnings;
  ($warnings, $erich_warnings) = &test_locus_for_errors($locus);
  print LOG "$warnings" if(defined($warnings));
  #Erich Schwarz wants some of these - emsch@its.caltech.edu
  print ERICHLOG "$erich_warnings" if(defined($erich_warnings));
}

# Look for loci in current_DB not in geneace
print "\nLooking for new loci in /wormsrv2/current_DB:\n\n";
&find_new_loci_in_current_DB($db);

print LOG "\nThere were $errors errors in $size loci.\n";


##########################################
# Tidy up and mail relevant log files
##########################################

$db->close;
close(LOG);
close(ERICHLOG);

my $maintainer = "All";
&mail_maintainer($0,$maintainer,$log);

my $interested ="krb\@sanger.ac.uk, emsch\@its.caltech.edu, ck1\@sanger.ac.uk";
&mail_maintainer($0,"$interested",$erichlog);
exit(0);

################################################
sub find_new_loci_in_current_DB{
  my $db = shift;
  my $warnings;

  # open a database connection to current_DB and grab all loci names (excluding polymorphisms)
  my $new_db = Ace->connect(-path  => '/wormsrv2/current_DB',
		    -program =>$tace) || do { print LOG "Connection failure: ",Ace->error; die();};
  my @current_DB_loci = $db->fetch(-query=>'Find Locus;!Polymorphism');
  $new_db->close;

  #cross reference in geneace
  foreach my $loci(@current_DB_loci){
    my $new_loci = $db->fetch(-class=>'Locus',-name=>"$loci");
    unless(defined($new_loci)){
      $warnings .= "ERROR 19: $new_loci in current_DB is not in /wormsrv1/geneace\n";
      $errors++;
    }
  }
  print LOG "\n$warnings\n" if $warnings;;
}



###############################################

sub test_locus_for_errors{
  my $locus = shift;
  my $warnings;
  my $erich_warnings;

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
  
  # test for no Gene tag AND Genomic_sequence tag
  if(!defined($locus->at('Type.Gene')) && defined($locus->at('Molecular_information.Genomic_sequence'))){
    $warnings .= "ERROR 6: $locus has 'Genomic_sequence' tag but no 'Gene' tag\n";
    $errors++;
  }

  # test for Genomic_sequence tag but no value   
  if(defined($locus->at('Molecular_information.Genomic_sequence')) && !defined($locus->Genomic_sequence)){
    $warnings .= "ERROR 7: $locus has 'Genomic_sequence' tag but no associated value\n";
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
	  $warnings .= "ERROR 8: $locus has multiple 'Genomic_sequence' tags (see $problem)\n";
	  $errors++;
	}
      }
    }

    # can also test to see if there are attached annotations to these sequences which Erich should
    # know about, i.e. they should now be attached to the loci object instead
    foreach my $seq (@genomic_sequences){
      my ($newseq) = $db->fetch(-class=>'Sequence',-name=>"$seq");
      if(defined($newseq->at('Visible.Provisional_description')) || 
	 defined($newseq->at('Visible.Concise_description')) ||
	 defined($newseq->at('Visible.Detailed_description'))){  
	$erich_warnings .= "$seq has attached functional annotation which should now be attached to $locus\n";
	$errors++;
      }
    }
  }
  
  # test for Polymorphism AND Gene tags both present
  if(defined($locus->at('Type.Gene')) && defined($locus->at('Type.Polymorphism'))){  
    $warnings .= "ERROR 9: $locus has a 'Polymorphism' tag AND 'Gene' tag\n";
    $errors++;
  }

  # test for Enzyme tag (doesn't fit in with new Gene model)
  if(defined($locus->at('Molecular_information.Enzyme'))){  
    $warnings .= "ERROR 10: $locus has an 'Enzyme' tag\n";
    $errors++;
  }

  # test for Canonical_gene tag present
  if(defined($locus->at('Molecular_information.Canonical_gene'))){  
    $warnings .= "ERROR 11: $locus has a 'Canonical_gene' tag\n";
    $errors++;
  }
  
  # test for New_name tag and other information apart from Gene and Species
  if(defined($locus->at('Name.New_name'))){
    if(defined($locus->at('Species')) && defined($locus->at('Type.Gene'))){
      # check for any other info in object
      my @tags = $locus->tags();
      if(scalar(@tags)>3){
	$warnings .= "ERROR 12: $locus has New_name tag + extra info. Transfer into new gene?\n";
	$errors++;
      }
    }
    else{
      $warnings .= "ERROR 13: $locus has no species and/or Gene tag present\n";
      $errors++;
    }
  }
  
  # Test for Polymorphisms with no P in their title
  if(defined($locus->at('Type.Polymorphism'))){
    if($locus !~ /P/){
      $warnings .= "ERROR 14: $locus has no 'P' in title\n";
      $errors++;
    }
  }

  # Look for Gene_class tag in non-gene objects 
  if(!defined($locus->at('Type.Gene'))){
    if(defined($locus->at('Name.Gene_class'))){
      $warnings .= "ERROR 15: $locus has Gene_class tag but it is not a gene!\n";
      $errors++;
    }
  }
  
  # test for Other_name tag but no value   
  if(defined($locus->at('Name.Other_name')) && !defined($locus->Other_name)){
    $warnings .= "ERROR 16: $locus has 'Other_name' tag but no associated value\n";
    $errors++;
  }

  # test for Other_name value which is also a Locus name in its own right
  # N.B. This is not always a bad thing
  if(defined($locus->at('Name.Other_name'))){
    my @other_names = $locus->Other_name;
    foreach my $other_name (@other_names){
      my ($newloci) = $db->fetch(-class=>'Locus',-name=>"$other_name");
      if(defined($newloci)){
	unless($newloci->at('Name.New_name')){
	  $warnings .= "ERROR 17: $locus has Other_name: $other_name which is also a separate locus object\n";
	  $errors++;
	}   
      }
    }
  }
  # Remind of outstanding CGC_unresolved tags
  if(defined($locus->CGC_unresolved)){
    my ($unresolved_details) = $locus->at('Type.Gene.CGC_unresolved');
    $warnings .= "ERROR 18: $locus has CGC_unresolved tag: \"$unresolved_details\"\n";
    $errors++;
  }

  return($warnings, $erich_warnings);

}

#####################################################################

sub create_log_files{
  my $rundate    = `date +%y%m%d`; chomp $rundate;
  $log = "/wormsrv2/logs/geneace_check.log.$rundate.$$";
  $erichlog = "/wormsrv2/logs/geneace_check.erichlog.$rundate.$$";
  open(LOG,">$log") || die "cant open $log";
  print LOG "$0 started at ",`date`,"\n";
  print LOG "=============================================\n";

  open(ERICHLOG,">$erichlog") || die "cant open $erichlog";
  print ERICHLOG "$0 started at ",`date`,"\n";
  print ERICHLOG "=============================================\n";
  print ERICHLOG "This mail is generated automatically.\n";
  print ERICHLOG "If you have any queries please email ar2\@sanger.ac.uk or krb\@sanger.ac.uk\n\n";
}

