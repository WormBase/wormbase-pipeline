#!/usr/local/bin/perl5.8.0 -w
# 
# geneace_check.pl
#
# Initiated by Keith Bradnam
#
# Script to run consistency checks on the geneace database
#
# Last updated by: $Author: ck1 $
# Last updated on: $Date: 2003-09-08 14:18:55 $


use strict;
use lib "/wormsrv2/scripts/"; 
use Wormbase;
use Ace;
use Carp;
use Getopt::Long;



###################################################
# Miscellaneous important variables               # 
###################################################

my $tace = &tace;                                          # tace executable path
my $curr_db = "/nfs/disk100/wormpub/DATABASES/current_DB"; # Used for some cross checking with geneace
my $def_dir = "/wormsrv1/geneace/wquery";                  # where lots of table-maker definitions are kept
my $rundate = `date +%y%m%d`; chomp $rundate;              # Used by various parts of script for filename creation
my $maintainers = "All";                                   # Default for emailing everyone logfile
my $caltech_errors = 0;                                    # counter for tracking errors going into Caltech email
my $jah_errors = 0;                                        # counter for tracking errors going into Caltech email
my $log;                                                   # main log file for most output
my $caltech_log;                                           # Additional log file for problems that need to be sent to Caltech
my $jah_log;                                               # Additional log file for problems to be sent to Jonathan Hodgkin at CGC
my (%L_name_F_WBP, %L_name_F_M_WBP);                       # hashes for checking Person and Author merging?

# list of hard-coded loci with other-name(s) as valid independent loci (exceptions for main name / other_name merging) 
# @exceptions and %exceptions are made global as they are used for checking both Locus and Strain classes
my @exceptions = 
  qw (aka-1 bar-1 cas-1 clh-2 clh-3 ctl-1 ctl-2 egl-13 evl-20 gst-4 mig-1 sle-1 slo-1 rap-1 rpy-1 dmo-1 mod-1 
      old-1 plk-1 ptp-3 rab-18 rsp-1 rsp-2 rsp-4 rsp-5 rsp-6 sca-1 sus-1 twk-1 twk-2 twk-3 twk-4 twk-5 twk-6 twk-7 twk-8
      twk-9 twk-10 twk-11 twk-12 twk-13 twk-14 twk-16 twk-17 twk-18 twk-20 twk-21 twk-22 twk-23 twk-24 twk-25 twk-26 
      twk-32 unc-58 sup-9
     );

# load elements of @exceptions into hash for fast checking
my %exceptions;
foreach (@exceptions){$exceptions{$_}++}; 


###################################################
# command line options                            # 
###################################################

my ($help, $debug, $class, @classes, $database, $ace, $verbose);

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "class=s"     => \@classes,
	    "database=s"  => \$database,
            "ace"         => \$ace, 
	    "verbose"     => \$verbose
           );


# Display help if required
if ($help){&usage("Help")}


# Use debug mode?
if($debug){
  print "\nDEBUG = \"$debug\"\n";
  ($maintainers = "$debug" . '\@sanger.ac.uk');
}


# choose database to query: default is /wormsrv1/geneace 
my $default_db = "/wormsrv1/geneace";
($default_db = $database) if ($database);
print "\nUsing $default_db as default database.\n\n";



# Open file for acefile output?
if ($ace){
  my $acefile = "/wormsrv1/geneace/CHECKS/geneace_check.$rundate.$$.ace";
  open(ACE, ">>$acefile") || croak $!;
  system("chmod 777 $acefile");
}


# Set up other log files                         
&create_log_files;




######################################################################################################
######################################################################################################
##                                                                                                  ##
##                                       MAIN PART OF SCRIPT                                        ##
##                                                                                                  ##
######################################################################################################
######################################################################################################

# open a connection to database
my $db = Ace->connect(-path  => $default_db,
		      -program =>$tace) || do { print LOG "Connection failure: ",Ace->error; die();};


# Process separate classes if specified on the command line else process all classes
@classes = ("locus","laboratory","allele","strain","rearrangement","sequence","mapping","evidence") if (!@classes);

foreach $class (@classes){
  if ($class =~ m/locus/i)           {&process_locus_class}
  if ($class =~ m/(laboratory)/i)    {&process_laboratory_class}
  if ($class =~ m/allele/i)          {&process_allele_class}
  if ($class =~ m/strain/i)          {&process_strain_class}
  if ($class =~ m/(rearrangement)/i) {&process_rearrangement}
  if ($class =~ m/(sequence)/i)      {&process_sequence}
  if ($class =~ m/(mapping)/i)       {&check_genetics_coords_mapping}
  if ($class =~ m/(evidence)/i)      {&check_evidence}
}  


#######################################
# Tidy up and mail relevant log files #
#######################################

$db->close;
close(LOG);
close(CALTECHLOG);
close(JAHLOG);
close(ACE) if ($ace);


# Always mail to $maintainers (which might be a single user under debug mode)
mail_maintainer($0,$maintainers,$log);


# Also mail to Erich unless in debug mode or unless there is no errors
my $interested ="krb\@sanger.ac.uk, emsch\@its.caltech.edu, kimberly\@minerva.caltech.edu, ck1\@sanger.ac.uk";
mail_maintainer($0,"$interested",$caltech_log) unless ($debug || $caltech_errors == 0); 


# Email Jonathan Hodgkin subset of errors that he might be able to help with unless 
# in debug mode or no errors
mail_maintainer($0,"cgc\@wormbase.org",$jah_log) unless ($debug || $jah_errors == 0);


exit(0);



##############################################################################################
#
#
#                                    geneace_check subroutines
#
#
#
##############################################################################################


sub process_locus_class{

  # get all loci
  my @loci = $db->fetch(-class => 'Locus',
		        -name  => '*');
  

  # Loop through loci checking for various potential errors in the Locus object
  print "Checking Locus class for errors:\n" if ($verbose);
  print LOG "\nChecking Locus class for errors:\n--------------------------------\n";

  foreach my $locus (@loci){
    print "$locus " if ($verbose);  # useful to see where you are in the script if running on command line
    my $warnings;
    ($warnings) = &test_locus_for_errors($locus);
    print "\n" if ($verbose);
    print LOG "$warnings" if(defined($warnings));
    undef($locus);
  }


  # Look for loci in current_DB that are not in geneace
  print "\nLooking for loci in /nfs/disk100/wormpub/DATABASES/current_DB that are not in $default_db\n" if ($verbose);
  &find_new_loci_in_current_DB($db);

   
  # Find sequences connected to multiple loci
  print "Looking for sequences connected to multiple loci\n" if ($verbose);
  &find_sequences_with_multiple_loci($db);

}

##############################

sub test_locus_for_errors{
  my $locus = shift;
  my $warnings;

  #test for Map AND !NEXT
  if($locus->at('Map')){
    my $map = $locus->Map;
    if (!defined($map)){
      $warnings .= "ERROR 1: $locus has a 'Map' tag but that tag has no map value!\n";
      print "." if ($verbose);
    }
  }

  # test for more than one Map
  if(defined($locus->at('Map'))){
    my @maps = $locus->at('Map');
    if(scalar(@maps) > 1){
      $warnings .= "ERROR 2: $locus has more than one Map object\n";
      print "." if ($verbose);
    }    
  }


  # test for more than one gene class
  if(defined($locus->at('Name.Gene_class'))){
    my @gene_classes = $locus->at('Name.Gene_Class');
    if(scalar(@gene_classes) > 1){
      $warnings .= "ERROR 3: $locus has more than one Gene_class\n";
      print "." if ($verbose);
    }    
  }

  
  # test for no Type tag (every valid Locus object should have a Type tag)
  if(!defined($locus->at('Type'))){  

    # Look for this Locus name in ?Gene_name class
    my ($gene_name) = $db->fetch(-class=>'Gene_name',-name=>"$locus");

    # Does this exist in ?Gene_name class, and is there an Other_name_for tag?
    if(($gene_name) && defined($gene_name->at('Gene.Other_name_for'))){
      my $other_gene = $gene_name->Other_name_for;
      $warnings .= "ERROR 4: $locus has no 'Type' tag but this Locus is listed as Other_name for $other_gene\n";
      print "." if ($verbose);
    }
    else{
      $warnings .= "ERROR 5: $locus has no 'Type' tag present and is not an Other_name for something else\n";
      print "." if ($verbose);
      if ($ace){
	# Basic info to promote into a proper 'Gene'
	print ACE "\n\nLocus : \"$locus\"\n";
	print ACE "Gene\n";
	print ACE "Non_CGC_name \"$locus\"\n";
      }
    }
  }


  # test for Gene AND !Species 
  if(defined($locus->at('Type.Gene')) && !defined($locus->at('Species'))){
    $warnings .= "ERROR 6: $locus has a 'Gene' tag but not a 'Species' tag\n";
    print "." if ($verbose);
    if ($ace){
      print ACE "\n\nLocus : \"$locus\"\n";
      if ($locus !~ /Cb-|Cr-|Cv-/){
	print ACE "Species \"Caenorhabditis elegans\"\n";
      }
      if ($locus =~/^Cb-.+/){
	print ACE "Species \"Caenorhabditis briggsae\"\n";
      }
      if ($locus =~/^Cr-.+/){
	print ACE "Species \"Caenorhabditis remanei\"\n";	
      }
      if ($locus =~/^Cv-.+/){
	print ACE "Species \"Caenorhabditis vulgaris\"\n";
      }	
    }
  }


  # test for Gene AND !CGC_name AND !Non_CGC_name
  if(defined($locus->at('Type.Gene')) && !defined($locus->at('Name.CGC_name'))
     && !defined($locus->at('Name.Non_CGC_name'))){					       
    $warnings .= "ERROR 7: $locus has a 'Gene' tag but not a 'CGC_name' or 'Non_CGC_name' tag\n";
    print "." if ($verbose);
  }


  # test for CGC_name AND Non_CGC_name
  if(defined($locus->at('Name.CGC_name')) && defined($locus->at('Name.Non_CGC_name'))){					       
    $warnings .= "ERROR 8: $locus has a 'CGC_name' tag *AND* a 'Non_CGC_name' tag\n";
    if ($ace){
      print ACE "\nLocus : \"$locus\"\n";
      print ACE "-D Non_CGC_name\n";
    }  
    print "." if ($verbose);
  }


  # test for CGC_name AND !CGC_approved
  if(defined($locus->at('Name.CGC_name')) && !defined($locus->at('Type.Gene.CGC_approved'))){		    
    $warnings .= "ERROR 9: $locus has a 'CGC_name' tag but not a 'CGC_approved' tag\n";
    if ($ace){
      print ACE "\nLocus : \"$locus\"\n";
      print ACE "CGC_approved\n";
    } 
    print "." if ($verbose);
  }

  # test for CGC_name AND CGC_unresolved
  if(defined($locus->at('Name.CGC_name')) && defined($locus->at('Type.Gene.CGC_unresolved'))){		    
    $warnings .= "ERROR 10: $locus has a 'CGC_name' tag but also has a 'CGC_unresolved' tag\n";
    if ($ace){
      print ACE "\nLocus : \"$locus\"\n";
      print ACE "-D CGC_unresolved\n";
    } 
    print "." if ($verbose);
  }

  # test for Non_CGC_name AND CGC_approved
  if(defined($locus->at('Name.Non_CGC_name')) && defined($locus->at('Type.Gene.CGC_approved'))){		    
    $warnings .= "ERROR 11: $locus has a 'Non_CGC_name' tag *AND* a 'CGC_approved' tag\n";
    if ($ace){
      print ACE "\nLocus : \"$locus\"\n";
      print ACE "-D CGC_approved\n";
    } 
    print "." if ($verbose);
  }

  # test for !Gene AND Gene_class 
  if(!defined($locus->at('Type.Gene')) && defined($locus->at('Name.Gene_class'))){
    $warnings .= "ERROR 12: $locus has a 'Gene_class' tag but not a 'Gene' tag\n";
    print "." if ($verbose);
    if ($ace){
      print ACE "\n\nLocus : \"$locus\"\n";
      print ACE "Gene\n";
    }
  }

  # test for CGC_approved but !Gene_class 
  if(defined($locus->at('Type.Gene.CGC_approved')) && !defined($locus->at('Name.Gene_class'))){
    $warnings .= "ERROR 13: $locus has a 'CGC_approved' tag but not a 'Gene_class' tag\n";
    print "." if ($verbose);
    my $gene_class = substr ($locus, 0, 3);
    if($ace){
      print ACE "\n\nLocus : \"$locus\"\n";
      print ACE "Gene_class $gene_class\n";
    }
  }
  

  # test for no Gene tag AND Genomic_sequence tag
  if(!defined($locus->at('Type.Gene')) && defined($locus->at('Molecular_information.Genomic_sequence'))){
    $warnings .= "ERROR 14: $locus has 'Genomic_sequence' tag but no 'Gene' tag\n";
    print "." if ($verbose);
    if ($ace){
      print ACE "\n\nLocus : \"$locus\"\n";  
      print ACE "Gene\n";
    }
  }


  # test for no CGC_approved tag AND Genomic_sequence/Transcript/Pseudogene tag
  if(!defined($locus->at('Type.Gene.CGC_approved'))){ 
    my $species = $locus->Species;
    if ($species eq "Caenorhabditis elegans"){
      if(defined($locus->at('Molecular_information.Genomic_sequence')) || 
	 (defined($locus->at('Molecular_information.Transcript'))) ||
	 (defined($locus->at('Molecular_information.Pseudogene')))){
	$warnings .= "ERROR 15: $locus has 'Genomic_sequence', 'Transcript', or 'Pseudogene' tag but no 'CGC_approved' tag\n" if $locus ne "arl-query";
	$warnings .= "ERROR 15: $locus has 'Genomic_sequence', 'Transcript', or 'Pseudogene' tag but no 'CGC_approved' tag\n" if $verbose;
	print JAHLOG "ERROR: $locus has 'Genomic_sequence', 'Transcript', or 'Pseudogene' tag but no 'CGC_approved' tag\n";
	$jah_errors++;
	print "." if ($verbose);
      }
    }
  }


  # test for Genomic_sequence tag but no value   
  if(defined($locus->at('Molecular_information.Genomic_sequence')) && !defined($locus->Genomic_sequence)){
    $warnings .= "ERROR 16: $locus has 'Genomic_sequence' tag but no associated value\n";
    print "." if ($verbose);
    if ($ace){
      print ACE "\n\nLocus : \"$locus\"\n"; 
      print ACE "-D Genomic_sequence\n";
    } 
  }

  # test for more than one Genomic_sequence tag, but need to allow for splice variants. A bit tricky this
  # and I think my RE (which just looks for word.number.letter) might allow some errors to creep through.  
  if(defined($locus->at('Molecular_information.Genomic_sequence'))){
    my @genomic_sequences = $locus->Genomic_sequence;

    if(scalar(@genomic_sequences)>1){
      my @problems = $locus->at('Molecular_information.Genomic_sequence');
      foreach my $problem (@problems){
	if ($problem !~ m/[\w\d]+\.\d+[a-z]/){
	  $warnings .= "ERROR 17: $locus has multiple 'Genomic_sequence' tags (see $problem)\n";
	  print "." if ($verbose);
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
	print CALTECHLOG "$seq has attached functional annotation which should now be attached to $locus\n";
	$caltech_errors++;
      }
    }
  }
  

  # test for Genomic_sequence AND Transcript tags both present
  if(defined($locus->at('Molecular_information.Genomic_sequence')) && defined($locus->at('Molecular_information.Transcript'))){  
    $warnings .= "ERROR 18: $locus has a 'Genomic_sequence' tag AND 'Transcript' tag\n";
    print "." if ($verbose);
  }


  # test for Genomic_sequence AND Pseudogene tags both present
  if(defined($locus->at('Molecular_information.Genomic_sequence')) && defined($locus->at('Molecular_information.Pseudogene'))){  
    $warnings .= "ERROR 19: $locus has a 'Genomic_sequence' tag AND 'Pseudogene' tag\n";
    print "." if ($verbose);
  }


  # test for Pseudogene AND Transcript tags both present
  if(defined($locus->at('Molecular_information.Pseudogene')) && defined($locus->at('Molecular_information.Transcript'))){  
    $warnings .= "ERROR 20: $locus has a 'Pseudogene' tag AND 'Transcript' tag\n";
    print "." if ($verbose);
  }


  # test for Genomic_sequence AND !Sequence_name
  if(defined($locus->at('Molecular_information.Genomic_sequence')) && !defined($locus->at('Name.Sequence_name'))){  
    my $seq = $locus->Genomic_sequence;
    $warnings .= "ERROR 21: $locus has a 'Genomic_sequence' tag but not a  'Sequence_name' tag\n";
    print "." if ($verbose);
    if ($ace){
      print ACE "\n\nLocus : \"$locus\"\n"; 
      print ACE "Sequence_name $seq\n";
    } 
  }



 # test for Genomic_sequence = Sequence_name
  if(defined($locus->at('Molecular_information.Genomic_sequence')) && defined($locus->at('Name.Sequence_name'))){  
    my @seq = $locus->Genomic_sequence;
    my @name = $locus->Sequence_name;

    # compare @seq with @name for difference (order does not matter)
    my @comp_result = array_comp(\@seq, \@name);
    my @diff = @{$comp_result[0]}; 
    
    if (@diff){
      $warnings .= "ERROR 22: $locus has a 'Genomic_sequence' name different to 'Sequence_name'\n";
      print "." if ($verbose);
    }
  }


  # test for Transcript AND !Transcript_name
  if(defined($locus->at('Molecular_information.Transcript')) && !defined($locus->at('Name.Transcript_name'))){  
    my $seq = $locus->Transcript;
    $warnings .= "ERROR 23: $locus has a 'Transcript' tag but not a 'Transcript_name' tag\n";
    print "." if ($verbose);
    if ($ace){
      print ACE "\n\nLocus : \"$locus\"\n"; 
      print ACE "Transcript_name $seq\n";
    } 
  }


  # test for Transcript = Transcript_name
  if(defined($locus->at('Molecular_information.Transcript')) && defined($locus->at('Name.Transcript_name'))){  
    my @seq = $locus->Transcript;
    my @name = $locus->Transcript_name;

    # compare @seq with @name for difference (order does not matter)
    my @comp_result = array_comp(\@seq, \@name);
    my @diff = @{$comp_result[0]}; 
    
    if (@diff){ 
      $warnings .= "ERROR 24: $locus has a 'Transcript' name different to 'Transcript_name'\n";
      print "." if ($verbose);
    }
  }



  # test for Pseudogene AND !Pseudogene_name
  if(defined($locus->at('Molecular_information.Pseudogene')) && !defined($locus->at('Name.Pseudogene_name'))){  
    my $seq = $locus->Pseudogene;
    $warnings .= "ERROR 25: $locus has a 'Pseudogene' tag but not a 'Pseudogene_name' tag\n";
    print "." if ($verbose);
    if ($ace){
      print ACE "\n\nLocus : \"$locus\"\n"; 
      print ACE "Pseudogene_name $seq\n";
    } 
  }


  # test for Pseudogene = Pseudogene_name
  if(defined($locus->at('Molecular_information.Pseudogene')) && defined($locus->at('Name.Pseudogene_name'))){  
    my @seq = $locus->Pseudogene;
    my @name = $locus->Pseudogene_name;

    # compare @seq with @name for difference (order does not matter)
    my @comp_result = array_comp(\@seq, \@name);
    my @diff = @{$comp_result[0]}; 
    
    if (@diff){
      $warnings .= "ERROR 26: $locus has a 'Pseudogene' name different to 'Pseudogene_name'\n";
      print "." if ($verbose);
    }
  }


  # test for Polymorphism AND Gene tags both present
  if(defined($locus->at('Type.Gene')) && defined($locus->at('Type.Polymorphism'))){  
    $warnings .= "ERROR 27: $locus has a 'Polymorphism' tag AND 'Gene' tag\n";
    print "." if ($verbose);
  }


  # Look for Gene_class tag in non-gene objects 
  if(!defined($locus->at('Type.Gene'))){
    if(defined($locus->at('Name.Gene_class'))){
      $warnings .= "ERROR 28: $locus has Gene_class tag but it is not a gene!\n";
      print "." if ($verbose);
      if ($ace){
        print ACE "\n\nLocus : \"$locus\"\n"; 
        print ACE "Gene\n";
      } 
    }
  }
  

  # test for Other_name tag but no value   
  if(defined($locus->at('Name.Other_name')) && !defined($locus->Other_name)){
    $warnings .= "ERROR 29: $locus has 'Other_name' tag but no associated value\n";
    print "." if ($verbose);
  }


  # Remind of outstanding CGC_unresolved tags
  if(defined($locus->CGC_unresolved)){
    my ($unresolved_details) = $locus->at('Type.Gene.CGC_unresolved');
    $warnings .= "ERROR 30: $locus has CGC_unresolved tag: \"$unresolved_details\"\n" if $locus ne "arl-qury";
    $warnings .= "ERROR 30: $locus has CGC_unresolved tag: \"$unresolved_details\"\n" if $verbose;
  }


  # Look for loci with no Positive_clone info but with information that could be added from
  # Genomic_sequence, Transcript, or Pseudogene tags
  if(!defined($locus->at('Positive.Positive_clone'))){
    my $seq;
    ($seq = $locus->Genomic_sequence) if(defined($locus->at('Molecular_information.Genomic_sequence')));
    ($seq = $locus->Transcript) if(defined($locus->at('Molecular_information.Transcript')));
    ($seq = $locus->Pseudogene) if(defined($locus->at('Molecular_information.Pseudogene')));

    # If there is a sequence/transcript/pseudogene
    if($seq){
      # need to chop off the ending to just get clone part
      $seq =~ s/\..*//;
      $warnings .= "ERROR 31: $locus has no 'Positive_clone' tag but info is available from Genomic_sequence, Transcript, Pseudogene tag\n";
      print "." if ($verbose);
      # must use Inferred_automatically from Evidence hash for this type of info
      print ACE "\n\nLocus : \"$locus\"\nPositive_clone \"$seq\" Inferred_automatically \"From sequence, transcript, pseudogene data\"\n" if ($ace);
    }
  }


  # test for Other_name value which is also a Locus name in its own right
  # Also can test for where a Locus has an other name but this fact hasn't been added to parent Gene_class object
  # as a Remark
  if(defined($locus->at('Name.Other_name'))){
    
    foreach (@exceptions){$exceptions{$_}++};  # @exceptions and %exceptions are defined as global

    # get other names of locus
    my @other_names = $locus->Other_name;
    foreach my $other_name (@other_names){
      
      # Does other_name exist as separate loci?
      if($db->fetch(-class=>'Locus',-name=>"$other_name")){
	# Is it on exceptions list?
	if($exceptions{$locus}) {
	  # no need for warning
	  if ($verbose){
	    print LOG "INFO: $locus has $other_name as Other_name...$other_name is still a separate Locus object (exception)\n";
	  }  
	}
	else{
	  # can warn that this is potential problem
	  print LOG "WARNING: $locus has $other_name as Other_name...$other_name is still a separate Locus object\n";
	  print "." if ($verbose);
	  if ($ace){
	    print ACE "\n-R Locus : \"$other_name\" \"$locus\"\n";
	    print ACE "\nLocus : \"$locus\"\n";
	    print ACE "Other_name \"$other_name\"\n";
	  }
	}	
	# Does the Parent Gene_class have a Remark for this Other_name?
	my $gene_class = $locus->Gene_class;	
	my ($obj) = $db->fetch(-class=>'Gene_class',-name=>"$gene_class") if(defined($locus->at('Name.Gene_class')));
	if(!defined($obj->at('Remark'))){
	  $warnings .= "ERROR 32: $locus has Other_name but no Remark in Gene_class object\n";
	  print "." if ($verbose);      
	  if($ace){
	    print ACE "\n\nGene_Class : \"$gene_class\"\n";
	    print ACE "Remark \"$locus is also an unofficial other name of $other_name\" CGC_data_submission\n";
	    
	  }         
	}
      }
    }
  }


  return($warnings);

}


############################################

sub find_new_loci_in_current_DB{
  my $db = shift;
  my $warnings;

  # open a database connection to current_DB and grab all loci names (excluding polymorphisms)
  my $new_db = Ace->connect(-path  => "$curr_db",
		    -program =>$tace) || do { print LOG "Connection failure: ",Ace->error; die();};
  my @current_DB_loci = $db->fetch(-query=>'Find Locus;!Polymorphism');
  $new_db->close;


  #cross reference current_DB loci against those in geneace
  foreach my $loci(@current_DB_loci){
    my $new_loci = $db->fetch(-class=>'Locus',-name=>"$loci");
    unless(defined($new_loci)){
      $warnings .= "ERROR 33: $new_loci in current_DB is not in $default_db\n";
    }
  }
  print LOG "\n$warnings\n" if $warnings;

}


#############################

sub find_sequences_with_multiple_loci {
  my $db = shift;

  # Find sequences that have more than one locus_genomic_seq tag
  my @sequences = $db->fetch(-query=>'Find Sequence COUNT Locus_genomic_seq >1');
  foreach my $seq (@sequences){
    my @loci = $seq->Locus_genomic_seq;

    print LOG "ERROR 34: Multiple loci (@loci) are connected to the same sequence ($seq)\n";
    print JAHLOG "ERROR: Multiple loci (@loci) are connected to the same sequence ($seq)\n";

    $jah_errors++;
  }
} 


##############################################################################################################################


#######################################
# Check misuse of Evidence in 8 classes 
# Convert Author to Person
#######################################

sub check_evidence {


my $WBPerson_F_M_L_names=<<EOF;
Table-maker -p "$def_dir/WBperson_first_middle_last_names.def" quit
EOF

  my @WBPerson = &process_WBPerson_names($WBPerson_F_M_L_names, $curr_db);
 
  print LOG "\n\nChecking misuse of Evidence and converting Author to Person / Non-Person to Author:\n";
  print LOG "-----------------------------------------------------------------------------------\n\n";

# dump flat files with time stamps

my $command=<<END;
find locus * 
show -a -T -f /tmp/locus_dump.ace

find allele *
show -a -T -f /tmp/allele_dump.ace

find strain *
show -a -T -f /tmp/strain_dump.ace

find gene_class *
show -a -T -f /tmp/geneclass_dump.ace

find 2_point_data *
show -a -T -f /tmp/2_pt_dump.ace

find Multi_pt_data *
show -a -T -f /tmp/multi_pt_dump.ace

find Pos_neg_data *
show -a -T -f /tmp/posneg_dump.ace

find Laboratory *
show -a -T -f /tmp/lab_dump.ace

quit
END

  open (DUMP, "| $tace $default_db") || die "Failed to connect to Geneace";
  print DUMP $command;
  close DUMP;

  system ("cat /tmp/locus_dump.ace /tmp/allele_dump.ace /tmp/strain_dump.ace /tmp/2_pt_dump.ace /tmp/multi_pt_dump.ace /tmp/posneg_dump.ace /tmp/geneclass_dump.ace /tmp/lab_dump.ace> /tmp/class_dump.ace");

  my @dumps = qw (locus_dump.ace allele_dump.ace strain_dump.ace geneclass_dump.ace 2_pt_dump.ace multi_pt_dump.ace posneg_dump.ace class_dump.ace lab_dump.ace);  
 
  foreach (@dumps){system ("chmod 777 /tmp/$_")}

  open(IN, "/tmp/class_dump.ace") || die $!;

# look for person/author names that needs to be converted in 8 classes (regardless of which tag as the scripts greps from string pattern

  my $evid_errors = 0;
  my $updates = 0;
  my $info_num = 0;
  my (@counters, $class_obj, $class, $obj, $tag, $ori, $b4_evi, $name, $paper, $author, $last_name, $initials);

  while (<IN>){
    chomp;
    if ($_ =~ /^(Locus) : \"(.+)\" -O .+/){$class = $1; $obj = $2}
    if ($_ =~ /^(Allele) : \"(.+)\" -O .+/){$class = $1; $obj = $2}
    if ($_ =~ /^(Strain) : \"(.+)\" -O .+/){$class = $1; $obj = $2}
    if ($_ =~ /^(Gene_Class) : \"(.+)\" -O .+/){$class = $1; $obj = $2}
    if ($_ =~ /^(2_point_data) : \"(.+)\" -O .+/){$class = $1; $obj = $2}     
    if ($_ =~ /^(Multi_pt_data) : \"(.+)\" -O .+/){$class = $1; $obj = $2}  
    if ($_ =~ /^(Pos_neg_data) : \"(.+)\" -O .+/){$class = $1; $obj = $2}  
    if ($_ =~ /^(Laboratory) : \"(.+)\" -O .+/){$class = $1; $obj = $2}  

    if ($_ =~ /((\w+)\s+.+)Person_evidence -O .+\"(.+)\" -O.+/){
      $ori = $_;
      $b4_evi = $1;	
      $tag = $2;
      $name = $3;
      if ($name !~ /WBPerson\d+/){
	#$evid_errors++;
	print LOG "\nERROR: $class $obj has non-Person $name under main tag $tag\n";

	@counters = get_WBPerson_ID($name, $class, $obj, $tag, $ori, $b4_evi, "PtoA");  
	$evid_errors += $counters[0];
	$updates += $counters[1];
	$info_num += $counters[2]; 
      }  
    }
    if ($_ =~ /((\w+)\s+.+)Paper_evidence -O .+\"(.+)\" -O.+/){
      $ori = $_;
      $b4_evi = $1;
      $tag = $2;
      $paper = $3;
      $class_obj = $class." : "."\"$obj\"";
      if ($paper !~ /\[.+\]/){
        $evid_errors++;
	print LOG "\nERROR: $class $obj has Paper $paper under main tag $tag\n";
	if ($ace){
	  print ACE "\n$class_obj\n";
	  print ACE "-D $ori\n\n";
	  print ACE "\n$class_obj\n";
	  $paper =~ s/\[|\]//g;
	  print ACE "$b4_evi Paper_evidence \"\[$paper$\]\"\n";
        }
      }
     # Likely to be deleted, when Caltech is clear about rules of using cgc and pmid
     # if ($paper !~ /\[cgc\d+\]/ && $paper =~ /\[pmid(\d+)\]/){
#	$evid_errors++;
#	$paper = $1;
#	print LOG "\n: $class $obj has Paper $paper under main tag $tag\n";
#	if ($ace){
#	  print ACE "\n$class_obj\n";
#	  print ACE "-D $ori\n\n";
#	  print ACE "\n$class_obj\n";
#	  print ACE "$b4_evi PMID_evidence \"\[$paper$\]\"\n"; 
#        }
#      }
    }  

    if ($_ =~ /((\w+)\s+.+)Author_evidence -O .+\"(.+)\" -O.+/){
      $ori = $_;
      $b4_evi = $1;
      $tag = $2;
      $author = $3;

      @counters = get_WBPerson_ID($author, $class, $obj, $tag, $ori, $b4_evi); 
      $evid_errors += $counters[0];
      $updates += $counters[1];
      $info_num += $counters[2]; 
    }
    if ($_ =~ /(Unregistered_lab_members)\s+-O\s\"\d+.+_\w+\" \"(.+)\" -O .+/){
    
      $ori = $_;
      $tag = $1;
      $author = $2;
      #print $author, "#\n";
      # print "$author, $class, $obj, $tag, $ori@\n";
      @counters = get_WBPerson_ID($author, $class, $obj, $tag, $ori); 
      $evid_errors += $counters[0];
      $updates += $counters[1];
      $info_num += $counters[2]; 
    }
    
  }
  print LOG "\n\nThere are $evid_errors Evidence errors in 8 classes checked\n";
  print LOG "\n$updates Authors can be converted to Persons\n";
  print LOG "\n$info_num Authors are not Persons\n" if $verbose;
  system ("rm -f /tmp/*_dump.ace");

  #################################################################
  # subroutines for checking evidence & converting Author to Person
  #################################################################
  
  # make hashes for last name and its corresponding first/middle name and WBPersonID
  # these hashes are used by get_WBPerson_ID subroutine

  sub process_WBPerson_names {
    my ($def, $db)=@_;
    my ($WBPerson, $F_name, $M_name, $L_name, $F_char, $M_char);
    open (FH, "echo '$def' | $tace $db | ") || die "Couldn't access current_DB\n";
    while (<FH>){
      chomp($_);
      if ($_ =~ /^\"(WBPerson\d+)\"\s+\"(\w+)\"\s+\"(\w+|\w+.)\"\s+\"(\w+|\w+-\w+)\"$/){ 
	$WBPerson = $1;
	$F_name = $2;
	$M_name = $3;
	$L_name = $4;
	$F_char = substr ($F_name, 0, 1);
	$M_char = substr ($M_name, 0, 1);
	push (@{$L_name_F_M_WBP{$L_name}}, $F_char.$M_char, $F_char, $WBPerson);
      }
      if ($_ =~ /^\"(WBPerson\d+)\"\s+\"(\w+)\"\s+\"(\w+|\w+-\w+)\"$/){ 
	$WBPerson = $1;
	$F_name = $2;
	$L_name = $3;
	$F_char = substr ($F_name, 0, 1);
	push (@{$L_name_F_WBP{$L_name}}, $F_char, $WBPerson);
      } 
      
    }
    close FH;
  }

  # last name (key) is checked against hashes created from process_WBPerson_names subroutine
  # the values of the key (first name, middle name) are checked to assign identify
  
  # convert authors that have WBPersonID to person 
  # (or unregistered_lab_members to registered_lab_members for laboratory class)
  
  # move names in person_evidence that are not person to author_evidence 

  sub get_WBPerson_ID {
    my ($name, $class, $obj, $tag, $ori, $b4_evi, $conversion) = @_;
    my ($last_name, $initials, $num);
    my $class_obj = $class." : "."\"$obj\"";
    my $convert = 0;
    my $info_count = 0;
    my $evid_errors = 0;

    if ($name =~ /\w+\,\w+/ || $name =~ /\w+\,\s+\w+/){
      $evid_errors++;
      if (!defined $conversion){
	($last_name, $initials) = split(/ /, $name);
	$last_name =~ s/,//;
	print LOG "\nERROR: $class $obj has $name (Author) under main tag $tag\n";
	print LOG"=====>Correct $name as $last_name $initials\n";
	if ($ace){
	  print ACE "\n$class_obj\n";
	  print ACE "-D $ori\n\n";
	  print ACE "\n$class_obj\n";
	  if (defined $b4_evi){
	    print ACE "$b4_evi Author_evidence \"$last_name $initials\"\n";
	  }
	  else {
	    print ACE "$tag \"$last_name $initials\"\n";
	  }
	}
      }
    }
    else {
      ($last_name, $initials) = split(/ /, $name);
    }
    
    if (!exists $L_name_F_WBP{$last_name} && !exists $L_name_F_M_WBP{$last_name}){
      if (defined $conversion){
	print LOG "=====>Move $name under Author_evidence as NO corresponding WBPersonID exists\n";
      }
      if ($verbose && !defined $conversion){
	$info_count++; 
	print LOG "\nINFO: $class $obj has $name (Author under $tag tag): NOT yet a Person\n";
      }
      if ($ace && defined $conversion){
	print ACE "\n$class_obj\n";
	print ACE "-D $ori\n\n";
	print ACE "\n$class_obj\n";
	if (defined $b4_evi){
	  print ACE "$b4_evi Author_evidence \"$last_name $initials\"\n";
	}
      }
    }

    if (exists $L_name_F_WBP{$last_name}){
      $num = scalar @{$L_name_F_WBP{$last_name}};
      
      for (my $i=0; $i< $num; $i=$i+2){
	if ($initials eq @{$L_name_F_WBP{$last_name}}->[$i]){
	  $convert++;
	  if (!defined $conversion){
	    print LOG "\nUPDT: $class $obj has $name (Author) under $tag tag\n";
          }
	  if ($num == 2){
	    print LOG "=====>$name can now be Person @{$L_name_F_WBP{$last_name}}->[$i+1]\n";   
	    if ($ace){
              print ACE "\n$class_obj\n";
              print ACE "-D $ori\n\n";
	      print ACE "\n$class_obj\n";
	      if (defined $b4_evi){
		print ACE "$b4_evi"." Person_evidence \"@{$L_name_F_WBP{$last_name}}->[$i+1]\"\n";
	      }
	      else {
		print ACE "Registered_lab_members \"@{$L_name_F_WBP{$last_name}}->[$i+1]\"\n";
	      }	
            }
	  }
	  else {
	    print LOG "=====>$name might be Person @{$L_name_F_WBP{$last_name}}->[$i+1]\n"; 
	    if ($ace){
              print ACE "\n$class_obj\n";
              print ACE "-D $ori\n\n";
	      print ACE "\n$class_obj\n";
	      if (defined $b4_evi){
		print ACE "$b4_evi"." Person_evidence \"@{$L_name_F_WBP{$last_name}}->[$i+1]\"\n";
	      }
	      else {
		print ACE "Registered_lab_members \"@{$L_name_F_WBP{$last_name}}->[$i+1]\"\n";
	      }	
            }
	  }
	}
      }
    }

    if (exists $L_name_F_M_WBP{$last_name}){
      $num = scalar @{$L_name_F_M_WBP{$last_name}};
      
      for (my $i=0; $i< $num; $i=$i+3){
	if ($initials eq @{$L_name_F_M_WBP{$last_name}}->[$i] || 
            $initials eq @{$L_name_F_M_WBP{$last_name}}->[$i+1] ){
          $convert++;
          if (!defined $conversion){
	    print LOG "\nUPDT: $class $obj has $name (Author) under $tag tag\n";
          }
       	  if ($num == 3){
	    print LOG "=====>$name can now be Person @{$L_name_F_M_WBP{$last_name}}->[$i+2]\n";
            if ($ace){
              print ACE "\n$class_obj\n";
              print ACE "-D $ori\n\n";
	      print ACE "\n$class_obj\n";
	      if (defined $b4_evi){
		print ACE "$b4_evi"." Person_evidence \"@{$L_name_F_M_WBP{$last_name}}->[$i+2]\"\n";
	      }
	      else {
		print ACE "Registered_lab_members \"@{$L_name_F_M_WBP{$last_name}}->[$i+2]\"\n";
	      } 
            }
	  }
	  else {
	    print LOG "=====>$name might be Person @{$L_name_F_M_WBP{$last_name}}->[$i+2]\n";
            if ($ace){
              print ACE "\n$class_obj\n";
              print ACE "-D $ori\n\n";
	      print ACE "\n$class_obj\n";
	      if (defined $b4_evi){
		print ACE "$b4_evi"." Person_evidence \"@{$L_name_F_M_WBP{$last_name}}->[$i+2]\"\n";
	      }
	      else {
		print ACE "Registered_lab_members \"@{$L_name_F_M_WBP{$last_name}}->[$i+2]\"\n";
	      } 
            }
	  }
	}
      }
    } 
    $conversion =();
    return $evid_errors, $convert, $info_count;
  }
}	


##############################################################################################################################


sub process_laboratory_class{

  print "\n\nChecking Laboratory class for errors:\n";
  print LOG "\n\nChecking Laboratory class for errors:\n";
  print LOG "-------------------------------------\n";
  #grab lab details
  my @labs = $db->fetch(-class => 'Laboratory',
		        -name  => '*');

  # test for Allele_designation and Representative tags
  foreach my $lab (@labs){
    if($verbose && !defined($lab->at('CGC.Allele_designation')) && $lab =~ /\w{3,3}/){  
      print LOG "INFO: $lab has no Allele_designation tag present (exception)\n";
    }    
    
    if(!defined($lab->at('CGC.Representative')) && $lab ne "CGC"){  
      if ($lab ne "XA"){
	print LOG "WARNING: $lab has no Representative tag present\n";
      }
    }  
    undef($lab);
  }
}


##############################################################################################################################
 

sub process_allele_class{

  print"\n\nChecking Allele class for errors:\n";
  print LOG "\n\nChecking Allele class for errors:\n";
  print LOG "---------------------------------\n";

  my @alleles = $db->fetch('Allele','*');
  my ($allele, $desig, $desig2, $main_allele, %allele_gene, $gene);
  my $allele_error = 0;

  my $allele_designation_to_LAB=<<EOF;
  Table-maker -p "$def_dir/allele_designation_to_LAB.def" quit
EOF
  
  my %location;
  # get a hash of allele-desig (key) and lab-desig (value)
  if ($ace){%location=allele_location($allele_designation_to_LAB, $default_db)};
  
  foreach $allele (@alleles){
   
    # check allele has not location tag
    if(!defined($allele->at('Location'))){
      # catch non-standard upper-case allele name
      if ($allele =~ /^[A-Z].+/){ 
	$allele_error++;
	print LOG "ERROR: $allele has no Location tag present (no info available)\n";
      }
      else {
	$allele_error++;
	print LOG "ERROR: $allele has no Location tag present\n";

	# acefile output for normal or double alleleles having no location tag
	if ($ace){
	  if ($allele =~ /^([a-z]{1,})\d+$/){
	    $desig = $1;
	    print  ACE "\n\nAllele : \"$allele\"\n";
	    print  ACE "Location \"$location{$desig}\"\n";
	    next;
	  }

	  ######################
	  # catch double alleles
          ######################

	  # double allele has same designation
	  if ($allele =~ /^([a-z]{1,})\d+([a-z]{1,})\d+$/){
	    $desig = $1; $desig2 = $2;
	    if ("$desig" eq "$desig2"){
	      print  ACE "\n\nAllele : \"$allele\"\n";
	      print  ACE "Location \"$location{$desig}\"\n";
	    }
	    else {

	      # double allele has diff designation
	      $desig = $1; $desig2 = $2;
	      print  ACE "\n\nAllele : \"$allele\"\n";
	      print  ACE "Location \"$location{$desig}\"\n";
	      print  ACE "Location \"$location{$desig2}\"\n";
	    }
	    next;
          }
	}
      }
    }

    if($allele -> Gene){
      my @loci=$allele->Gene(1);
      if (scalar @loci > 1){
	$allele_error++;
	print LOG "ERROR: $allele is connected to more than one Loci: @loci\n";
      }
    }
  }
  
  my $allele_has_flankSeq_and_no_seq=<<EOF;
  Table-maker -p "$def_dir/allele_has_flankSeq_and_no_seq.def" quit
EOF
  
  &allele_has_flankSeq_and_no_seq($allele_has_flankSeq_and_no_seq, $default_db);
  
  my $allele_has_predicted_gene_and_no_seq=<<EOF;
  Table-maker -p "$def_dir/allele_has_predicted_gene_and_no_seq.def" quit
EOF
  
  &allele_has_predicted_gene_and_no_seq($allele_has_predicted_gene_and_no_seq, $default_db);
  
  my $allele_has_no_methods=<<EOF;
  Table-maker -p "$def_dir/allele_methods.def" quit
EOF
  
  check_missing_allele_method($allele_has_no_methods, $default_db);

  print LOG "\nThere are $allele_error errors in Allele class\n";

  ##################################################
  # subroutines for process_allele_class routine
  ##################################################
  
  sub allele_location {
    my ($def, $dir)=@_;
    my %location_desig;
    open (FH, "echo '$def' | $tace $dir | ") || die "Couldn't access geneace\n";
    while (<FH>){
      chomp($_);
      if ($_ =~ /^\"(.+)\"\s+\"(.+)\"/){
	$location_desig{$2} = $1  # $2 is allele_designation $1 is LAB	
      }
    }
    return %location_desig;
  }
  
  sub allele_has_flankSeq_and_no_seq {
    
    my ($def, $dir) = @_;
    open (FH, "echo '$def' | $tace $dir | ") || die "Couldn't access geneace\n";
    while (<FH>){
      chomp($_);
      if ($_ =~ /^\"/){
	$_ =~ s/\"//g;
	$allele_error++;
	print LOG "ERROR: Allele $_ has flanking sequences but has no parent sequence\n";
      }
    }
  }
  

  sub allele_has_predicted_gene_and_no_seq {
          
    my ($def, $dir) = @_;
    my ($allele, $seq, $parent, $cds);
    
    open (FH, "echo '$def' | $tace $dir | ") || die "Couldn't access geneace\n";
    while (<FH>){
      chomp($_);
      if ($_ =~ /^\"(.+)\"\s+\"(.+)\"\s$/) {
	$allele = $1;
	$cds = $2;
	$allele_error++;
	print LOG  "ERROR: Allele $allele has predicted gene but has no parent sequence\n";
	if ($ace){
	  get_parent_seq($cds, $allele);
	}
      }
      if ($_ =~ /\"(.+)\"\s+\"(.+)\"\s+\"(.+)\"/){
	$allele = $1;  
	$cds = $2;
	$seq = $3;
	if ($seq eq $cds){
	  $allele_error++;
	  # temporarily commented out to see if suitable to automate this part
	  print LOG "ERROR: $allele has a different parent sequence ($seq) based on its predicted gene ($cds) (overlapped clones?)\n";
	  if ($ace){
	    #print ACE "\n\nAllele : \"$allele\"\n";
	    #print ACE "-D Sequence \"$seq\"\n";
	    &get_parent_seq($cds, $allele);
	  }
	}
	if ($seq ne $cds && $seq !~ /SUPERLINK.+/){
	  $parent=get_parent_seq($cds, $allele, "getparent");
	  if ($parent ne $seq){
	    $allele_error++;
	    print LOG "ERROR: $allele has a different parent sequence ($seq) based on its predicted gene ($cds) (overlapped clones?)\n";
	    if ($ace){
	      # temporarily commented out to see if suitable to automate this part
	      #print ACE "\nAllele : \"$allele\"\n";
	      #print ACE "-D Sequence \"$seq\"\n";
	      #print ACE "Sequence \"$parent\"\n";
	    }
	  } 
	}
      }
    }
    sub get_parent_seq {
      my ($predict, $allele, $get_parent) = @_;   
      my ($parent, $cds);
      if ($predict =~ /(.+)\.(\d+)[a-z]/ || $predict =~ /(.+)\.(\d+)/){
	$parent =  $1;
	if (!$get_parent){
	  print ACE "\n\nAllele : \"$allele\"\n";
	  print ACE "Sequence \"$parent\"\n";
	}
      }
      return $parent;
    }
  }
  
  sub check_missing_allele_method {
    my ($def, $dir) = @_;
    my ($allele, $tag);
    
    open (FH, "echo '$def' | $tace $dir | ") || die "Couldn't access geneace\n";
    while (<FH>){
      chomp($_);
      print $_, "\n";
      if ($_ =~ /^\"(.+)\"\s+\"(.+)\"/){
	$allele = $1;
	$tag = $2;
	# check type of method to use
	if ($tag eq "KO_consortium_allele"){$tag = "Knockout_allele"}
	if ($tag eq "Transposon_insertion"){$tag = "Transposon_insertion"}
	if ($ace){output($allele, $tag, "ace")}
	else {output($allele, $tag)}		    
      }
      if ($_ =~ /^\"(.+)\"/){
	print $1, "\n";
	$allele = $1;
	if ($ace){output($allele, "Allele", "ace")}
	else {$allele_error++; print LOG "ERROR: Allele $allele has no Method \"Allele\"\n"}		    
      }
    }
    sub output {
      my ($allele, $tag, $ace) = @_;
      $allele_error++;
      print LOG "ERROR: Allele $allele has no Method $tag\n";
      if ($ace ne ""){
	# output acefile for correct type of method of an allele
	print ACE "\n\nAllele : \"$allele\"\n";
	print ACE "Method \"$tag\"\n";
      }
    }
  }
}

##############################################################################################################################


sub process_strain_class {

  # Check if sequence name of a strain genotype is now connected to a locus

  my (@strains, @seqs, $genotype, @genes, $seq, $strain, $extract);

  print"\n\nChecking Strain class for errors:\n";
  print LOG "\n\nChecking Strain class for errors:\n";
  print LOG "---------------------------------\n";

  @strains = $db->fetch('Strain','*');
  @seqs = $db->fetch('Sequence','*');

  my %seqs;
  foreach (@seqs){
    $seqs{$_}++;
  }

  my $cgc_approved_loci=<<EOF;
  Table-maker -p "$def_dir/cgc_approved_loci.def" quit
EOF
  
  my %cgc_loci=cgc_loci($cgc_approved_loci, $default_db);
  my $e;

  foreach $strain (@strains){
    if (!$strain->Location){
      print LOG "WARNING: Strain $strain has no location tag\n";
      if ($ace){
	$strain =~ /([A-Z]+)\d+/;
	print ACE "\n\nStrain : \"$strain\"\n";
	print ACE "Location \"$1\"\n";
      }
    }
    else { 
      my $cgc=$strain->Location;
      if ($strain->Genotype){
	$genotype = $strain->Genotype(1);
	print "$genotype\n" if ($strain eq "RB992");
	$extract = $genotype;
	$extract =~ s/\(|\)|\/|\+|;|\?|\{|\}|,|\=|\.$/ /g;
	$extract =~ s/ I | II | III | IV | V | X / /g;
	$extract =~ s/ I | II | III | IV | V | X | f / /g; # further tidying up of chromosomes
	$extract =~ s/^\s|\w{3}-\s| f | A //g;
	$extract =~ s/\s{1,}/ /g;
	@genes=split(/ /,$extract);
	print @genes, "\n" if ($strain eq "C52E12.2");
	foreach (@genes){
	  if($seqs{$_}){
	    my $seq = $db->fetch('Sequence', $_);
	    if ($seq->Locus_genomic_seq){
	      my @loci=$seq->Locus_genomic_seq(1);
	      if ($cgc eq "CGC"){
		foreach $e (@loci){
		  if ($cgc_loci{$e}){
		    print LOG "WARNING: CGC Strain $strain has sequence_name $_ in Genotype, which can now become $e\n";
		  }
                }  
	      }
	      else {
		foreach $e (@loci){
		  if ($cgc_loci{$e}){  
		    print LOG "WARNING: Non_CGC Strain $strain has sequence_name $_ in Genotype, which can now become $e\n";
                  }
		}  
	      }  
	    }
	  }
	}
      }
    }
  }
  
  my ($locus, %allele_locus, %strain_genotype, $cds, %locus_cds, $main, $other_name, %other_main, $allele);

  my $get_genotype_in_strain=<<EOF;
  Table-maker -p "$def_dir/strain_genotype.def" quit 
EOF
  my $allele_to_locus=<<EOF;
Table-maker -p "$def_dir/allele_to_locus.def" quit
EOF
    

  open (FH1, "echo '$get_genotype_in_strain' | $tace $default_db | ") || die $!;
  open (FH5, "echo '$allele_to_locus' | $tace $default_db | ") || die $!;

  while(<FH1>){
    chomp;
 
     if ($_ =~ /\"(.+)\"\s+\"(.+)\"/){
      $strain_genotype{$1} = $2;
    }
  }   
  while (<FH5>){
    chomp $_;
    if ($_ =~ /\"(.+)\"\s+\"(.+)\"/){
      $locus = $2;
      $locus =~ s/\\//;
      if (!$exceptions{$locus}){
	push(@{$allele_locus{$1}}, $locus);
      }
    }
  } 

  foreach my $strain (keys %strain_genotype){
    my @matches = ($strain_genotype{$strain} =~  /((Cb-\w{3,3}-\d+|Cr-\w{3,3}-\d+|\w{3,3}-\d+)\(\w+\d+\))/g);
    foreach (@matches){
      my @la = split(/\s+/, $_);
      foreach (@la){
	if ($_ =~ /(Cb-\w{3,3}-\d+|Cr-\w{3,3}-\d+|\w{3,3}-\d+)\((\w+\d+)\)/){
	  $locus = $1; $allele = $2; 
	  
	  # diff allele->locus link in geneace and allele->locus in strain genotype
	  # if diff, the print error LOG
	  if ((defined @{$allele_locus{$allele}}) && ("@{$allele_locus{$allele}}" ne "$locus")){
  	    print LOG "ERROR: Strain $strain has $locus($allele) in genotype: ";
	    print LOG "change each $locus to @{$allele_locus{$allele}}\n";
          }
	}
      }
    }
  }
}

#################################################################

sub cgc_loci {
  my ($def, $db) = @_;
  my (@cgc_loci, %cgc_loci);
  open (FH, "echo '$def' | $tace $db | ") || die "Couldn't access $db\n"; 
   
  while (<FH>){
    chomp $_;
    if ($_ =~ /^\"(.+)\"/){
      push(@cgc_loci, $1);
    }
  }
  foreach (@cgc_loci){$cgc_loci{$_}++}
  return %cgc_loci;
} 

###############################
# Process Rearrangement class #
###############################

sub process_rearrangement {
 
  print"\n\nChecking Rearrangement class for errors:\n";
  print "----------------------------------------\n";	
  print LOG "\n\nChecking Rearrangement class for errors:\n";
  print LOG "----------------------------------------\n";
  # checks presence of non-rearrangement object 
  # as objects of Rearrangement class

  my @rearr;
 
  @rearr = $db -> fetch('Rearrangement','*'); 
  foreach (@rearr){
    if ($_ !~/\w+(Df|Dp|Ex|T|In|C|D)\d*/){
      print LOG "WARNING: $_ is NOT an object of Rearrangement\n";
    }
  }  
} 

##########################
# Process Sequence class #
##########################

sub process_sequence {

  print"\n\nChecking Sequence class for errors:\n";
  print LOG "\n\nChecking sequence class for errors:\n";
  print LOG "-----------------------------------\n";

  my $get_seqs_with_multiple_loci=<<EOF;
  Table-maker -p "$def_dir/get_seq_has_multiple_loci.def" quit 
EOF

  my (%Seq_loci, $seq, $locus);
  my $dir = "/wormsrv1/geneace";
  
  open (FH, "echo '$get_seqs_with_multiple_loci' | $tace $dir | ") || die "Couldn't access geneace\n";
  while (<FH>){
    chomp($_);
    if ($_ =~ /^\"/){
      $_ =~ s/\"//g;
      ($seq, $locus)=split(/\s+/, $_);
      $Seq_loci{$seq}=$locus;
    }
  }
  foreach (keys %Seq_loci){
    print LOG "ERROR: $_ has multiple loci attached.\n";
  }
}   


#############################################

sub check_genetics_coords_mapping {

  print "\nChecking discrepancies in genetics/coords mapping of each CDS/Transcript:\n\n";
  print LOG "\nChecking discrepancies in genetics/coords mapping:\n\n";
  print JAHLOG "\nChecking discrepancies in genetics/coords mapping:\n\n";
  system ("/wormsrv2/scripts/get_interpolated_gmap.pl -db /wormsrv1/geneace -diff");
  my $map_diff = "/wormsrv2/logs/mapping_diff.".$rundate;
  open(IN, $map_diff) || die $!;
  while(<IN>){
    print LOG $_;
    print JAHLOG $_;
  }
}



#############################

sub loci_as_other_name {

  my ($def, $dir, $db) = @_;
  my ($main, $other_name);

  open (FH, "echo '$def' | $tace $dir | ") || die "Couldn't access geneace\n";
  while (<FH>){
    chomp $_;
    if ($_ =~ /^\"/){
      $_ =~ s/\"//g;
      $_ =~ /(.+)\s+(.+)/;     
      $main = $1;
      $other_name = $2;
      $other_name =~ s/\\//g;
      $other_name = $db->fetch('Locus', $other_name); 
      
      if ($other_name){
	if ($exceptions{$main}){  # see global @exceptions and %exceptions
	  print LOG "INFO: $main has $other_name as Other_name...$other_name is still a separate Locus object (exception)\n";
        }
        else {
	  print LOG "WARNING: $main has $other_name as Other_name...$other_name is still a separate Locus object\n";
        }
	if ($ace && !$exceptions{$main}){
	  print ACE "\n-R Locus : \"$other_name\" \"$main\"\n";
	  print ACE "\nLocus : \"$main\"\n";
          print ACE "Other_name \"$other_name\"\n";
        }
      }  
    }
  }
}



##############################

sub create_log_files{

  # Create history logfile for script activity analysis
  $0 =~ m/\/*([^\/]+)$/; system ("touch /wormsrv2/logs/history/$1.`date +%y%m%d`");

  # create main log file using script name for
  my $script_name = $1;
  $script_name =~ s/\.pl//; # don't really need to keep perl extension in log name

  $log = "/wormsrv2/logs/$script_name.$rundate.$$";
  open (LOG, ">$log") or die "cant open $log";
  print LOG "$script_name\n";
  print LOG "started at ",`date`,"\n";
  print LOG "=============================================\n";
  print LOG "\n";

  $jah_log = "/wormsrv2/logs/$script_name.jahlog.$rundate.$$";
  open(JAHLOG, ">>$jah_log") || die "Can't open $jah_log\n";
  print JAHLOG "This mail is generated automatically for CGC on $rundate\n"; 
  print JAHLOG "If you have any queries please email ck1\@sanger.ac.uk or krb\@sanger.ac.uk\n\n";
  print JAHLOG "=========================================================================\n";
  
  # create separate log with errors for Erich
  $caltech_log = "/wormsrv2/logs/geneace_check.caltech_log.$rundate.$$";
  open(CALTECHLOG,">$caltech_log") || die "cant open $caltech_log";
  print CALTECHLOG "$0 started at ",`date`,"\n";
  print CALTECHLOG "This mail is generated automatically for Caltech\n";
  print CALTECHLOG "If you have any queries please email ck1\@sanger.ac.uk or krb\@sanger.ac.uk\n\n";
  print CALTECHLOG "================================================================================================\n";

}


##############################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

sub array_comp{

  my(@union, @isect, @diff, %union, %isect, %count, $e);
  my ($ary1_ref, $ary2_ref)=@_;
  @union=@isect=@diff=();
  %union=%isect=();
  %count=();
  foreach $e(@$ary1_ref, @$ary2_ref){
    $count{$e}++;
  }
  foreach $e (keys %count){
    push (@union, $e);
    if ($count{$e}==2){push @isect, $e;}
    else {push @diff, $e;}
  } 
  return \@diff, \@isect, \@union;  # all returned into one array  
}


__END__

=head2 NAME - geneace_check.pl  

=head3 <USAGE> 
 

=head2 Options: [h or help] [d or debug] [c or class] [a or ace]

            All options have single letter or wordy aliases 

B<-help:>     
            Displays usage of the script: surprised?

B<-debug:>   
            Debug mode, follow by user name, eg. -d ck1 or -debug ck1

            This allows checking results mailed only to the user, 
            otherwise results would email to the gangs of Sanger 
            Wormbase group and Caltech (Erich) 

B<-class:>   
            Allows checking only specified classes in Geneace. 
            All classes will be checked without this option.
            Choosing multiple classes is supported. 
            Class names are case insensitive.
            Current valid class names: 

               locus
               allele
               laboratory
               strain
               rearrangement
               sequence
               evidence
               mapping

            For example: -c allele -c locus OR -class sequence -class rearrangment

B<-databse:> 
            Allows specifying path to a specific database.
            Default database path is /wormsrv1/geneace without this option.

            For example: -db /wormsrv2/autoace or -database /wormsrv1/someone/Test_DB
 

B<-ace:>     
            Allows generating ace file for fixing erros spotted by 
            this checking script.
            Default location and filename of ace file:
            /wormsrv1/geneace/CHECKS/geneace_check.rundate.processid.ace
                      
            For example: -a or -ace


B<-verbose:>
            Toggles on extra output.   Useful when running on command line and not on crob job
            For the locus class it will display each locus name as it is processes it and show
            (next to the locus name) a dot for each error it had found in that locus


=head3 <RUN geneace_check.pl>

            Geneace check is now set to run on Sundays.
            ##Run Geneace check over weekend
            20 7 * * 0 /wormsrv2/scripts/geneace_check.pl -a &
