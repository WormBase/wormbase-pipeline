#!/usr/local/bin/perl5.6.1 -w
# 
# geneace_check.pl
#
# by Anthony Rogers based on a script by Keith Bradnam
#
# Script to run consistency checks on the geneace database
#
# Last updated by: $Author: ck1 $
# Last updated on: $Date: 2002-11-29 17:05:01 $

use Ace;
use lib "/wormsrv2/scripts/"; 
use Wormbase;
use strict;

our $log;
our $erichlog;
&create_log_files;

# track errors in each class
our $locus_errors = 0;
our $lab_errors = 0;
our $allele_errors = 0;
our $strain_update = 0;

# open a connection to geneace

our $tace = &tace;   # tace executable path

my $db = Ace->connect(-path  => '/wormsrv1/geneace/',
		      -program =>$tace) || do { print LOG "Connection failure: ",Ace->error; die();};


# Process various classes looking for errors
&process_locus_class;
&process_laboratory_class;
&process_allele_class;
&process_strain_class;
&rearrangement;


#######################################
# Tidy up and mail relevant log files #
#######################################

$db->close;
close(LOG);
close(ERICHLOG);

my $maintainer = "All";
#&mail_maintainer($0,$maintainer,$log);

my $interested ="krb\@sanger.ac.uk, emsch\@its.caltech.edu, ck1\@sanger.ac.uk";
#&mail_maintainer($0,"$interested",$erichlog);
exit(0);


#######################
# Process Locus class #
#######################

sub process_locus_class{

  my @loci = $db->fetch(-class => 'Locus',
		      -name  => '*');
  
  my $size =scalar(@loci);
  
  # Loop through loci checking for various potential errors in the Locus object
  print "\nChecking loci for errors:\n";
  print LOG "\nChecking Locus class for errors\n";
  foreach my $locus (@loci){
    #print "$locus\n";
    my $warnings;
    my $erich_warnings;
    ($warnings, $erich_warnings) = &test_locus_for_errors($locus);
    print LOG "$warnings" if(defined($warnings));
    #Erich Schwarz wants some of these - emsch@its.caltech.edu
    print ERICHLOG "$erich_warnings" if(defined($erich_warnings));
  }

  # Look for loci in current_DB not in geneace
  # Look for sequence in current_DB that is a pseudogene and has locus connection
  print "\nLooking for new loci in /wormsrv2/current_DB:\n\n";

  my $get_seg_with_pseudogene_locus=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/get_all_seq_with_pseudogene_and_locus.def" quit
EOF

  &find_new_loci_in_current_DB($db, $get_seg_with_pseudogene_locus);
  print LOG "\nThere were $locus_errors errors in $size loci.\n";
}

############################
# Process Laboratory class #
############################

sub process_laboratory_class{

  print "\n\nChecking Laboratory class for errors:\n";
  print LOG "\n\nChecking Laboratory class for errors:\n";
  #grab lab details
  my @labs = $db->fetch(-class => 'Laboratory',
		        -name  => '*');

  # test for Allele_designation and Representative tags
  foreach my $lab (@labs){
    if(!defined($lab->at('CGC.Allele_designation')) && $lab ne "CGC"){  
      print LOG "$lab has no Allele_designation tag present\n";
      #print  "$lab has no Allele_designation tag present\n"; 
      $lab_errors++;
    }    
    if(!defined($lab->at('CGC.Representative')) && $lab ne "CGC"){  
      print LOG "$lab has no Representative tag present\n";
      #print  "$lab has no Representative tag present\n";
      $lab_errors++;
    }  
  }
  print LOG "\nThere were $lab_errors errors in Laboratory class\n";
  #print "\nThere were $lab_errors errors in Laboratory class\n";
}
 
########################
# Process Allele class #
########################

sub process_allele_class{
 
  print"\n\nChecking Allele class for errors:\n";
  print LOG "\n\nChecking Allele class for errors:\n";

  my @alleles = $db->fetch('Allele','*');
  print scalar @alleles, "\n";
  my ($allele, %allele_gene, $gene, $seq_name, @seq1, @seq2, @seqs, $cdb);

  $cdb = Ace->connect(-path  => '/wormsrv2/current_DB/',
	              -program =>$tace) || do { print LOG "Connection failure: ",Ace->error; die();};
  
  @seqs=Table_maker();

  print scalar @seqs, "\n";
  my %seqs;
    foreach (@seqs){
      $seqs{$_}++;
  }

  # test for Location tag and if an allele is connected to multiple loci

  foreach $allele (@alleles){
    if(!defined($allele->at('Location'))){  
      print LOG "$allele has no Location tag present\n";
      $allele_errors++;
    }
   
    # checking if sequence name in Allele has now a locus name 

    if($allele -> Gene){
      my @loci=$allele->Gene(1);
      push(@{$allele_gene{$allele}}, @loci); 
      foreach $gene (@loci){
	if($seqs{$gene}){
	  my $seq = $cdb->fetch('Sequence', $gene);
	  if ($seq->Locus_genomic_seq){
	    my @LOCI=$seq->Locus_genomic_seq(1);
	    print LOG "Sequence tag of Allele $allele points to $seq, which can now become @LOCI.\n";
	    $allele_errors++;   
	  }
	}
      }
    }  
  }

  foreach (keys %allele_gene){
    if ((scalar @{$allele_gene{$_}}) > 1){
      print LOG "$_ is connected to more than one Loci: @{$allele_gene{$_}}\n";
      $allele_errors++; 
    }
  }
  print LOG "\nThere were $allele_errors errors in Allele class\n";
}

########################
# Process Strain class #  
########################

sub process_strain_class {

  # Check if sequence name of a strain genotype is now connected to a locus

  my (@strains, @seqs, $genotype, @genes, $seq, $strain, $extract);

  print"\n\nChecking Strain class for errors:\n";
  print LOG "\n\nChecking Strain class for errors:\n";
  @strains = $db->fetch('Strain','*');
  @seqs = $db->fetch('Sequence','*');

  my %seqs;
  foreach (@seqs){
    $seqs{$_}++;
  }
  foreach $strain (@strains){
    if ($strain->Genotype){
      $genotype = $strain->Genotype(1);
      $extract = $genotype;
      $extract =~ s/\(|\)|\/|\+|;|\?|\{|\}|,|\=|\.$/ /g;
      $extract =~ s/ I | II | III | IV | V | X / /g;
      $extract =~ s/ I | II | III | IV | V | X | f / /g; # further tidying up of chromosomes
      $extract =~ s/^\s|\w{3}-\s| f | A //g;
      $extract =~ s/\s{1,}/ /g;
      @genes=split(/ /,$extract);
      foreach (@genes){
	if($seqs{$_}){
	  my $seq = $db->fetch('Sequence', $_);
	  if ($seq->Locus_genomic_seq){
	    my @loci=$seq->Locus_genomic_seq(1);
	    print LOG "Strain $strain has sequence_name $_ in Genotype, which can now become @loci.\n";
	    $strain_update++; 
	  }  
	}
      }
    }
  }
  print LOG "\nThere are $strain_update genotypes to be updated in Strain class.\n";
  #print "\nThere are $strain_update genotypes to be updated in Strain class.\n";
} 

###############################
# Process Rearrangement class #
###############################

sub rearrangement {
 
  print"\n\nChecking Rearrangement class for errors:\n";
  print LOG "\n\nChecking Rearrangement class for errors:\n";
  # checks presence of non-rearrangement object 
  # as objects of Rearrangement class

  my (@rearr, $count);
 
  @rearr = $db -> fetch('Rearrangement','*'); 
  foreach (@rearr){
    if ($_ !~/\w+(Df|Dp|Ex|T|In|C|D)\d*/){
      $count++;
      print LOG "$_ is NOT an object of Rearrangement\n";
    }
  }  
  print LOG "\n\nThere are $count error(s) in Rearrangement class.\n";
} 

################################################

sub find_new_loci_in_current_DB{
  my ($db, $def) = @_;
  my $warnings;
  my @genes=();
  my $dir="/wormsrv2/current_DB";
  my $locus_errors=0;

  # open a database connection to current_DB and grab all loci names (excluding polymorphisms)
  my $new_db = Ace->connect(-path  => '/wormsrv2/current_DB',
		    -program =>$tace) || do { print LOG "Connection failure: ",Ace->error; die();};
  my @current_DB_loci = $db->fetch(-query=>'Find Locus;!Polymorphism');

  open (FH, "echo '$def' | tace $dir | ") || die "Couldn't access geneace\n";
  while (<FH>){
    chomp($_);
    if ($_ =~ /^\"/){
      $_ =~ s/\"|//g;
      $_ =~ s/\s+/ /g;
      my @items=split(/ /, $_);
      push (@genes, $items[2]);
    }
  }
    
  $new_db->close;

  #cross reference in geneace
  foreach my $loci(@current_DB_loci){
    my $new_loci = $db->fetch(-class=>'Locus',-name=>"$loci");
    unless(defined($new_loci)){
      $warnings .= "ERROR 19: $new_loci in current_DB is not in /wormsrv1/geneace\n";
      $locus_errors++;
    }
  }
  print LOG "\n$warnings\n" if $warnings;

  # check geneace locus w/o pseudogene tag
  foreach (@genes){
    my $gene = $db->fetch('Locus',$_);
    if (!$gene->at('Type.Gene.Pseudogene')){
      $locus_errors++;
      print LOG "$gene has no Pseudogene tag, but its corresponding seq does.\n";
    }
  }
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
      $locus_errors++;
    }
  }

  # test for more than one gene class
  if(defined($locus->at('Name.Gene_class'))){
    my @gene_classes = $locus->at('Name.Gene_Class');
    if(scalar(@gene_classes) > 1){
      $warnings .= "ERROR 2: $locus has more than one Gene_class\n";
      $locus_errors++;
    }    
  }
  
  # test for no Type tag
  if(!defined($locus->at('Type'))){  
    $warnings .= "ERROR 3: $locus has no Type tag present\n";
    $locus_errors++;
  }

  # test for Gene AND !Species 
  if(defined($locus->at('Type.Gene')) && !defined($locus->at('Species'))){
    $warnings .= "ERROR 4: $locus has a 'Gene' tag but not a 'Species' tag\n";;
    $locus_errors++;
  }

  # test for !Gene AND Gene_class 
  if(!defined($locus->at('Type.Gene')) && defined($locus->at('Name.Gene_class'))){
    $warnings .= "ERROR 5: $locus has a 'Gene_class' tag but not a 'Gene' tag\n";
    $locus_errors++;
  }
  
  # test for no Gene tag AND Genomic_sequence tag
  if(!defined($locus->at('Type.Gene')) && defined($locus->at('Molecular_information.Genomic_sequence'))){
    $warnings .= "ERROR 6: $locus has 'Genomic_sequence' tag but no 'Gene' tag\n";
    $locus_errors++;
  }

  # test for Genomic_sequence tag but no value   
  if(defined($locus->at('Molecular_information.Genomic_sequence')) && !defined($locus->Genomic_sequence)){
    $warnings .= "ERROR 7: $locus has 'Genomic_sequence' tag but no associated value\n";
    $locus_errors++;
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
	  $locus_errors++;
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
	$locus_errors++;
      }
    }
  }
  
  # test for Polymorphism AND Gene tags both present
  if(defined($locus->at('Type.Gene')) && defined($locus->at('Type.Polymorphism'))){  
    $warnings .= "ERROR 9: $locus has a 'Polymorphism' tag AND 'Gene' tag\n";
    $locus_errors++;
  }

  # test for Enzyme tag (doesn't fit in with new Gene model)
  if(defined($locus->at('Molecular_information.Enzyme'))){  
    $warnings .= "ERROR 10: $locus has an 'Enzyme' tag\n";
    $locus_errors++;
  }

  # test for Canonical_gene tag present
  if(defined($locus->at('Molecular_information.Canonical_gene'))){  
    $warnings .= "ERROR 11: $locus has a 'Canonical_gene' tag\n";
    $locus_errors++;
  }
  
  # test for New_name tag and other information apart from Gene and Species
  if(defined($locus->at('Name.New_name'))){
    if(defined($locus->at('Species')) && defined($locus->at('Type.Gene'))){
      # check for any other info in object
      my @tags = $locus->tags();
      if(scalar(@tags)>3){
	$warnings .= "ERROR 12: $locus has New_name tag + extra info. Transfer into new gene?\n";
	$locus_errors++;
      }
    }
    else{
      $warnings .= "ERROR 13: $locus has no species and/or Gene tag present\n";
      $locus_errors++;
    }
  }
=start  
  # Test for Polymorphisms with no P in their title
  if(defined($locus->at('Type.Polymorphism'))){
    if($locus !~ /P/){
      $warnings .= "ERROR 14: $locus has no 'P' in title\n";
      $locus_errors++;
    }
  }
=end
=cut
  # Look for Gene_class tag in non-gene objects 
  if(!defined($locus->at('Type.Gene'))){
    if(defined($locus->at('Name.Gene_class'))){
      $warnings .= "ERROR 15: $locus has Gene_class tag but it is not a gene!\n";
      $locus_errors++;
    }
  }
  
  # test for Other_name tag but no value   
  if(defined($locus->at('Name.Other_name')) && !defined($locus->Other_name)){
    $warnings .= "ERROR 16: $locus has 'Other_name' tag but no associated value\n";
    $locus_errors++;
  }
=start
  # test for Other_name value which is also a Locus name in its own right
  # N.B. This is not always a bad thing
  if(defined($locus->at('Name.Other_name'))){
    my @other_names = $locus->Other_name;
    foreach my $other_name (@other_names){
      my ($newloci) = $db->fetch(-class=>'Locus',-name=>"$other_name");
      if(defined($newloci)){
	unless($newloci->at('Name.New_name')){
	  $warnings .= "ERROR 17: $locus has Other_name: $other_name which is also a separate locus object\n";
	  $locus_errors++;
	}   
      }
    }
  }
=end
=cut
  # Remind of outstanding CGC_unresolved tags
  if(defined($locus->CGC_unresolved)){
    my ($unresolved_details) = $locus->at('Type.Gene.CGC_unresolved');
    $warnings .= "ERROR 18: $locus has CGC_unresolved tag: \"$unresolved_details\"\n";
    $locus_errors++;
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

sub Table_maker {

  my $get_predicted_genes=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/get_all_predicted_gene_names.def" quit 
EOF
  my $get_genome_seqs=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/get_all_genome_sequence_names.def" quit
EOF

  my $dir="/wormsrv2/current_DB";
  
  my @names=();
  open (FH1, "echo '$get_predicted_genes' | tace $dir | ") || die "Couldn't access current_DB\n";
  open (FH2, "echo '$get_genome_seqs' | tace $dir | ") || die "Couldn't access current_DB\n";
  
    while (<FH1>){
      chomp($_);
      if ($_ =~ /^\"/){ 
	$_ =~ s/"//g;
        #print $_, "\n";
        push (@names, $_);
      }
    }
    while (<FH2>){
      chomp($_);
      if ($_ =~ /^\"/){ 
	$_ =~ s/"//g;
        #print $_, "\n";
        push (@names, $_);
      }
    }
  return @names;
}

