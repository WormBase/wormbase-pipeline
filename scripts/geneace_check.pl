#!/usr/local/bin/perl5.6.1 -w
# 
# geneace_check.pl
#
# by Keith Bradnam
#
# Script to run consistency checks on the geneace database
#
# Last updated by: $Author: ck1 $
# Last updated on: $Date: 2003-02-03 09:32:43 $


use strict;
use lib "/wormsrv2/scripts/"; 
use Wormbase;
use Ace;
use Getopt::Long;


###################################################
# variables and command-line options with aliases # 
###################################################

my ($help, $debug, $database, $class, @class, $ace);
my $maintainers = "All";
our $tace = &tace;   # tace executable path
our ($log, $erichlog, $jahlog);

my $rundate = `date +%y%m%d`; chomp $rundate;
my $acefile = "/wormsrv2/logs/geneace_check_ACE.$rundate.$$";

GetOptions ("h|help"        => \$help,
            "d|debug=s"     => \$debug,
	    "c|class=s"     => \@class,
	    "db|database=s" => \$database,
            "ace|a"         => \$ace,  
           );

################################################ 
# choose database to query: default is Geneace #
################################################

my $default_db;

if ($database eq ""){
  $default_db = "/wormsrv1/geneace"; 
  print "\nUsing Geneace as default database\n\n";
}

else {
  $default_db = $database;
  print "\nUsing $database as default database path\n\n";
}

# Display help if required
&usage("Help") if ($help);

# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = "$debug" . '\@sanger.ac.uk');
}

&create_log_files;


# open a connection to geneace
my $db = Ace->connect(-path  => $default_db,
		      -program =>$tace) || do { print LOG "Connection failure: ",Ace->error; die();};


############## 
# MAIN LOOPS #
##############

# track errors in each class
our $locus_errors = 0;
our $lab_errors = 0;
our $allele_errors = 0;
our $strain_errors = 0;
our $rearrangement_errors = 0;
our $sequence_errors = 0;


####################################################
# Process various classes looking for errors:      #
# choose class to check - multiple classes allowed #
####################################################

if(!@class){
  print "Checking all classes in Geneace.....\n\n";
  &process_locus_class;
  &process_laboratory_class;
  &process_allele_class;
  &process_strain_class;
  &process_rearrangement;
  &process_sequence;
 }
 
else{
  foreach $class (@class){
    $class = lc($class);  # makes command line option case-insensitive
    if ($class eq "") {
      &process_locus_class;
      &process_laboratory_class;
      &process_allele_class;
      &process_strain_class;
      &process_rearrangement;
      &process_sequence;
    }
    if ($class =~ /locus/)                 {&process_locus_class}
    if ($class =~ /(laboratory|lab)/)      {&process_laboratory_class}
    if ($class =~ /allele/)                {&process_allele_class}
    if ($class =~ /strain/)                {&process_strain_class}
    if ($class =~ /(rearrangement|rearr)/) {&process_rearrangement}
    if ($class =~ /(sequence|seq)/)        {&process_sequence}
  }  
}

#######################################
# Tidy up and mail relevant log files #
#######################################

$db->close;
close(LOG);
close(ERICHLOG);
close(JAHLOG);

# Always mail to $maintainers (which might be a single user under debug mode)
mail_maintainer($0,$maintainers,$log);


# Also mail to Erich unless in debug mode
my $interested ="krb\@sanger.ac.uk, emsch\@its.caltech.edu, ck1\@sanger.ac.uk";
mail_maintainer($0,"$interested",$erichlog) unless $debug; 


# Email to Jonathan for problematic loci
my $CGC = "ck1\@sanger.ac.uk, krb\@sanger.ac.uk"; 
mail_maintainer($0,$CGC,$jahlog) unless $debug;

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
  print LOG "\nChecking Locus class for errors:\n";
  print LOG "--------------------------------\n";

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
 
  &find_new_loci_in_current_DB($get_seg_with_pseudogene_locus, $db);
   
  #Look for loci that are other_names and still are obj of ?Locus -> candidate for merging
  my $locus_has_other_name=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/locus_has_other_name.def" quit
EOF
  
  loci_as_other_name($locus_has_other_name, $default_db, $db);

  my $locus_to_CDS=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/locus_to_CDS.def" quit
EOF

  #loci_point_to_same_CDS($locus_to_CDS, $default_db);

  # check CGC_approved loci is XREF to existing Gene_Class   
  # check locus in geneace that are connected to CDS but not CGC_approved
  
  
  my $cgc_loci_not_linked_to_geneclass=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/CGC_loci_not_linked_to_geneclass.def" quit
EOF
  my $get_all_gene_class=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/get_all_gene_class.def" quit
EOF
  my $locus_to_CDS_but_not_CGC_approved=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/locus_to_CDS_but_not_CGC_approved.def" quit
EOF
    
  locus_CGC($cgc_loci_not_linked_to_geneclass, $get_all_gene_class, $locus_to_CDS_but_not_CGC_approved, $default_db);

  my $cds_of_each_locus=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/cds_of_each_locus.def" quit
EOF
  my $seq_name_of_each_locus=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/seq_name_of_each_locus.def" quit
EOF
  
  cds_name_to_seq_name($cds_of_each_locus, $seq_name_of_each_locus, $default_db);
    
  print LOG "\nThere are $locus_errors errors in $size loci.\n";

}

sub cds_name_to_seq_name {
  if ($ace){open (CDS, ">>$acefile") || die "Can't write to file!"}	
  my ($def1, $def2, $db) = @_;
  my ($locus, $cds, %locus_cds, %locus_seq_name);
  open (FH1, "echo '$def1' | tace $db | ") || die "Couldn't access $db\n";
  open (FH2, "echo '$def2' | tace $db | ") || die "Couldn't access $db\n";
  
  while (<FH1>){
    chomp($_);
    if ($_ =~ /\"(.+)\"\s+\"(.+)\"/){push (@{$locus_cds{$1}}, $2)}
  }
  while (<FH2>){
    chomp($_);
    if ($_ =~ /\"(.+)\"\s+\"(.+)\"/){push (@{$locus_seq_name{$1}}, $2)}
  }
  
  my (%ary1, @ary1, %ary2, @ary2, $ea1, $ea2);
  foreach (keys %locus_cds){
    my $ary2 = \@{$locus_cds{$_}}; my $ary1 = \@{$locus_seq_name{$_}};
    foreach (@$ary2){$ary2{$_}++}
    foreach (@$ary1){$ary1{$_}++}
    foreach $ea1 (@$ary1){
      if (!$ary2{$ea1}){	
        print  LOG "ERROR: Locus $_ has incorrect Sequence_name $ea1\n";
        $locus_errors++;
	if ($ace){
	  print CDS "\n\nLocus : \"$_\"\n";
          print CDS "-D Sequence_name \"$ea1\"\n";
        }
      }
    }
    foreach $ea2 (@$ary2){
      if (!$ary1{$ea2}){
        print LOG "ERROR: Locus $_ is missing Sequence_name $ea2\n";
        $locus_errors++;
        if ($ace){
	  print CDS "\n\nLocus : \"$_\"\n";   
          print CDS "Sequence_name \"$ea2\"\n";
        }
      }	 	
    }
  }  
} 

############################
# Process Laboratory class #
############################

sub process_laboratory_class{

  print "\n\nChecking Laboratory class for errors:\n";
  print LOG "\n\nChecking Laboratory class for errors:\n";
  print LOG "-------------------------------------\n";
  #grab lab details
  my @labs = $db->fetch(-class => 'Laboratory',
		        -name  => '*');

  # test for Allele_designation and Representative tags
  foreach my $lab (@labs){
    if(!defined($lab->at('CGC.Allele_designation')) && $lab ne "CGC"){  
      print LOG "WARNING: $lab has no Allele_designation tag present\n";
      $lab_errors++;
    }    
    if(!defined($lab->at('CGC.Representative')) && $lab ne "CGC"){  
      print LOG "WARNING: $lab has no Representative tag present\n";
      $lab_errors++;
    }  
  }
  print LOG "\nThere are $lab_errors errors in Laboratory class\n";
}
 
########################
# Process Allele class #
########################

sub process_allele_class{
 
  print"\n\nChecking Allele class for errors:\n";
  print LOG "\n\nChecking Allele class for errors:\n";
  print LOG "---------------------------------\n";

  my @alleles = $db->fetch('Allele','*');
  my ($allele, %allele_gene, $gene, $seq_name, @seq1, @seq2, @seqs, $cdb);

  $cdb = Ace->connect(-path  => '/wormsrv2/current_DB/',
	              -program =>$tace) || do { print LOG "Connection failure: ",Ace->error; die();};
  
  @seqs=Table_maker();

  #print scalar @seqs, "\n";
  my %seqs;
    foreach (@seqs){
      $seqs{$_}++;
  }

  # test for Location tag and if an allele is connected to multiple loci

  foreach $allele (@alleles){
    if(!defined($allele->at('Location'))){  
      print LOG "ERROR: $allele has no Location tag present\n";
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
	    print LOG "WARNING: Sequence tag of Allele $allele points to $seq, which can now become @LOCI.\n";
	    $allele_errors++;   
	  }
	}
      }
    }  
  }

  $cdb->close;

  foreach (keys %allele_gene){
    if ((scalar @{$allele_gene{$_}}) > 1){
      print LOG "ERROR: $_ is connected to more than one Loci: @{$allele_gene{$_}}\n";
      $allele_errors++; 
    }
  }

  my $allele_has_flankSeq_and_no_seq=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/allele_has_flankSeq_and_no_seq.def" quit
EOF

  allele_has_flankSeq_and_no_seq($allele_has_flankSeq_and_no_seq, $default_db);

  my $allele_has_predicted_gene_and_no_seq=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/allele_has_predicted_gene_and_no_seq.def" quit
EOF

  allele_has_predicted_gene_and_no_seq($allele_has_predicted_gene_and_no_seq, $default_db);

  print LOG "\nThere are $allele_errors errors in Allele class\n";
}

sub allele_has_flankSeq_and_no_seq {
  
  my ($def, $dir, $db) = @_;
  open (FH, "echo '$def' | tace $dir | ") || die "Couldn't access geneace\n";
  while (<FH>){
    chomp($_);
    if ($_ =~ /^\"/){
      $_ =~ s/\"//g;
      print LOG "WARNING: Allele $_ has flanking sequences but is NOT connected to sequence\n"; 
      $allele_errors++; 
    }
  }
}

sub allele_has_predicted_gene_and_no_seq {
  
  my ($def, $dir, $db) = @_;
  open (FH, "echo '$def' | tace $dir | ") || die "Couldn't access geneace\n";
  while (<FH>){
    chomp($_);
    if ($_ =~ /^\"/){
      $_ =~ s/\"//g;
       print LOG "WARNING: Allele $_ has predicted gene but is NOT connected to sequence\n"; 
      $allele_errors++; 
    }
  }
}

########################
# Process Strain class #  
########################

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
  foreach $strain (@strains){
    if (!$strain->Location){
      print LOG "WARNING: Strain $strain has no location tag\n";
      $strain_errors++; 
    }
    else { 
      my $cgc=$strain->Location;
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
	      if ($cgc eq "CGC"){	
		print LOG "WARNING: CGC Strain $strain has sequence_name $_ in Genotype, which can now become @loci.\n";
		$strain_errors++; 
	      }
	      else {
		print LOG "WARNING: Non_CGC Strain $strain has sequence_name $_ in Genotype, which can now become @loci.\n";
		$strain_errors++; 
	      }  
	    }
	  }
	}
      }
    }
  }
  print LOG "\nThere are $strain_errors errors in Strain class.\n";
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
      $rearrangement_errors++;
      print LOG "WARNING: $_ is NOT an object of Rearrangement\n";
    }
  }  
  print LOG "\n\nThere are $rearrangement_errors errors in Rearrangement class.\n";
} 

##########################
# Process Sequence class #
##########################

sub process_sequence {

  print"\n\nChecking Sequence class for errors:\n";
  print LOG "\n\nChecking sequence class for errors:\n";
  print LOG "-----------------------------------\n";

  my $get_seqs_with_multiple_loci=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/get_seq_has_multiple_loci.def" quit 
EOF

  my (%Seq_loci, $seq, $locus);
  my $dir = "/wormsrv1/geneace";
  
  open (FH, "echo '$get_seqs_with_multiple_loci' | tace $dir | ") || die "Couldn't access geneace\n";
  while (<FH>){
    chomp($_);
    if ($_ =~ /^\"/){
      $_ =~ s/\"//g;
      ($seq, $locus)=split(/\s+/, $_);
      $Seq_loci{$seq}=$locus;
    }
  }
  foreach (keys %Seq_loci){
    $sequence_errors++;
    print LOG "ERROR: $_ has multiple loci attached.\n";
  }
  print LOG "\n\nThere are $sequence_errors errors in Sequence class\n";
}   

##############################

sub find_new_loci_in_current_DB{
  my ($def, $db) = @_;
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
      print LOG "WARNING: $gene has no Pseudogene tag, but its corresponding seq does.\n";
    }
  }
}

#############################

sub loci_as_other_name {

  my ($def, $dir, $db) = @_;
  my ($main, $other_name);

  open (FH, "echo '$def' | tace $dir | ") || die "Couldn't access geneace\n";
  while (<FH>){
    chomp($_);
    if ($_ =~ /^\"/){
      $_ =~ s/\"//g;
      $_ =~ /(\w+-\d+|\w+-\d+\.[\d\w]+|.+)\s+(.+)/;
     # print "$1 -> $2\n";
      $main = $1;
      $other_name = $2;
      $other_name =~ s/\\//g;
      $other_name = $db->fetch('Locus', "$other_name"); 
      if ($other_name){
	$locus_errors++;
	print LOG "WARNING: $main has $other_name as Other_name...$other_name is still a separate Locus object\n";
      }
    }
  }
}

#############################

sub loci_point_to_same_CDS {
 
   my ($def, $dir) = @_;
   my (%CDS_loci);

   open (FH, "echo '$def' | tace $dir | ") || die "Couldn't access geneace\n";
   while (<FH>){
     chomp($_);
     if ($_ =~ /^\"/){
       $_ =~ s/\"//g;
       $_ =~  /(\w+-\d+|\w+-\d+\.[\d\w]+|.+)\s+(.+)/;
       push(@{$CDS_loci{$2}}, $1);
     }
   }
   foreach (keys %CDS_loci){
     if (scalar @{$CDS_loci{$_}} > 1){
       $locus_errors++;
       print LOG "ERROR: $_ is connected to @{$CDS_loci{$_}}: case of merging?\n";
       print JAHLOG "ERROR: $_ is connected to @{$CDS_loci{$_}}: case of merging?\n";
     }
   }

} 
  
##############################

sub locus_CGC {

  if ($ace){open (LC, ">>$acefile") || die "Can't write to file!"}
  my ($def1, $def2, $def3, $dir) = @_;
  my (@gc, %gc, $gc);

  open (FH, "echo '$def2' | tace $dir | ") || die "Couldn't access geneace\n"; 
  while (<FH>){
    chomp($_);
    if ($_ =~ /^\"/){
      $_ =~ s/\"//g;
      push(@gc, $_);
    }
  }
  foreach(@gc){$gc{$_}++}
  
  
  open (FH, "echo '$def1' | tace $dir | ") || die "Couldn't access geneace\n";
  while (<FH>){
    chomp($_);
    if ($_ =~ /^\"/){
      $_ =~ s/\"//g;
      $gc = $_;
      $gc =~ s/-\d+//;
      if ($gc{$gc}){ 
	print LOG "ERROR: $_ is CGC_approved but not XREF to existing $gc Gene_Class\n"; 
	$locus_errors++; 
        if ($ace){
	  print LC "\n\nGene_Class : \"$gc\"\n";
          print LC "Loci \"$_\"\n";
        }
      }
    }
  }
  open (FH, "echo '$def3' | tace $dir | ") || die "Couldn't access geneace\n";
  while (<FH>){
    chomp($_);
    if ($_ =~ /^\"/){
      $_ =~ s/\"//g;
      if ($_ !~ /^[A-Z]/){
	print LOG "WARNING: $_ is linked to coding sequence but not CGC_approved\n"; 
	print JAHLOG "WARNING: $_ is linked to coding sequence but not CGC_approved\n"; 	
	$locus_errors++; 
      }
    }
  }
}

##############################

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


# This test no longer needed as we don't have New_name tag now  
  # test for New_name tag and other information apart from Gene and Species
#  if(defined($locus->at('Name.New_name'))){
#    if(defined($locus->at('Species')) && defined($locus->at('Type.Gene'))){
#      # check for any other info in object
#      my @tags = $locus->tags();
#      if(scalar(@tags)>3){
#	$warnings .= "ERROR 12: $locus has New_name tag + extra info. Transfer into new gene?\n";
#	$locus_errors++;
#      }
#    }
#    else{
#      $warnings .= "ERROR 13: $locus has no species and/or Gene tag present\n";
#      $locus_errors++;
#    }
#  }

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

  $jahlog = "/wormsrv2/logs/$script_name.jahlog.$rundate.$$";
  open(JAHLOG, ">$jahlog") || die "Can't open $jahlog\n";
  print JAHLOG "Data created on $rundate\n";
  print JAHLOG "=============================================\n";
  
  # create separate log with errors for Erich
  $erichlog = "/wormsrv2/logs/geneace_check.erichlog.$rundate.$$";
  open(ERICHLOG,">$erichlog") || die "cant open $erichlog";
  print ERICHLOG "$0 started at ",`date`,"\n";
  print ERICHLOG "=============================================\n";
  print ERICHLOG "This mail is generated automatically.\n";
  print ERICHLOG "If you have any queries please email ar2\@sanger.ac.uk or krb\@sanger.ac.uk\n\n";
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

##############################

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

