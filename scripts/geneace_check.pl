#!/usr/local/bin/perl5.6.1 -w
# 
# geneace_check.pl
#
# by Keith Bradnam
#
# Script to run consistency checks on the geneace database
#
# Last updated by: $Author: krb $
# Last updated on: $Date: 2003-06-19 09:07:55 $

use strict;
use lib "/wormsrv2/scripts/"; 
use Wormbase;
use Ace;
use Getopt::Long;


###################################################
# variables and command-line options with aliases # 
###################################################

my ($help, $debug, $database, $class, @class, $ace, $verbose);
my $maintainers = "All";

# hashes for checking Person and Author merging?
my (%L_name_F_WBP, %L_name_F_M_WBP);


our $tace = &tace;   # tace executable path
our ($log, $erichlog, $jahlog, $JAHmsg, $Emsg, $caltech, @CGC, $cgc, $reverse_log, $map_diff);

my $rundate = `date +%y%m%d`; chomp $rundate;
my $acefile = "/wormsrv2/logs/geneace_check_ACE.$rundate.$$";
my $curr_db = "/nfs/disk100/wormpub/DATABASES/current_DB"; 

GetOptions ("h|help"        => \$help,
            "d|debug=s"     => \$debug,
	    "c|class=s"     => \@class,
	    "db|database=s" => \$database,
            "a|ace"         => \$ace, 
	    "v|verbose"        => \$verbose
           );



# Display help if required
&usage("Help") if ($help);


##########################################################
# choose database to query: default is /wormsrv1/geneace #
##########################################################

my $default_db;

if ($database eq ""){
  $default_db = "/wormsrv1/geneace"; 
}
else {
  $default_db = $database;
}
print "\nUsing $default_db as default database.\n\n";


# Use debug mode?
if($debug){
  print "DEBUG = \"$debug\"\n\n";
  ($maintainers = "$debug" . '\@sanger.ac.uk');
}

# Open output ace file if specified
if ($ace){open (ACE, ">>$acefile") || die "Can't write to file!\n"}

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
  print "Checking all classes in $default_db.....\n\n";
  &process_locus_class;
  &process_laboratory_class;
  &process_allele_class;
  &process_strain_class;
  &process_rearrangement;
  &process_sequence;
  &check_genetics_coords_mapping;
  #&chech_reverse_physicals;
  &check_evidence;
 }
 
else{
  foreach $class (@class){
    $class = lc($class);  # makes command line option case-insensitive
    if ($class =~ /locus/)                 {&process_locus_class}
    if ($class =~ /(laboratory|lab)/)      {&process_laboratory_class}
    if ($class =~ /allele/)                {&process_allele_class}
    if ($class =~ /strain/)                {&process_strain_class}
    if ($class =~ /(rearrangement|rearr)/) {&process_rearrangement}
    if ($class =~ /(sequence|seq)/)        {&process_sequence}
    if ($class =~ /(mapping|map)/)         {&check_genetics_coords_mapping}
    if ($class =~ /(reverse|rev)/)         {&chech_reverse_physicals}
    if ($class =~ /(evidence|evi)/)        {&check_evidence}
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

open(MAIL1, "$erichlog") || die "Can't read in file $erichlog";

my @caltech=<MAIL1>;
$caltech=join('', @caltech);
if ($caltech ne $Emsg){      
  mail_maintainer($0,"$interested",$erichlog) unless $debug; 
}

# Email to Jonathan for problematic loci 
my $CGC = "ck1\@sanger.ac.uk, krb\@sanger.ac.uk"; 

open(MAIL2, "$jahlog") || die "Can't read in file $erichlog";

@CGC=<MAIL2>;
$cgc=join('', @CGC);

if ($cgc ne $JAHmsg){   
  mail_maintainer($0,$CGC,$jahlog) unless $debug;
}

chdir "/wormsrv2/logs";

my @files = qw ($acefile  $log $jahlog $JAHmsg $erichlog $Emsg $reverse_log $map_diff);
foreach(@files){
  system("$_") if $_;
}

exit(0);


##############################################################################################
#
#
#                                    geneace_check subroutines
#
#
#
##############################################################################################



#######################################
# Check misuse of Evidence in 8 classes 
# Convert Author to Person
#######################################

sub check_evidence {


my $WBPerson_F_M_L_names=<<EOF;
Table-maker -p "/wormsrv1/geneace/wquery/WBperson_first_middle_last_names.def" quit
EOF

  my @WBPerson = &process_WBPerson_names($WBPerson_F_M_L_names, $curr_db);
 
  print LOG "\nChecking misuse of Evidence and converting Author to Person / Non-Person to Author:\n";
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

# look for person/author names that needs to be converted

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
	$evid_errors++;
	print LOG "\nERROR: $class $obj has non-Person $name under main tag $tag\n";

	@counters = get_WBPerson_ID($name, $class, $obj, $tag, $ori, $b4_evi, "PtoA");  
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
      if ($paper !~ /\[.+]/){
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
  print LOG "\n\nThere are $evid_errors Evidence errors in 7 classes checked\n";
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
    open (FH, "echo '$def' | tace $db | ") || die "Couldn't access current_DB\n";
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
    undef($locus);
  }

  # Look for loci in current_DB not in geneace
  # Look for sequence in current_DB that is a pseudogene and has locus connection
  print "\nLooking for new loci in /nfs/disk100/wormpub/DATABASES/current_DB:\n\n";

  my $get_seg_with_pseudogene_locus=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/get_all_seq_with_pseudogene_and_locus.def" quit
EOF
 
  &find_new_loci_in_current_DB($get_seg_with_pseudogene_locus, $db);
   
  #Look for loci that are other_names and still are obj of ?Locus -> candidate for merging
 
  my $locus_has_other_name=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/locus_has_other_name.def" quit
EOF
  
  &loci_as_other_name($locus_has_other_name, $default_db, $db);

  my $locus_to_CDS=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/locus_to_CDS.def" quit
EOF

  &loci_point_to_same_CDS($locus_to_CDS, $default_db);

  my $cgc_approved_and_non_cgc_name=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/cgc_approved_and_non_cgc_name.def" quit
EOF
  my $cgc_approved_has_no_cgc_name=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/cgc_approved_has_no_cgc_name.def" quit
EOF
  
  &gene_name_class($cgc_approved_and_non_cgc_name, $cgc_approved_has_no_cgc_name, $default_db); 

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
    
  &locus_CGC($cgc_loci_not_linked_to_geneclass, $get_all_gene_class, $locus_to_CDS_but_not_CGC_approved, $default_db);

  my $cds_of_each_locus=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/cds_of_each_locus.def" quit
EOF
  my $seq_name_of_each_locus=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/seq_name_of_each_locus.def" quit
EOF
  
  &cds_name_to_seq_name($cds_of_each_locus, $seq_name_of_each_locus, $default_db);

 my $no_remark_in_geneclass_for_merged_loci=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/no_remark_in_geneclass_for_merged_loci.def" quit
EOF

  &add_remark_for_merged_loci_in_geneclass($no_remark_in_geneclass_for_merged_loci, $default_db);

    
  print LOG "\nThere are $locus_errors errors in $size loci.\n";

}

##############################################################

sub add_remark_for_merged_loci_in_geneclass {
  my ($def, $db)=@_;
  my ($locus, $other, $gc);
  open (FH, "echo '$def' | tace $db | ") || die "Couldn't access $db\n";
  while (<FH>){
    chomp $_;
    if ($_ =~ /^\"(.+)\"\s+\"(.+)\"\s+\"(.+)\"/){
      $locus_errors++;
      $locus = $1;
      $other = $2;
      $gc = $3;
      $locus =~ s/\\//g;  # get rid of \ in locus like AAH\/1 from table maker
      print LOG "$locus became an other_name of $other and no remark is added in Gene_class $gc\n";
      if ($ace){
	print ACE "\n\nGene_Class : \"$gc\"\n";
	print ACE "Remark \"$locus is also an unofficial other name of $other\" CGC_data_submission\n";
      }
    }
  }
}



##############################################################


sub gene_name_class {
  my ($def1, $def2, $db) = @_;
  my $locus;
  open (FH1, "echo '$def1' | tace $db | ") || die "Couldn't access $db\n"; 
  open (FH2, "echo '$def2' | tace $db | ") || die "Couldn't access $db\n";
  
  while (<FH1>){
    chomp $_;
    if ($_ =~ /^\"(.+)\"/){
      $locus = $1;
      $locus_errors++;
      print LOG "ERROR: $locus is CGC_approved but still has NON_CGC_name tag\n";
      if ($ace){
	print ACE "\n\nLocus : \"$locus\"\n";
	print ACE "-D Non_CGC_name\n";
      }	
    }
  }
  
  while (<FH2>){
    chomp $_;
    if ($_ =~ /^\"(.+)\"/){
      $locus = $1;
      print LOG "ERROR: $locus is CGC_approved but has no CGC_name tag\n";
      $locus_errors++;
      if ($ace){
	print ACE "\n\nLocus : \"$locus\"\n";
	print ACE "CGC_name \"$locus\"\n";
      }	
    }
  }    
}

##############################################################

sub cds_name_to_seq_name {

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
	  print ACE "\n\nLocus : \"$_\"\n";
          print ACE "-D Sequence_name \"$ea1\"\n";
        }
      }
    }
    foreach $ea2 (@$ary2){
      if (!$ary1{$ea2}){
        print LOG "ERROR: Locus $_ is missing Sequence_name $ea2\n";
        $locus_errors++;
        if ($ace){
	  print ACE "\n\nLocus : \"$_\"\n";   
          print ACE "Sequence_name \"$ea2\"\n";
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
    if($verbose && !defined($lab->at('CGC.Allele_designation')) && $lab =~ /\w{3,3}/){  
      print LOG "INFO: $lab has no Allele_designation tag present (exception)\n";
    }    
    
    if(!defined($lab->at('CGC.Representative')) && $lab ne "CGC"){  
      if ($lab ne "XA"){
	print LOG "WARNING: $lab has no Representative tag present\n";
	$lab_errors++;
      }
    }  
    undef($lab);
  }
  print  LOG "\nThere are $lab_errors errors in Laboratory class\n";
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

  $cdb = Ace->connect(-path  => '/nfs/disk100/wormpub/DATABASES/current_DB/',
	              -program =>$tace) || do { print LOG "Connection failure: ",Ace->error; die();};
  
  @seqs=Table_maker();

  #print scalar @seqs, "\n";
  my %seqs;
  foreach (@seqs){$seqs{$_}++}

  # check if an allele has no Location tag 
  #          an allele is connected to multiple loci
  #          the sequence tag of an allele has a locus name

  my $allele_designation_to_LAB=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/allele_designation_to_LAB.def" quit
EOF
  
  my %location;
  if ($ace){%location=allele_location($allele_designation_to_LAB, $default_db)};
 
  foreach $allele (@alleles){
    if(!defined($allele->at('Location'))){
      if ($allele =~ /^[A-Z].+/){
        $allele_errors++;
	print LOG "ERROR: $allele has no Location tag present (no info available)\n";
      }
      else {
	print LOG "ERROR: $allele has no Location tag present\n";
	$allele_errors++;
	if ($ace){
	  my $desig = $allele;        
	  $desig =~ s/\d+//;
	  if (exists $location{$desig}){
	    print  ACE "\n\nAllele : \"$allele\"\n";
	    print  ACE "Location \"$location{$desig}\"\n";
	  }
	}	
      }
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
    undef($allele);
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

  #&allele_has_flankSeq_and_no_seq($allele_has_flankSeq_and_no_seq, $default_db);

  my $allele_has_predicted_gene_and_no_seq=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/allele_has_predicted_gene_and_no_seq.def" quit
EOF

  #&allele_has_predicted_gene_and_no_seq($allele_has_predicted_gene_and_no_seq, $default_db);

  my $allele_methods=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/allele_methods.def" quit
EOF

  check_missing_allele_method($allele_methods, $default_db);

  print LOG "\nThere are $allele_errors errors in Allele class\n";
}


#################################################################

sub allele_location {
  my ($def, $dir)=@_;
  my %location_desig;
  open (FH, "echo '$def' | tace $dir | ") || die "Couldn't access geneace\n";
  while (<FH>){
    chomp($_);
    if ($_ =~ /^\"(.+)\"\s+\"(.+)\"/){
      $location_desig{$2} = $1  # $2 is allele_designation $1 is LAB	
    }
  }
  return %location_desig;
}
  
#################################################################

sub allele_has_flankSeq_and_no_seq {
      
  my ($def, $dir) = @_;
  open (FH, "echo '$def' | tace $dir | ") || die "Couldn't access geneace\n";
  while (<FH>){
    chomp($_);
    if ($_ =~ /^\"/){
      $_ =~ s/\"//g;
      print LOG "WARNING: Allele $_ has flanking sequences but is NOT connected to parent sequence\n";
      $allele_errors++; 
    }
  }
}


#################################################################

sub allele_has_predicted_gene_and_no_seq {
          
  my ($def, $dir) = @_;
  my ($allele, $seq, $parent, $cds);
      
  open (FH, "echo '$def' | tace $dir | ") || die "Couldn't access geneace\n";
  while (<FH>){
    chomp($_);
    if ($_ =~ /^\"(.+)\"\s+\"(.+)\"\s$/) {
      $allele = $1;
      $cds = $2;
      print LOG  "WARNING: Allele $allele has predicted gene but is NOT connected to parent sequence\n";
      $allele_errors++;
      if ($ace){
        get_parent_seq($cds, $allele);
      }
    }
    if ($_ =~ /\"(.+)\"\s+\"(.+)\"\s+\"(.+)\"/){
      $allele = $1;  
      $cds = $2;
      $seq = $3;
      if ($seq eq $cds){
        print LOG "ERROR: Allele $allele has an incorrect parent sequence ($seq) with respect to its predicted gene ($cds)\n";
        $allele_errors++;
        if ($ace){
          print ACE "\n\nAllele : \"$allele\"\n";
          print ACE "-D Sequence \"$seq\"\n";
          &get_parent_seq($cds, $allele);
        }
      }
      if ($seq ne $cds && $seq !~ /SUPERLINK.+/){
        $parent=get_parent_seq($cds, $allele, "getparent");
        if ($parent ne $seq){
          print LOG "ERROR: Allele $allele has an incorrect parent sequence ($seq) with respect to its predicted gene ($cds)\n";
          $allele_errors++;
          if ($ace){
            print ACE "\nAllele : \"$allele\"\n";
            print ACE "-D Sequence \"$seq\"\n";
            print ACE "Sequence \"$parent\"\n";
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

#################################################################

sub check_missing_allele_method {
  my ($def, $dir) = @_;
  my ($allele, $tag);
      
  open (FH, "echo '$def' | tace $dir | ") || die "Couldn't access geneace\n";
  while (<FH>){
    chomp($_);
    print $_, "\n";
    if ($_ =~ /^\"(.+)\"\s+\"(.+)\"/){
      $allele_errors++;
      $allele = $1;
      $tag = $2;
      if ($tag eq "KO_consortium_allele"){$tag = "Knockout_allele"}
      if ($tag eq "Transposon_insertion"){$tag = "Transposon_insertion"}
      if ($ace){output($allele, $tag, "ace")}
      else {output($allele, $tag)}		    
    }
    if ($_ =~ /^\"(.+)\"/){
      print $1, "\n";
      $allele = $1;
      if ($ace){output($allele, "Allele", "ace")}
      else {print LOG "ERROR: Allele $allele has no Method \"Allele\"\n"}		    
    }
  }
  sub output {
    my ($allele, $tag, $ace) = @_;
    print LOG "ERROR: Allele $allele has no Method $tag\n";
    if ($ace ne ""){
      print ACE "\n\nAllele : \"$allele\"\n";
      print ACE "Method \"$tag\"\n";
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

  my $cgc_approved_loci=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/cgc_approved_loci.def" quit
EOF
  
  my %cgc_loci=cgc_loci($cgc_approved_loci, $default_db);
  my $e;

  foreach $strain (@strains){
    if (!$strain->Location){
      print LOG "WARNING: Strain $strain has no location tag\n";
      $strain_errors++; 
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
		foreach $e (@loci){
		  if ($cgc_loci{$e}){
		    print LOG "WARNING: CGC Strain $strain has sequence_name $_ in Genotype, which can now become $e\n";
		    $strain_errors++; 
		  }
                }  
	      }
	      else {
		foreach $e (@loci){
		  if ($cgc_loci{$e}){  
		    print LOG "WARNING: Non_CGC Strain $strain has sequence_name $_ in Genotype, which can now become $e\n";
		    $strain_errors++; 
                  }
		}  
	      }  
	    }
	  }
	}
      }
    }
  }
  
  my ($locus, %locus_strain, $cds, %locus_cds, $main, $other_name, %other_main);
 
  my $get_loci_in_strain=<<EOF;
  Table-maker -p "/wormsrv1/geneace/wquery/locus_in_strain.def" quit 
EOF

  my $locus_has_other_name_to_cds=<<EOF;
Table-maker -p "/wormsrv1/geneace/wquery/locus_has_other_name_to_cds.def" quit
EOF
  my $locus_to_CDS=<<EOF;
Table-maker -p "/wormsrv1/geneace/wquery/locus_to_CDS.def" quit
EOF

  my $locus_to_trans=<<EOF;
Table-maker -p "/wormsrv1/geneace/wquery/locus_to_Transcripts.def" quit
EOF
  
  open (FH1, "echo '$get_loci_in_strain' | tace $default_db | ") || die $!;
  open (FH2, "echo '$locus_has_other_name_to_cds' | tace $default_db | ") || die $!;
  open (FH3, "echo '$locus_to_CDS' | tace $default_db | ") || die $!;
  open (FH4, "echo '$locus_to_trans' | tace $default_db | ") || die $!;

  while (<FH1>){
    chomp;
    if ($_ =~ /\"(.+)\"\s+\"(.+)\"/){
      $strain = $1; $locus = $2;
      push(@{$locus_strain{$locus}}, $strain);
    }
  }   

 
  while (<FH2>){
    chomp $_;
    if ($_ =~ /\"(.+)\"\s+\"(.+)\"\s+\"(.+)\"/){
      $main = $1;
      $other_name = $2;
      $other_name =~ s/\\//;
      if ($3){
	$cds = $3;
	push(@{$other_main{$other_name}}, $main, $cds);
      }
      else {
	push(@{$other_main{$other_name}}, $main);
      }
    }
  }

  while (<FH3>){
    chomp $_;
    if ($_ =~ /\"(.+)\"\s+\"(.+)\"/){
      $locus = $1;
      $locus =~ s/\\//;
      $cds = $2;
      push(@{$locus_cds{$locus}}, $cds);
    }
  } 
  
  
  while (<FH4>){
    chomp $_;
    if ($_ =~ /\"(.+)\"\s+\"(.+)\"/){
      $locus = $1;
      $locus =~ s/\\//;
      $cds = $2;
      push(@{$locus_cds{$locus}}, $cds);
    }
  } 
  
  foreach (keys %locus_strain){
    if(exists $other_main{$_} ){ 
      $strain_errors++;
      print LOG "WARNING: $_ (@{$other_main{$_}}->[1]) in genotype of strain @{$locus_strain{$_}} can now be ";
      print LOG "@{$other_main{$_}}->[0] (@{$locus_cds{@{$other_main{$_}}->[0]}}) -> main name\n";
    }
  } 
  print LOG "\nThere are $strain_errors errors in Strain class.\n";
}

#################################################################

sub cgc_loci {
  my ($def, $db) = @_;
  my (@cgc_loci, %cgc_loci);
  open (FH, "echo '$def' | tace $db | ") || die "Couldn't access $db\n"; 
   
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

##############################################

sub chech_reverse_physicals {

  print "\nChecking reverse physicals of gmap marker loci in Geneace:\n\n";
  system ("/wormsrv2/scripts/get_interpolated_gmap.pl -db /wormsrv1/geneace -reverse");
  my $reverse_log = `echo /wormsrv2/logs/reverse_physicals_WS*.$rundate.*`;
  open(IN, $reverse_log) || die $!;
  while(<IN>){
    print LOG $_;
    print JAHLOG $_;
  }
}

##############################

sub find_new_loci_in_current_DB{
  my ($def, $db) = @_;
  my $warnings;
  my @genes=();
  my $dir="/nfs/disk100/wormpub/DATABASES/current_DB";
  my $locus_errors=0;

  # open a database connection to current_DB and grab all loci names (excluding polymorphisms)
  my $new_db = Ace->connect(-path  => '/nfs/disk100/wormpub/DATABASES/current_DB',
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
      if ($ace){ 
        print ACE "\n\nLocus : \"$gene\"\n"; 
        print ACE "Pseudogene\n";
      }
    }
  }
}

#############################

sub loci_as_other_name {

  my ($def, $dir, $db) = @_;
  my ($main, $other_name, @exceptions, %exceptions);

  open (FH, "echo '$def' | tace $dir | ") || die "Couldn't access geneace\n";
  while (<FH>){
    chomp $_;
    if ($_ =~ /^\"/){
      $_ =~ s/\"//g;
      $_ =~ /(.+)\s+(.+)/;     
      $main = $1;
      $other_name = $2;
      $other_name =~ s/\\//g;
      $other_name = $db->fetch('Locus', $other_name); 
      
      #######################################################    
      # hard coded loci for no main name / other_name merging 
      #######################################################
      @exceptions = 
      qw (aka-1 cas-1 clh-2 clh-3 ctl-1 ctl-2 egl-13 evl-20 gst-4 mig-1 sle-1 slo-1 rap-1 rpy-1 dmo-1 mod-1
          old-1 plk-1 ptp-3 rab-18 rsp-1 rsp-2 rsp-4 rsp-5 rsp-6 sca-1 sus-1);

      foreach (@exceptions){$exceptions{$_}++};  

      if ($other_name){
	if ($exceptions{$main}){
	  print LOG "INFO: $main has $other_name as Other_name...$other_name is still a separate Locus object (exception)\n";
        }
        else {
	  $locus_errors++;
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

#############################

sub loci_point_to_same_CDS {
 
   my ($def, $dir) = @_;
   my %CDS_loci;

   open (FH, "echo '$def' | tace $dir |") || die "Couldn't access geneace\n";
   while (<FH>){
     chomp $_;
     if ($_ =~ /\"(.+)\"\s+\"(.+)\"/){
       #$_ =~ s/\"//g;
       #$_ =~ /(.+)\s(.+)/;
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

  #if ($ace){open (LC, ">>$acefile") || die "Can't write to file!"}
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
	  print ACE "\n\nGene_Class : \"$gc\"\n";
          print ACE "Loci \"$_\"\n";
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
    if ($ace){
      print ACE "\n\nLocus : \"$locus\"\n";
      print ACE "Gene\n";
    }
    $locus_errors++;
  }

  # test for Gene AND !Species 
  if(defined($locus->at('Type.Gene')) && !defined($locus->at('Species'))){
    $warnings .= "ERROR 4: $locus has a 'Gene' tag but not a 'Species' tag\n";
    if ($ace){
      print ACE "\n\nLocus : \"$locus\"\n";
      if ($locus !~ /Cb-|Cr-|Cv/){
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
    $locus_errors++;
  }

  # test for Gene AND !CGC_name AND !Non_CGC_name
  if(defined($locus->at('Type.Gene')) && !defined($locus->at('Name.CGC_name'))
     && !defined($locus->at('Name.Non_CGC_name'))){					       
    $warnings .= "ERROR: $locus has a 'Gene' tag but not a 'CGC_name' or 'Non_CGC_name' tag\n";
    $locus_errors++;
  }


  # test for !Gene AND Gene_class 
  if(!defined($locus->at('Type.Gene')) && defined($locus->at('Name.Gene_class'))){
    $warnings .= "ERROR 5: $locus has a 'Gene_class' tag but not a 'Gene' tag\n";
    print ACE "\n\nLocus : \"$locus\"\n";
    print ACE "Gene\n";
    $locus_errors++;
  }
  
  # test for no Gene tag AND Genomic_sequence tag
  if(!defined($locus->at('Type.Gene')) && defined($locus->at('Molecular_information.Genomic_sequence'))){
    $warnings .= "ERROR 6: $locus has 'Genomic_sequence' tag but no 'Gene' tag\n";
    if ($ace){
      print ACE "\n\nLocus : \"$locus\"\n";  
      print ACE "Gene\n";
    }
    $locus_errors++;
  }

  # test for Genomic_sequence tag but no value   
  if(defined($locus->at('Molecular_information.Genomic_sequence')) && !defined($locus->Genomic_sequence)){
    $warnings .= "ERROR 7: $locus has 'Genomic_sequence' tag but no associated value\n";
    if ($ace){
      print ACE "\n\nLocus : \"$locus\"\n"; 
      print ACE "-D Genomic_sequence\n";
    } 
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


  # Test for Polymorphisms with no P in their title
  #if(defined($locus->at('Type.Polymorphism'))){
   # if($locus !~ /P/){
    #  $warnings .= "ERROR 14: $locus has no 'P' in title\n";
     # $locus_errors++;
    #}
  #}
  # Look for Gene_class tag in non-gene objects 
  if(!defined($locus->at('Type.Gene'))){
    if(defined($locus->at('Name.Gene_class'))){
      $warnings .= "ERROR 15: $locus has Gene_class tag but it is not a gene!\n";
      $locus_errors++;
      if ($ace){
        print ACE "\n\nLocus : \"$locus\"\n"; 
        print ACE "Gene\n";
      } 
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
	  $warnings .= "ERROR 17: $locus has an Other_name: $other_name which is also a separate locus object\n";
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
  open(JAHLOG, ">>$jahlog") || die "Can't open $jahlog\n";
  print JAHLOG "This mail is generated automatically for CGC on $rundate\n"; 
  $JAHmsg = "This mail is generated automatically for CGC on $rundate\n"; 
  print JAHLOG "If you have any queries please email ck1\@sanger.ac.uk or krb\@sanger.ac.uk\n\n";
  $JAHmsg .= "If you have any queries please email ck1\@sanger.ac.uk or krb\@sanger.ac.uk\n\n";
  print JAHLOG "=========================================================================\n";
  $JAHmsg .=  "=========================================================================\n";
  
  # create separate log with errors for Erich
  $erichlog = "/wormsrv2/logs/geneace_check.erichlog.$rundate.$$";
  open(ERICHLOG,">$erichlog") || die "cant open $erichlog";
  print ERICHLOG "$0 started at ",`date`,"\n";
  $Emsg = "$0 started at ".`date`."\n";
  print ERICHLOG "This mail is generated automatically for Caltech\n";
  $Emsg .= "This mail is generated automatically for Caltech\n";
  print ERICHLOG "If you have any queries please email ck1\@sanger.ac.uk or krb\@sanger.ac.uk\n\n";
  $Emsg .= "If you have any queries please email ck1\@sanger.ac.uk or krb\@sanger.ac.uk\n\n";   
  print ERICHLOG "================================================================================================\n";
  $Emsg .= "================================================================================================\n";
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

  my $dir="/nfs/disk100/wormpub/DATABASES/current_DB";
  
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
               lab or laboratory
               strain
               rearr or rearrangement
               seq or sequence

            For example: -c allele -c locus OR -class seq -class rearr

B<-databse:> 
            Allows specifying path to a specific database.
            Default database path is /wormsrv1/geneace without this option.

            For example: -db /wormsrv2/autoace or -database /wormsrv1/someone/Test_DB
 

B<-ace:>     
            Allows generating ace file for fixing erros spotted by 
            this checking script.
            Default location and filename of ace file:
            /wormsrv2/logs/geneace_check_ACE.rundate.processid
                      
            For example: -a or -ace


=head3 <RUN geneace_check.pl>

            Geneace check is now set to run on Sundays.
            ##Run Geneace check over weekend
            20 7 * * 0 /wormsrv2/scripts/geneace_check.pl -a &
