#!/software/bin/perl -w
# 
# variation_auto_fix.pl
#
# Script to fix consistencies highlighted by geneace_check.pl
#

use strict;
use lib $ENV{"CVS_DIR"};
use Wormbase;
use Ace;
use Ace::Object;
use Carp;
use Getopt::Long;
use GENEACE::Geneace;
use File::Path;

###################################################
# command line options                            # 
###################################################
my ($help, $debug, $test, $class, @classes, $database, $ace, $verbose);
my @skip_methods;

GetOptions ('help'          => \$help,
            'debug=s'       => \$debug,
	    'class=s'       => \@classes,
	    'database=s'    => \$database,
            'ace'           => \$ace,
	    'verbose'       => \$verbose,
            'test'          => \$test,
            'skipmethod=s@' =>  \@skip_methods,
	   );

###################################################
# Miscellaneous important variables               # 
###################################################
my $wb = Wormbase->new(-test => $test, -debug => $debug );
# choose database to query: default is geneace
$database = '/nfs/production/flicek/wormbase/wb/DATABASES/geneace' unless $database;
print "Using database $database.\n\n";

my $tace = $wb->tace;          # tace executable path
my $curr_db = '/nfs/production/flicek/wormbase/wb/DATABASES/current_DB'; # Used for some cross checking with geneace
my $def_dir = "$ENV{CVS_DIR}/../wquery/geneace";                          # where lots of table-maker definitions are kept

my $rundate = $wb->rundate;                                # Used by various parts of script for filename creation
my $maintainers = 'hinxton\@wormbase.org';
my $log_dir = "$database/logs";                            # some of the many output files are put here (ar2)
my $log;                                                   # main log file for most output
my (%L_name_F_WBP, %L_name_F_M_WBP);                       # hashes for checking Person and Author merging?


##################################################
# other variables                                #
##################################################

# Display help if required
if ($help){&usage('Help')}

# Use debug mode?
if($debug){
  print "\nDEBUG = \"$debug\"\n";
  ($maintainers = "$debug" . '\@sanger.ac.uk');
}

# Open file for acefile output?
my $acefile;
if ($ace){
  mkdir "$database/CHECKS" unless ( -e "$database/CHECKS");
  $acefile = "$database/CHECKS/Variation_fixes.$rundate.$$.ace";
  open(ACE, ">>$acefile") || croak $!;
  system("chmod 777 $acefile");
}

my $next_build_ver = $wb->get_wormbase_version() + 1 ; # next build number

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
my $db = Ace->connect(-path  => $database,
		      -program =>$tace) || do { print LOG "Connection failure: ",Ace->error; die();};


my $ga = init Geneace($wb);
# hash for converting locus/seq. names <-> gene_id
my $Gene_info = $ga -> gene_info($database, "seq2id");
my %Gene_info = %{$Gene_info};

print  "Checking Variations and attempting to fix things:\n";
print  LOG "Checking Variations and attempting to fix things:\n";


# check that there is a Version tag
foreach my $variation ($db->fetch(-query=>'Find Variation WHERE NOT Laboratory')){
    print LOG "ERROR: $variation ($variation has no Laboratory\n";
    unless (defined $variation->Public_name){
	next;
    }
    my $allele_name = $variation->Public_name;
    if (($allele_name->name =~ /([[:alpha:]]+)\d+/) || ($allele_name->name =~ /([[:alpha:]]+)\d+/)) {
	next if $1 =~ /^\S{2}P/;
	next if $1 =~ /^CE/;
	next if $1 =~ /WBVar/;
	my $prefix = $1;
	print "$prefix\n";
	my @lab = $db->fetch(-query=>'Find Laboratory WHERE Allele_designation = $prefix');
	if (defined $lab[0]) {
	    my $lab_fix=$lab[0]->name;
	    print "$variation->name should be connected to $lab_fix\n";
	    print ACE "Variation $variation->name\nLaboratory $lab_fix\n\n" if ($ace);
	}
	else {
	    print "Cannot auto fix Lab connection for $variation\n";
	}
    }
}

#######################################
# Tidy up and mail relevant log files #
#######################################

$db->close;
print LOG "\nEnded at ",`date`,"\n";
close(LOG)
    or warn $! ? "Error closing log: $!"
    : "Exit status $? log close";
print "Exit status $? log close\n\n" if ($debug);


close(ACE) if ($ace);
print "Output here $acefile\n" if ($ace);
# email everyone specified by $maintainers
$wb->mail_maintainer('geneace_check: SANGER',$maintainers,$log);

exit(0);
#--------------------------------------------------------------------------------------------------------------------









                          #############################################################
                          #         SUBROUTINES FOR -class evidence option            #
                          #############################################################

sub check_evidence {

  my $WBPerson_F_M_L_names="Table-maker -p \"$def_dir/WBperson_first_middle_last_names.def\"\nquit\n";
  my @WBPerson = &process_WBPerson_names($WBPerson_F_M_L_names, $curr_db);

  print "\n\nChecking misuse of Evidence and converting Author to Person / Non-Person to Author:\n\n";
  print LOG "\n\nChecking misuse of Evidence and converting Author to Person / Non-Person to Author:\n";
  print LOG "-----------------------------------------------------------------------------------\n";


# dump flat files with time stamps
my $dump_dir = "$database/CHECKS";

my $command=<<END;
Find Gene * 
show -a -T -f $dump_dir/gene_dump.ace

Find Variation *
show -a -T -f $dump_dir/variation_dump.ace

Find strain *
show -a -T -f $dump_dir/strain_dump.ace

Find gene_class *
show -a -T -f $dump_dir/geneclass_dump.ace

Find 2_point_data *
show -a -T -f $dump_dir/2_pt_dump.ace

Find Multi_pt_data *
show -a -T -f $dump_dir/multi_pt_dump.ace

Find Pos_neg_data *
show -a -T -f $dump_dir/posneg_dump.ace

Find Laboratory *
show -a -T -f $dump_dir/lab_dump.ace

quit
END

  open (DUMP, "| $tace $database") || warn "Failed to connect to Geneace";
  print DUMP $command;
  close DUMP;

  `cat $dump_dir/gene_dump.ace $dump_dir/variation_dump.ace $dump_dir/strain_dump.ace $dump_dir/2_pt_dump.ace $dump_dir/multi_pt_dump.ace $dump_dir/posneg_dump.ace $dump_dir/geneclass_dump.ace $dump_dir/lab_dump.ace > $dump_dir/class_dump.ace`;

  my @dumps = qw (gene_dump.ace allele_dump.ace strain_dump.ace geneclass_dump.ace 2_pt_dump.ace multi_pt_dump.ace posneg_dump.ace class_dump.ace lab_dump.ace);

  foreach (@dumps){system ("chmod 777 $dump_dir/$_")}

  open(IN, "$dump_dir/class_dump.ace") || warn $!;

# look for person/author names that needs to be converted in 8 classes (regardless of which tag as the scripts greps from string pattern

  my $evid_errors = 0;
  my $updates = 0;
  my $info_num = 0;
  my (@counters, $class_obj, $class, $obj, $tag, $ori, $b4_evi, $name, $paper, $author, $last_name, $initials);

  while (<IN>){
    chomp;
    if ($_ =~ /^(Gene) : \"(.+)\" -O .+/){$class = $1; $obj = $2}
    if ($_ =~ /^(Variation) : \"(.+)\" -O .+/){$class = $1; $obj = $2}
    if ($_ =~ /^(Strain) : \"(.+)\" -O .+/){$class = $1; $obj = $2}
    if ($_ =~ /^(Gene_class) : \"(.+)\" -O .+/){$class = $1; $obj = $2}
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
    if ($_ =~ /(\w+)\s+.+Paper_evidence -O .+\"(.+)\" -O.+/){
      $ori = $_;
      $tag = $1;
      $paper = $2;
      if ($paper !~ /WBPaper\d+/){
        $evid_errors++;
	print LOG "\nERROR: $class $obj has Paper $paper under main tag $tag\n";
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
  print LOG "\n\nThere are $evid_errors Evidence errors in 8 classes checked\n";
  print LOG "\n$updates Authors can be converted to Persons\n";
  print LOG "\n$info_num Authors are not Persons\n" if $verbose;
  system ("rm -f $dump_dir/*_dump.ace");
}




                          ###########################################################
                          #         SUBROUTINES FOR -class allele option            #
                          ###########################################################


sub process_allele_class{

  print"\n\nChecking Allele class for errors:\n";
  print LOG "\n\nChecking Allele class for errors:\n";
  print LOG "---------------------------------\n";

  # make hash of allele to lab connections
  my %allele2lab;
  my $def="Table-maker -p \"$def_dir/allele_designation_to_LAB.def\"\nquit\n";

  open (FH, "echo '$def' | $tace $database | ") || warn "Couldn't access $db\n";
  while (<FH>) {
    print;
    chomp;
    if ($_ =~ /^\"(.+)\"\s+\"(.+)\"/) {
      $allele2lab{$2} = $1	# $2 is allele_designation $1 is lab designation	
    }
  }

  my $query = "find Variation Allele";
  foreach my $meth (@skip_methods) {
    $query .= " AND Method != \"$meth\"";
  }
  my $alleles_it = $db->fetch_many(-query => "$query");

  # now loop through all alleles looking for problems
#  foreach my $allele (@alleles) {
  while( my $allele = $alleles_it->next){
    print "$allele\n" if ($verbose);

    # check allele has no laboratory tag
    if (!defined $allele->Laboratory ) {
      
      if (!$ace) {
	print LOG "ERROR: $allele has no Laboratory tag\n";
      } else {
	print LOG "ERROR(a): $allele has no Laboratory tag\n";

	# try to find lab designation for alleles with CGC-names (1 or 2 letters and then numbers)
	if ($allele =~ /^([a-z]{1,2})\d+$/) {
	  my $allele_name = $1;	  

	  print ACE "\nVariation : \"$allele\"\n";
	  print ACE "-D Laboratory\n";
	  print ACE "\nVariation : \"$allele\"\n";
	  print ACE "Laboratory \"$allele2lab{$allele_name}\"\n";
	  next;
	}	
      }
    }
  

    # check that Lab designation (when present) is correct (for CGC named alleles)
    if (defined $allele->Laboratory) {
      my $lab = $allele->Laboratory;
      if ($allele =~ m/^([a-z]{1,2})\d+$/) {
	my $allele_prefix = $1;
	if ( $allele2lab{$allele_prefix} ne $lab ) {
	  print LOG "ERROR: $allele is connected to $lab lab which is incorrect according to $allele2lab{$allele_prefix} lab info\n";
	}
      }
    }


    # check for Missense tag but no value
    if (!defined $allele->Missense) {
      my $missense = $allele->at('Description.Missense');
      if (defined($missense)) {
	print LOG "ERROR: $allele has a missense tag but does not have an associated value\n";
      }
    }
    # now check for structure of Missense field
 #   else {
 #     my $missense = $allele->Missense;
 #     if ($missense !~ m/^[A-Z]\(\d+\) to [A-Z]$/) {
#	print LOG "ERROR: $allele has an incorrect Missense value ($missense)\n";
 #     }
 #   }

    # Check for Public_name tag missing
    if (!defined($allele->Public_name)) {
	print LOG "ERROR: $allele has no Public_name tag\n";
    }

    # Check for Status tag missing
    if (!defined($allele->Status)) {
      print LOG "ERROR: $allele has no Status tag\n";
    }

    # Check for SeqStatus tag missing
    if (!defined($allele->SeqStatus)) {
	print LOG "ERROR: $allele has no SeqStatus tag\n" if ($allele->Status eq "Live");
    }

    # Check for method tag missing
    if (!defined($allele->Method)) {
      print LOG "ERROR: $allele has no Method tag\n";
    }

    # check for Method tag which has no value
    if (!defined $allele->Method) {
      my $method = $allele->at('Method');
      if (defined($method)) {
	print LOG "ERROR: $allele has a Method tag but not associated value\n";
      }
    }

    # check that Allele_type is set for all alleles with flanking sequences
    if (defined($allele->Flanking_sequences) && !defined($allele->Type_of_mutation)) {
      print LOG "ERROR: $allele has flanking sequences but no 'Type_of_mutation' tag\n";  
    }
    

    # test allele has method that matches Allele_type, eg, Deletion tag -> Deletion_allele
    if (defined $allele->Type_of_mutation) {
      my $expected_method; 
      my $observed_method;


      my @mut_type = $allele->Type_of_mutation;
      if ( scalar @mut_type == 1 ) {
	if ($mut_type[0] eq "Deletion" ) {
	  $expected_method = "Deletion_allele";
	} elsif ($mut_type[0] eq "Insertion" && !defined $allele->Transposon_insertion ) {
	  $expected_method = "Insertion_allele";
	} elsif ($allele->Transposon_insertion && $allele->at('Sequence_details.Type_of_mutation.Insertion')) {
	  $expected_method = "Transposon_insertion";
	} elsif ($allele->Type_of_mutation eq "Substitution" ) {
	  $expected_method = "Substitution_allele";
	} elsif ( grep(/Deletion/, @mut_type) and grep(/Insertion/, @mut_type) ) {
	  $expected_method = "Deletion_and_insertion_allele";
	}

	($observed_method = $allele->Method) if (defined $allele->Method);
 
	# does $observed method tag agree with expected method tag (based on Type_of_mutation tag)?
	if ($expected_method ne $observed_method) {
	  if ($ace) {
	    print LOG "ERROR(a): $allele has wrong method ($observed_method): change to $expected_method\n";
	    print ACE "\nAllele : \"$allele\"\n";
	    print ACE "-D Method\n";
	    print ACE "\nAllele : \"$allele\"\n";
	    print ACE "Method \"$expected_method\"\n";
	  } else {
	    print LOG "ERROR: $allele has method $observed_method which might need to be $expected_method\n";
	  }
	}
      }

      # find alleles that have flanking_seqs but no SMAPPED sequence
      if ( $allele ->Flanking_sequences && ! defined $allele ->Mapping_target ) {
	print LOG "ERROR: Allele $allele has Flanking_sequences tag but has no Mapping_target tag\n";
      }
    }  
  }
}


sub usage {
  my $error = shift;
  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

sub create_log_files{

  $log = "$log_dir/geneace_check.$rundate.$$";

  open (LOG, ">$log") or warn "cant open $log";
  print LOG "geneace_check\n";
  print LOG "started at ",`date`,"\n";
  print LOG "=============================================\n";
  print LOG "\n";


  if($ace){
    print LOG "The (a) following ERROR, or UPDT, eg, indicates ace output \n$acefile for direct upload to correct problems.\n\n";
  }
  
}


__END__

=head2 NAME - variation_auto_fix.pl

=head3 <USAGE>


=head2 Options: [help] [debug] [class] [ace]


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
            Current valid options to check these classes are:

               gene
               allele
               laboratory
               strain
               rearrangement
               evidence
               mapping
               gmap
               multipoint

            For example: -class allele -class gene OR -class sequence -class rearrangment

B<-database:>
            Allows specifying path to a specific database.
            Default database path is /nfs/wormpub/DATABASES/geneace without this option.

            For example: -database /wormsrv2/autoace or -database /wormsrv1/someone/Test_DB


B<-ace:>
            Allows generating ace file for fixing erros spotted by
            this checking script.
            Default location and filename of ace file:
            /nfs/wormpub/DATABASES/geneace/CHECKS/geneace_check.rundate.processid.ace
            For example: -ace


B<-verbose:>
            Toggles on extra output.   Useful when running on command line and not on crob job
            For the locus class it will display each locus name as it is processes it and show
            (next to the locus name) a dot for each error it had found in that locus


