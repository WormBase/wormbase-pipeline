#!/software/bin/perl -w
#
# EMBL_Sequencefetch_species.pl
#
# Usage : EMBL_Sequencefetch.pl [-options]
#
# Last edited by: $Author: mh6 $
# Last edited on: $Date: 2012-07-02 11:53:27 $

my $script_dir = $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'};
use strict;
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Modules::Features;
use Species;

######################################
# variables and command-line options #
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $version, $organism, $output, $nolongtext, $dna, $database, $inc, $species, $dump_only, $input, $Type);

GetOptions ("help"       => \$help, #
            "debug=s"    => \$debug, # debug option, turns on more printing and only email specified user.
            "test"       => \$test, # only genomic RNAs will be used for a quick test.
            "verbose"    => \$verbose, # additional printing
            "store:s"    => \$store, # wormbase storable object
	    "organism=s" => \$organism, # Specify an organism
	    "dna"        => \$dna, # Create a seperate DNA file.
	    "database=s" => \$database, #database are you downloading sequence data for.
	    "inc"        => \$inc, #just query emblnew.
	    "dump_only"  => \$dump_only, # only dumps the Transcript data from the primary databases, not EMBL connection is made.
	    "nolongtext" => \$nolongtext, # don't dump longtext
	    "input=s"    => \$input, #EMBL flat file to be parse if doing manually.
	    "Type=s"     => \$Type, #If you know the type of sequence ie. EST, mRNA state this here (Used in conjunction with -input).
	   );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
			     -organism => $species,
			   );
}

#################################
# Set up some useful stuff      #
#################################
# establish log file.
my $log = Log_files->make_build_log($wormbase);
my $mfetch = "/software/bin/mfetch"; #mfetch script
my $tace = $wormbase->tace; # TACE PATH
# output dirs
my   $output_dir = $wormbase->database('camace')."/EMBL_sequence_info";
$wormbase->run_command ("mkdir $output_dir", $log) if (!-e $output_dir);

my %TSLs = $wormbase->TSL;
Features::set_tsl(\%TSLs);

# main stuff goes here
my (%acc_sv2sequence,%method2sequence,%feature2seq);
my $molecule;
my @molecules;
my @species;
my %acc2molecule;
my $new_sequence;
my $acefile;
my $dnafile;
my $Longtextfile;
my $miscdir = glob('~wormpub/BUILD_DATA/MISC_DYNAMIC');

#Molecules types
#molecule type to Type tag data hash. key:molecule value:type
my %molecule2rnatype = ( 'genomic RNA' => "ncRNA",
			 'mRNA' => "ESTormRNA",
			 'other RNA' => "ncRNA",
			 'pre-RNA' => "ESTorRNA",
			 'rRNA' => "rRNA",
			 'snRNA' => "snRNA",
			 'unassigned RNA' => "ncRNA",
			 'tRNA' => "tRNA",
			 'transcribed RNA' => "ncRNA",
		       );

my %species2taxonid = (
		       'elegans'           => "6239",
		       'briggsae'          => "6238",
		       'remanei'           => "31234",
		       'japonica'          => "281687",
		       'brenneri'          => "135651",
		       'pristionchus'      => "54126",
		       'heterorhabditis'   => "37862",
		       'brugia'            => "6279",
		       'cangaria'          => '96668',
		      );

if ($test) {
  @molecules  = ('mRNA');
  print "\n---------------------------------------------------------------------------------------------------------\n
WARNING: You are running in test mode, only molecule(s) [@molecules] will be returned!!!!!!!\n
---------------------------------------------------------------------------------------------------------\n";
}
else {
  @molecules  = ("genomic RNA","mRNA","other RNA","pre-RNA","rRNA","snRNA","unassigned RNA","tRNA","transcribed RNA");
}

# incremental mode
if ($inc) {
  print "\n---------------------------------------------------------------------------------------------------------\n
WARNING: You are retrieving EMBL sequences as an incremental update, to run a full re-sync, remove the -inc.\n
---------------------------------------------------------------------------------------------------------\n";
}
# Full mode
if ((!defined$inc) && (!defined$input)) { 
print "\n---------------------------------------------------------------------------------------------------------\n
WARNING: You are retrieving EMBL sequences as a full update (This will take a long time)!\nTo run an incremental update, add -inc to the command line.\n
---------------------------------------------------------------------------------------------------------\n";
}

##########################
# MAIN BODY OF SCRIPT
##########################

#######################################################################################
# Default option                                                                      #
# Update sequence data for all WormBase species (depending upon options) using mfetch #
#######################################################################################
if (!defined$input){
  my @organism;
  if (!defined($organism)) {
    @organism = (keys %species2taxonid);
  }
  
  elsif ($organism) {
    push (@organism, $organism);
  }


  foreach my $organism(@organism) {
    $log->write_to("============================================\nProcessing: $organism\n============================================\n");
    my $sourceDB;

    if ($database){
	   $sourceDB = $database;
	}
    elsif ($organism eq "elegans") {
      $sourceDB = $wormbase->database('camace');
    }
    elsif ($organism eq "heterorhabditis") {
      $sourceDB = glob "~wormpub/DATABASES/heterorhabditis";
    }
    elsif ($organism eq "pristionchus"){
      $sourceDB = glob "~wormpub/DATABASES/pristionchus";
    }
    else {$sourceDB = $wormbase->database($organism);}
    # Fetch sequence data from primary database.
    $log->write_to("Fetching sequence data from $sourceDB:\n");
    &fetch_database_info ($sourceDB) if ((-e $sourceDB) & (!defined $dump_only));
    foreach $molecule(@molecules) {
      $log->write_to("============================================\nProcessing: $molecule\n");
      my @entries;
      # fetch a list of the organism sequencae info from embl or emblnew.
      $log->write_to("Fetching Sequence data from EMBL: $organism - $molecule\n");
      &get_embl_data (\@entries, $molecule, $organism) if (!defined $dump_only);
      # Retrieve the full object and process.
      if (scalar@entries > 0) {
	$log->write_to("Fetching NEW EMBL data: $organism - $molecule\n");
	&get_new_data (\@entries, $molecule, $organism) if (!defined $dump_only);
	# Load the data back to the primary database.
	$log->write_to("Loading $organism data into primary database\n");
	&load_data ($sourceDB, $molecule, $organism) if (!defined $dump_only);
	$log->write_to("Processed: $molecule\n============================================\n");
      }
      # if there isn't any, just print this line in the log.
      if (scalar@entries < 1) {
	$log->write_to("No NEW EMBL data:$organism - $molecule\n");
      }
    }
    
    $log->write_to("============================================\nProcessed: $organism\n============================================\n\n");
  }
}

####################################
# Optionally just parse input file #
####################################

elsif ($input) {
  print "You are generating Sequence data and Features from a supplied EMBL flat file...is this correct?\n\n";
  
  if (!defined $Type) {
    print "ERROR: You need to specify a data Type (eg. EST mRNA) using the command line option -Type\n\nDo you want to assign one now and continue? (y or n)\n";
    my $answer=<STDIN>;
    if ($answer eq "y\n") {
      print "Please input a data Type:";
      $Type=<STDIN>;
    }
    if ($answer eq "n\n") {
      die "Failed to process $input as a data type was not specified\n To avoid this issue please use -Type in conjunction with -input on command line\n"
    }
  }
    &generate_data_from_flat ($input) ;
}


$log->write_to("\nRefreshed........ ;) \n\n");
$log->write_to ("WORK DONE------------FINISHED\n\n");
$log->mail();
exit(0);





                      ###################################################
                      #                 SUBROUTINES                     #
                      ###################################################

###############################################################################
# fetch_database_info                                                         #
# this subroutine connects to the acedb database and pulls out sequqnce info. #
###############################################################################
sub fetch_database_info {
  my $sub_sourceDB = shift;
  my $def_dir = $wormbase->database('camace')."/wquery";
  my $tablemaker_query =  "${def_dir}/SCRIPT:estfetch.def";
  my $command = "Table-maker -p $tablemaker_query\nquit\n";

  open (TACE, "echo '$command' | $tace $sub_sourceDB |");
  while (<TACE>) {
    my $status;
    chomp;
    s/\"//g;
    next if ($_ eq "");
    next if (/acedb\>/);

    if ((/^\S+\t(\S+).(\d+)\t(\S+)\t(\S+)/) or (/^\S+\t(\S+).(\d+)\t(\S+)\t/)) {
      unless($1 and $2 and $3) {
	$log->error("$_\n");
	next;
      }
      # split the line into various fields
      my $sequenceacc = $1;
      my $acc_sv = $2;
      my $method = $3;
      my $feature;

      if (defined$4) {
	$feature = $4;
#	print "Sequence: $sequenceacc has the feature $4\n\n"
      }
      if (defined$feature) {
	$status = "yes";
      } 
      else {
	$status = "No";
      }

      $acc_sv2sequence{$sequenceacc} = $acc_sv;
      $method2sequence{$sequenceacc} = $method;
      $feature2seq{$sequenceacc} = $status;
    }
  }
  close TACE;
  $log->write_to("DB: Got it!!\n\n");
}

################################################################################################
# get_embl_data                                                                                #
# This subroutine retrieves a list of accessions and compares to what is in the acedb database #
# New data is the retrieved and stored for processing.                                         #
################################################################################################

sub get_embl_data {
  my $entries = shift;
  my $molecule = shift;
  my $species = shift;
  if ($inc) {
    open (SEQUENCES, "$mfetch -d emblnew -i \"mol:$molecule&txi:$species2taxonid{$species}\" |");
    print "Incremental update selected <mfetch -d emblnew -i \"mol:$molecule&txi:$species2taxonid{$species}\">\n";
  }
  else {
    open (SEQUENCES, "$mfetch -d embl -i \"mol:$molecule&txi:$species2taxonid{$species}\" |");
    print "Default Full update selected <mfetch -d embl -i \"mol:$molecule&txi:$species2taxonid{$species}>\"\n";
    }
    while (<SEQUENCES>) {
      chomp;
      if (/no match/) {
	$log->write_to("No entries were retrieved for Species:$species in Molecule type:$molecule\n");
      }
      #eg. ID   AX254400; SV 1; linear; unassigned RNA; PAT; INV; 1161 BP.
      # drop the -f switch from mfetch gives AF273797.1
      elsif (my($accR,$svR) = /(\S+)\.(\d+)/) {
	#elsif (my($accR,$svR) = /^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+/) {
	print "$accR.$svR\n" if ($verbose);
	
	# Store data about the sequence is datahash.
	# add to hash. Accession is key, molecule is value
	$acc2molecule{$accR} = ($molecule2rnatype{$molecule});

	#Is this entry already in the database?
	if (defined $acc_sv2sequence{$accR}) { # If yes
	  print "PRIMARY:$acc_sv2sequence{$accR} EMBL:$svR\n" if ($verbose);
	  if ($acc_sv2sequence{$accR} == $svR) { # If the SV is the same don't do anything
	    $log->write_to("$accR Sequence Versions match PRIMARY:$acc_sv2sequence{$accR} EMBL:$svR\n") if $verbose;
	  } else {		# If the SV is different update it!!
	    push(@$entries,"$accR\.$svR");
	    $log->write_to("$accR - Sequence Versions don\'t match PRIMARY:$acc_sv2sequence{$accR} EMBL:$svR\n");
	    $log->write_to("Updated entry $accR fetch it!\n");
	  }
	  if ($feature2seq{$accR} eq "Yes") { # be careful of feature_data.
	    print "Warning $accR has associated Feature_data....this will need updating\n";
	    $log->write_to("Warning $accR has associated Feature_data....this will need updating\n");
	  }
	} else {
	  $log->write_to("New entry $accR fetch it!\n");
	  push(@$entries,"$accR\.$svR");
	}
      }
    }				# while SEQUENCE
    close (SEQUENCES);
    $log->write_to(scalar@$entries." entrie(s) in the $species $molecule subclass need to be retrieved from embl\n");
    $log->write_to("EMBL: Got it!!\n");
}

########################################################################################
# get_new_data                                                                         #
# This subroutine queries EMBL for all $subentries passed from previous subroutine and #
# writes out the associated output files.                                              #
########################################################################################

###################################################################################
# Let the real work commence - fetch the entire entry for new sequences from embl #
###################################################################################
# examples of problems/exceptions that need to be captured and pparse differently.
#D36756.1   ID   D36756;   SV 1; linear; mRNA;        EST; INV; 360  BP. yk36g11.5 exception
#CB400075.1 ID   CB400075; SV 1; linear; mRNA;        EST; INV; 559  BP. OST exception
#M89002.1   ID   M89002;   SV 1; linear; mRNA;        EST; INV; 396  BP. (CEL12H5) Chris Martin EST exception
#AU115113.2 ID   AU115113; SV 2; linear; mRNA;        EST; INV; 236  BP. yk726a3.3 exception
#CB393073.1 ID   CB393073; SV 1; linear; mRNA;        EST; INV; 331  BP. OSTR113E12_1 OST exception
#T00346.1   ID   T00346;   SV 1; linear; mRNA;        EST; INV; 318  BP. (CEESE35) Stratagene est exception
#AF273829.1 ID   AF273829; SV 1; linear; mRNA;        STD; INV; 1183 BP. Standard mRNA NDB
#L25083.1   ID   L25083;   SV 1; linear; genomic DNA; STD; INV; 412  BP. Standard DNA NDB

sub get_new_data {
  my $subentries = shift;
  my $submol = shift;
  my $suborganism = shift;
  my $submol_mod;
  if ($submol eq "genomic RNA") {$submol_mod = "genomic_RNA";}
  elsif ($submol eq "other RNA") {$submol_mod = "other_RNA";}
  elsif ($submol eq "unassigned RNA") {$submol_mod = "unassigned_RNA";}
  elsif ($submol eq "transcribed RNA") {$submol_mod = "transcribed_RNA";}
  else {$submol_mod = $submol;}
  #Output files
  $acefile = "$output_dir/new_${suborganism}_$submol_mod.ace";
  $wormbase->run_command ("rm $acefile", $log) if (-e $acefile);

  if ($dna) {
    $dnafile = "$output_dir/new_${suborganism}_$submol_mod.dna";
    $wormbase->run_command ("rm $dnafile", $log) if (-e $dnafile);
  }

  $Longtextfile = "$output_dir/new_${suborganism}_${submol_mod}_longtext.txt";
  $wormbase->run_command ("rm $Longtextfile", $log) if (-e $Longtextfile); 
  
  #open output file for full entries.
  open (OUT_ACE,  ">$acefile");
  print "$acefile opened!!\n" if $debug;
  open (OUT_DNA,  ">$dnafile") if defined($dna);
  open (OUT_LONG, ">$Longtextfile") unless (defined $nolongtext);

  foreach $new_sequence (@{$subentries}) {
    my ($seq,$status,$protid,$protver,@description,$idF,$svF,$idF2,$type,$def,$sequencelength);

    #     $new_sequence = "DQ342049"; #get 1 entry DQ342049.1
    open (NEW_SEQUENCE, "$mfetch -d embl -v full $new_sequence |");
    while (<NEW_SEQUENCE>) {
      #my $idF2;
      #Extract ID and Sequence Version.
      if ((/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+.+mRNA.+(EST)/) or (/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+.+mRNA.+(STD)/)) { #mRNA or EST
	#        if (/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+/) {
	print OUT_LONG "\nLongText : \"$1\"\n" unless (defined $nolongtext);
	$idF = $1;
	$svF = $2;
	$type = $3;
      } elsif (/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+/) { #Non-EST/mRNA entries
	print OUT_LONG "\nLongText : \"$1\"\n" unless (defined $nolongtext);
	$idF = $1;
	$svF = $2;
      }
      my $DEline1;
      if (/^DE/) {
	if (/^DE.+(yk\d+\S+\d+)\s+:\s+(\d)/) { #DE   Caenorhabditis elegans cDNA clone yk181f5 : 5' end, single read.
	  $idF2 = "$1".".$2";
	  $DEline1 = "1";
	} elsif (/^DE\s+(\S+)\s+Chris Martin.+/) { #DE   CEL01A1S1 Chris Martin sorted cDNA library Caenorhabditis elegans cDNA
	  $idF2 = $1;
	  $DEline1 = "1";
	} elsif (/DE\s+(OST\S+)\s+/) {
	  $idF2 =  $1;
	  $DEline1 = "1";
	} elsif (/DE\s+.+clone\s+(\S+)\s+similar/) { #DE   cDNA clone CEESE35 similar to ATP synthase lipid binding protein P1, mRNA
	  $idF2 =  $1;
	  $DEline1 = "1";
	} else {
	  $idF2 = $idF;		# if ($DEline1 eq "0");
	}
      }
      # grab various details out of EMBL entry
      if (/^DE\s+(.+)/) {
	$def = $1;
	# remove any offending '>' from def line. This is required by transcriptmasker.pl
	$def =~ s/\>//g;
	push(@description,"$def");
      }
      # Get protein ID info.
      if (/^FT\s+\/protein_id=\"(\S+)\.(\d+)\"/) {
	$protid=$1; $protver=$2;
      }
      # Start of DNA sequence in entry.
      if (/^SQ\s+\Sequence\s+(\d+)/) {
	$sequencelength = $1;
	print "$sequencelength\n" if ($verbose);
	print OUT_DNA "\n>$idF2\n" if $dna;
	print OUT_ACE "\nDNA \: \"$idF2\"\n";
	$status = 1;		#set status to 1 for identifying start of sequence.
	print "$status\n" if ($verbose);
      }
      print OUT_LONG "$_" unless (defined $nolongtext); #print out the embl entry line unmodified to the LongText file
      if (/^\s/) {
	s/\s+//g;
	s/\d+//g;
	s/[^acgtn]//ig;
	print OUT_DNA "$_\n" if (($status eq "1") && $dna);
	chomp;
	$seq .= "$_" if ($status eq "1");
	print OUT_ACE "$_\n" if ($status eq "1");
      }
      if (/^\/\//) {		#end of embl entry
	$status = 0;
	$DEline1 = 0;
	print OUT_LONG "***LongTextEnd***\n" unless (defined $nolongtext);
	print "$status\n" if ($verbose);
	print OUT_ACE "\nSequence : \"$idF2\"\n";
	#	print OUT_ACE "DNA $idF2 $sequencelength\n";
	print OUT_ACE "DNA $idF2\n";
	print OUT_ACE "Database EMBL NDB_AC $idF\n";
	print OUT_ACE "Database EMBL NDB_ID $idF\n";
	print OUT_ACE "Database EMBL NDB_SV $idF\.$svF\n";
	print OUT_ACE "DB_annotation EMBL $idF\n";
	print OUT_ACE "Protein_id $idF $protid $protver\n" if (defined$protid && defined$protver);
	if ($suborganism eq "pristionchus") {print OUT_ACE "Species \"Pristionchus pacificus\"\n";}
	elsif ($suborganism eq "heterorhabditis") { print OUT_ACE "Species \"Heterorhabditis bacteriophora\"\n";}
	else {print OUT_ACE "Species \"Caenorhabditis $suborganism\"\n";}
	print OUT_ACE "Title \"@description\"\n";

	# Properties depend on the molecule subdivision of embl that the sequence was fetched from.
	
	#ESTs
	if (defined$type) {
	  if ($type eq "EST") {
	    print OUT_ACE "Properties cDNA cDNA_EST\n";
	  }
	  #mRNAs
	  else {
	    #molecule type
	    my $rna_value = ($acc2molecule{$idF2});
	    print OUT_ACE "Properties RNA $rna_value\n" unless ($rna_value eq "ESTormRNA");
	    if ($rna_value eq "ESTormRNA") {
	      print OUT_ACE "Properties RNA mRNA\n";
	    }
	  }
	  # method depends on molecule type?? NDB or EST_elegans && will depend on organism.
	  if ($type eq "STD") {
	    print OUT_ACE "Method \"NDB\"\n";
	  }
	  elsif ($type eq "EST") {
	    print OUT_ACE "Method \"EST_$suborganism\"\n";
	  }
	} elsif (!defined $type) {
	  print OUT_ACE "Method \"NDB\"\n";
#	  print OUT_ACE "Properties RNA mRNA\n";
	}
      }				#end of record flag loop
   }                            #close returned entry loop and on to the next new sequence
    
    # Check for features on the retrieved DNA.
    &feature_finder ($seq,$idF2);
  }
  #close for each entries loop and on to the next entry.
  #close NEW_SEQUENCE file handle.
  close (NEW_SEQUENCE);
  #Close files
  close (OUT_ACE);
  close (OUT_DNA) if $dna;
  close (OUT_LONG) unless (defined $nolongtext);
  $log->write_to("\n\nOutput Files:\n");
  $log->write_to("Ace file for $suborganism $submol_mod => $acefile\n");
  $log->write_to ("DNA file for $suborganism $submol_mod => $dnafile\n") if ($dna);
  $log->write_to ("LongText file for $suborganism $submol_mod => $Longtextfile\n\n") unless (defined $nolongtext);
}

#####################################################################################################################
# generate_data_from_flat                                                                                           #
# subroutine to parse a file containing multiple EMBL flat files to produce .ace .dna and .longtext files for acedb.#
#####################################################################################################################

sub generate_data_from_flat {
  my $inputfile = shift;
  $output_dir = "/nfs/wormpub";
  #Output files
  $acefile = "$inputfile.ace";
  $wormbase->run_command ("rm $acefile", $log) if (-e $acefile);

  if ($dna) {
    $dnafile = "$inputfile.dna";
    $wormbase->run_command ("rm $dnafile", $log) if (-e $dnafile);
  }

  $Longtextfile = "$inputfile.txt";
  $wormbase->run_command ("rm $Longtextfile", $log) if (-e $Longtextfile); 
  
  #open output file for full entries.
  open (OUT_ACE,  ">$acefile");
  print "$acefile opened!!\n" if $debug;
  open (OUT_DNA,  ">$dnafile") if defined($dna);
  open (OUT_LONG, ">$Longtextfile") unless (defined $nolongtext);

  my (@description,$seq,$status,$protid,$protver,$start,$idF,$svF,$idF2,$type,$def,$sequencelength,$subspecies,$species_name);

 
  #     $new_sequence = "DQ342049"; #get 1 entry DQ342049.1
  open (NEW_SEQUENCE, "<$inputfile");
  
  while (<NEW_SEQUENCE>) {
      #my $idF2;
      #Extract ID and Sequence Version.
     
     if (/^ID\s+/) {
       $start = "1";
       @description = "";
       $seq = "";
       $start = "";
       $sequencelength = "";
     }
      
     if ((/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+.+mRNA.+(EST)/) or (/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+.+mRNA.+(STD)/)) { #mRNA or EST
	#        if (/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+/) {
	print OUT_LONG "\nLongText : \"$1\"\n" unless (defined $nolongtext);
	$idF = $1;
	$svF = $2;
	$type = $3;
      } elsif (/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+/) { #Non-EST/mRNA entries
	print OUT_LONG "\nLongText : \"$1\"\n" unless (defined $nolongtext);
	$idF = $1;
	$svF = $2;
      }
      # grab species info
      if (/^OS\s+(.+)/) {
	$subspecies = $1;
	if (/^OS\s+Caenorhabditis\s+(\S+)/) {$species_name = $1;}
	else {$species_name = "unknown";}
      }
     my $DEline1;
     if (/^DE/) {
       if (/^DE.+(yk\d+\S+\d+)\s+:\s+(\d)/) { #DE   Caenorhabditis elegans cDNA clone yk181f5 : 5' end, single read.
	 $idF2 = "$1".".$2";
	 $DEline1 = "1";
       } elsif (/^DE\s+(\S+)\s+Chris Martin.+/) { #DE   CEL01A1S1 Chris Martin sorted cDNA library Caenorhabditis elegans cDNA
	 $idF2 = $1;
	 $DEline1 = "1";
       } elsif (/DE\s+(OST\S+)\s+/) {
	 $idF2 =  $1;
	 $DEline1 = "1";
       } elsif (/DE\s+.+clone\s+(\S+)\s+similar/) { #DE   cDNA clone CEESE35 similar to ATP synthase lipid binding protein P1, mRNA
	 $idF2 =  $1;
	 $DEline1 = "1";
	} else {
	  $idF2 = $idF;		# if ($DEline1 eq "0");
	}
     }
     # grab various details out of EMBL entry
     if (/^DE\s+(.+)/) {
       $def = $1;
       # remove any offending '>' from def line. This is required by transcriptmasker.pl
       $def =~ s/\>//g;
       push(@description,"$def");
     }
     # Get protein ID info.
     if (/^FT\s+\/protein_id=\"(\S+)\.(\d+)\"/) {
       $protid=$1; $protver=$2;
     }
     # Start of DNA sequence in entry.
     if (/^SQ\s+\Sequence\s+(\d+)/) {
       $sequencelength = $1;
       print "$sequencelength\n" if ($verbose);
       print OUT_DNA "\n>$idF2\n" if $dna;
	print OUT_ACE "\nDNA \: \"$idF2\"\n";
	$status = 1;		#set status to 1 for identifying start of sequence.
	print "$status\n" if ($verbose);
      }
      print OUT_LONG "$_" unless (defined $nolongtext); #print out the embl entry line unmodified to the LongText file
      if (/^\s/) {
	s/\s+//g;
	s/\d+//g;
	s/[^acgtn]//ig;
	print OUT_DNA "$_\n" if (($status eq "1") && $dna);
	chomp;
	$seq .= "$_" if ($status eq "1");
	print OUT_ACE "$_\n" if ($status eq "1");
      }
      if (/^\/\//) {		#end of embl entry
	$status = 0;
	$DEline1 = 0;
	$start = 0
      }
      if ($start eq "0") {
	print OUT_LONG "***LongTextEnd***\n" unless (defined $nolongtext);
	print "$status\n" if ($verbose);
	print OUT_ACE "\nSequence : \"$idF2\"\n";
	#	print OUT_ACE "DNA $idF2 $sequencelength\n";
	print OUT_ACE "DNA $idF2\n";
	print OUT_ACE "Database EMBL NDB_AC $idF\n";
	print OUT_ACE "Database EMBL NDB_ID $idF\n";
	print OUT_ACE "Database EMBL NDB_SV $idF\.$svF\n";
	print OUT_ACE "DB_annotation EMBL $idF\n";
	print OUT_ACE "Protein_id $idF $protid $protver\n" if (defined$protid && defined$protver);
	print OUT_ACE "Species \"$subspecies\"\n";
	print OUT_ACE "Title \"@description\"\n";
	
	# Properties depend on the molecule subdivision of embl that the sequence was fetched from.
	
	#ESTs
	if (defined$Type) {
	  if ($Type eq "EST") {
	    print OUT_ACE "Properties cDNA cDNA_EST\n";
	  }
	  #mRNAs
	  else {
	    print OUT_ACE "Properties RNA mRNA\n";
	  }
	  # method depends on molecule type?? NDB or EST_elegans && will depend on organism.
	  if ($Type eq "mRNA") {
	    print OUT_ACE "Method \"NDB\"\n";
	  }
	  elsif ($Type eq "EST") {
	    print OUT_ACE "Method \"EST_${species_name}\"\n";
	  }
	} 
	elsif (!defined $Type) {
	  print OUT_ACE "Method \"NDB\"\n";
	}
	&feature_finder ($seq,$idF2),
      }
   }				#end of record flag loop

  close (NEW_SEQUENCE);
  close (OUT_ACE);
  close (OUT_DNA) if $dna;
  close (OUT_LONG) unless (defined $nolongtext);
  $log->write_to("\n\nOutput Files:\n");
  $log->write_to("Ace file => $acefile\n");
  $log->write_to ("DNA file => $dnafile\n") if ($dna);
  $log->write_to ("LongText file => $Longtextfile\n\n") unless (defined $nolongtext);
}


#############################################################################
# feature_finder                                                            #
# subroutine used to identify features contained within the transcript data #
# and to identify poor quality not clipped.                                 #
#############################################################################

    sub feature_finder {
      my $subseq = shift;
      my $subidF2 = shift;
      my $feature=Features::annot($subseq,$subidF2);
      if ($feature) {
	chomp $feature;
	print OUT_ACE "\n",$feature;
      }
    }

#########################################
# Load New Data Into Primary            #
# Loads new data into primary database. #
#########################################

sub load_data {
  my $sub_sourceDB = shift;
  my $submol = shift;
  my $suborganism = shift;
  my $submol_mod;
  if ($submol eq "genomic RNA") {$submol_mod = "genomic_RNA";}
  elsif ($submol eq "other RNA") {$submol_mod = "other_RNA";}
  elsif ($submol eq "unassigned RNA") {$submol_mod = "unassigned_RNA";}
  elsif ($submol eq "transcribed RNA") {$submol_mod = "transcribed_RNA";}
  else {$submol_mod = $submol;}

  $acefile = "$output_dir/new_${suborganism}_$submol_mod.ace";
  $wormbase->load_to_database($sub_sourceDB, $acefile, "EMBL_sequence_fetch.pl", $log) if (-e $acefile);
  $log->write_to("Loading $acefile into $sub_sourceDB\n\n") if (-e $acefile);

  unless (defined $nolongtext) {
    $Longtextfile = "$output_dir/new_${suborganism}_${submol_mod}_longtext.txt";
    $wormbase->load_to_database($sub_sourceDB, $Longtextfile, "EMBL_sequence_fetch.pl", $log) if (-e $Longtextfile);
    $log->write_to("Loading $Longtextfile into $sub_sourceDB\n\n") if (-e $Longtextfile);
  }
}
################################################################################
# dump data for BLATing                                                        #
# OBSOLETE SUBROUTINE - BLAT DATA IS NOW ONLY DUMPED WHEN A BUILD IS INITIATED #
################################################################################

sub dump_BLAT_data {
  my $dbdir = shift;
  my $subspecies = shift;
  my $EST_dir = $wormbase->basedir."_DATA/cDNA/$subspecies";
  $log->write_to("Dumping $subspecies from $dbdir\n\n");

  # Remove stale data if it exists on disk.
  my @types = ('mRNA','ncRNA','EST','OST','tc1','RST');
  foreach my $type (@types) {
    $wormbase->run_command ("rm $EST_dir/${type}", $log) if (-e $EST_dir."/${type}");
    $log->write_to("Removed $EST_dir/${type}\n\n")  if (-e $EST_dir."/${type}" && $debug);
  }

  my $command=<<END;
query find Sequence where method = NDB & RNA AND NEXT = mRNA & !Ignore\n
Dna -mismatch $EST_dir/mRNA\n
clear\n
query find Sequence where method = NDB & RNA AND NEXT != mRNA & !Ignore\n
Dna -mismatch $EST_dir/ncRNA\n
clear\n
query find Sequence where method = EST_$subspecies & !OST* & !Ignore & !RST*\n
Dna -mismatch $EST_dir/EST\n
clear\n
query find Sequence where method = EST_$subspecies & OST* & !Ignore\n
Dna -mismatch $EST_dir/OST\n
clear\n
query find Sequence where method = EST_$subspecies & RST* & !Ignore\n
Dna -mismatch $EST_dir/RST\n
clear\n
query find Sequence TC*\n
Dna -mismatch $EST_dir/tc1\n
clear\n
query find Sequence where method = "*EST*"\n
Write $miscdir/${subspecies}_EST_data_tmp.ace\n
clear\n
query find Sequence where method = "*NDB*"\n
Write $miscdir/${subspecies}_mRNA_data_tmp.ace\n
clear\n
quit\n
END
  print $command if ($debug);
  open (DB, "| $tace $dbdir") || die "Couldn't open $dbdir\n";
  print DB $command;
  close DB;
  $log->write_to("Finished $subspecies\n\n");
}


__END__

=pod

=head2 NAME - EMBL_Sequencefetch_species.pl

=head1 USAGE

=over 4

=item EMBL_Sequencefetch_species.pl  [-options]

=back

This script gets sequence data for all species (currently hacked to store the taxid within the script)

EMBL_Sequencefetch_species.pl

MANDATORY arguments:

=over 4

=item None at present.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4

=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.

=back

=over 4

=item -test, Test mode.

=back

=over 4

=item -species, only performs syncronisation for given species

=back

=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Paul Davis (pad@sanger.ac.uk)

=back

=cut
