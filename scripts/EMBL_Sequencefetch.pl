#!/software/bin/perl -w
#
# EMBL_Sequencefetch_species.pl
#
# Usage : EMBL_Sequencefetch.pl [-options]
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2015-07-02 10:49:33 $

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
BEGIN { $ENV{http_proxy}="http://wwwcache.sanger.ac.uk:3128";}
use warnings;
use LWP::Simple;

######################################
# variables and command-line options #
######################################

my ($help,$debug,$test,$verbose,$store,$wormbase,$organism,$output,$nolongtext,$dna,$database,$species,$dump_only,$input,$repull,$noload, $fasta, $source);

GetOptions ("help"       => \$help,       #
            "debug=s"    => \$debug,      # debug option, turns on more printing and only email specified user.
            "test"       => \$test,       # only genomic RNAs will be used for a quick test.
            "verbose"    => \$verbose,    # additional printing
            "store:s"    => \$store,      # wormbase storable object
	    "species=s"  => \$species,    # Specify a species if you have data from just a single species.
	    "dna"        => \$dna,        # Create a seperate DNA file.
	    "database=s" => \$database,   # Database are you downloading sequence data for.
	    "dump_only"  => \$dump_only,  # Only dumps the Transcript data from the primary databases.
	    "nolongtext" => \$nolongtext, # Don't dump longtext
	    "input=s"    => \$input,      # EMBL flat file to be parse if doing manually.
	    "fasta=s"    => \$fasta,      # Fasta flat file to parse if doing manually (e.g. a Trinity contig file of RNASeq reads).
	    "source=s"   => \$source,     # Specify the source that that the fasta file comes from. E.G. 'trinity', ....
	    "repull"     => \$repull,     # overrides the sequence data stored locally and re-fetches.
	    "noload"     => \$noload,     # Causes files to be generated but not loaded.
	    "output:s"   => \$output,     # Allows the user to specify the directory you want to save the data in.
	   );


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} 
else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
			     -species  => $species,
			   );
}

#################################
# Set up some useful stuff      #
#################################
# establish log file.
my $log = Log_files->make_build_log($wormbase);
my $tace = $wormbase->tace; # TACE PATH
# output dirs
my $output_dir;

if (defined $output) {
  $output_dir = $output;
} else {
  $output_dir = $wormbase->database('camace')."/EMBL_sequence_info";
}
$wormbase->run_command ("mkdir $output_dir", $log) if (!-e $output_dir);

# Set up the user agent and proxy data for the LWP sequence fetcher.
# Using the RobotUA so that a 1 second delay can be implemented (nicer to the web server if there are 1000s of sequences to fetch.)

my $ua = LWP::UserAgent->new;
print "UserAgent = $ua\n" if ($verbose);


#$ua->delay(1/60); # 1 second delay between requests
$ua->timeout(10);
$ua->env_proxy;

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

# http://www.ebi.ac.uk/ena/data/view/Taxon:6239&portal=sequence_coding&offset=1&length=1000&limit=1000&display=txt
# http://www.ebi.ac.uk/ena/data/view/Taxon:6239&portal=sequence&dataclass=EST&display=txt&header=true
# http://www.ebi.ac.uk/ena/data/warehouse/search?query="tax_eq(6239) AND mol_type="mRNA""&domain=sequence

my %species2taxonid = (
		       'elegans'           => "6239",
		       'briggsae'          => "6238",
		       'remanei'           => "31234",
		       'japonica'          => "281687",
		       'brenneri'          => "135651",
		       'pristionchus'      => "54126",
		       'brugia'            => "6279",
		       'cangaria'          => "96668",
		       'sratti'            => "34506",
                       'tmuris'            => "70415",
		       'ovolvulus'         => "6282",
		      );
my %taxon2species = (
		     '6239' => "elegans",
		     '6238' => "briggsae",
		     '31234' => "remanei",
		     '281687' => "japonica",
		     '135651' => "brenneri",
		     '54126' =>  "pristionchus",
		     '6279' => "brugia",
		     '96668' => "cangaria",
		     '34506' => "sratti",
   		     '70415' => "tmuris",
		     '6282' => "ovolvulus",
		    );

@molecules = ("EST","mRNA","STD");



print "\n---------------------------------------------------------------------------------------------------------
Retrieving EMBL sequences as a full update.
---------------------------------------------------------------------------------------------------------\n";


##########################
# MAIN BODY OF SCRIPT
##########################

#######################################################################################
# Default option                                                                      #
# Update sequence data for all WormBase species (depending upon options) using mfetch #
#######################################################################################
#Counters
my $retrieved = 0;
my $retrievedseen = 0;

####################################
# Optionally just parse input file #
####################################

if (defined $fasta) {
  &generate_data_from_fasta($fasta);


} elsif (defined $input) {
  print "You are generating Sequence data and Features from a supplied EMBL flat file...is this correct?\n\n";
#  if (!defined $Type) {
#    print "ERROR: You need to specify a data Type (eg. EST mRNA) using the command line option -Type\n\nDo you want to assign one now and continue? (y or n)\n";
#    my $answer=<STDIN>;
#    if ($answer eq "y\n") {
#      print "Please input a data Type ('EST' or 'mRNA' or '?'):";
#      $Type=<STDIN>;
#    }
#    if ($answer eq "n\n") {
#      die "Failed to process $input as a data type was not specified\n To avoid this issue please use -Type in conjunction with -input on command line\n"
#    }
#  }
  &generate_data_from_flat($input) ;



#########################################
# Get the data from the ENA's warehouse #
#########################################

} else { 
  my @organism;
  if (!defined($species)) {
    @organism = (keys %species2taxonid);
  }
  elsif ($species) {
    push (@organism, $species);
  }

  my @entries;
  foreach my $organism(@organism) {
    $log->write_to("============================================\nProcessing: $organism\n============================================\n");
    $species = $organism;
    my $sourceDB;
    if ($database){
      $sourceDB = $database;
    }
    elsif ($species eq "elegans") {
      $sourceDB = $wormbase->database('camace');
    }
    else {$sourceDB = $wormbase->database($species);}
    
    my $taxid = $species2taxonid{$species};
    my $full_species_name = $wormbase->full_name;

#Fetch existing sequence data from our canonical database
# Fetch sequence data from primary database.

    $log->write_to("Fetching sequence data from $sourceDB:\n");
    &fetch_database_info ($sourceDB) if ((-e $sourceDB) && !((defined $dump_only) || (defined $repull)));

    foreach $molecule (@molecules) {
#Fetch EST and mRNA data from ENA
      $log->write_to("Fetching $molecule Sequence data from EMBL: $species \n");
      my $work = &get_embl_data ($molecule, $taxid, $full_species_name, $species) if (!defined $dump_only);
      

# Process the new data.
      if ($work eq "1") {
	$log->write_to("Fetching NEW EMBL data: $species - $molecule\n");
	$log->write_to("Process $molecule Sequence data from EMBL: $species \n");
	&generate_data_from_flat ("$output_dir/${species}_${taxid}_${molecule}.txt", $molecule, $taxid, $full_species_name, $species);  
#Load data back in to the canonical databases if required
	if ($noload) {
	  $log->write_to("You need to look at the file $output_dir/new_${species}_${molecule}.ace as it hasn't been loaded\n");
	}
	else {
# Load the data back to the primary database.
	  $log->write_to("Loading $species data into primary database\n");
	  &load_data ($sourceDB, $molecule, $species) if (!defined $dump_only);
	}
      }
# if there isn't any, just print this line in the log.
      else {
	$log->write_to("No EMBL data:$species - $molecule\n");
      }
      $log->write_to("============================================\nProcessed: $molecule\nRetrieved: $retrieved\nSeen:$retrievedseen\n============================================\n\n");
    }
    $log->write_to("============================================\nFinished Processing : $organism\n============================================\n");
  }
}



$log->write_to ("WORK DONE\n");
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
    my $featstatus;
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
	$featstatus = "yes";
      } 
      else {
	$featstatus = "No";
      }

      $acc_sv2sequence{$sequenceacc} = $acc_sv;
      $method2sequence{$sequenceacc} = $method;
      $feature2seq{$sequenceacc} = $featstatus;
    }
  }
  close TACE;
  $log->write_to("DB: Got it.\n\n");
}

################################################################################################
# get_embl_data                                                                                #
# This subroutine retrieves a list of accessions and compares to what is in the acedb database #
# New data is the retrieved and stored for processing.                                         #
################################################################################################

sub get_embl_data {
  my $molecule = shift;
  my $taxid = shift;
  my $speciesfn = shift;
  my $species= shift;
  my $entries = shift;

  my $url;
  my $name;
  my $work;

  if ($molecule =~ /EST/) {
    $url = "http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22tax_eq%28$taxid%29%20AND%20dataclass=%22EST%22%22&result=sequence_release&display=fasta";
    if ($debug) {print "$url\n";}
    $name="$output_dir/${species}_${taxid}_EST.txt";
  }
  
  elsif ($molecule =~ /mRNA/) {
    $url = "http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22tax_eq%28$taxid%29%20AND%20dataclass!=%22STS%22%20and%20dataclass!=%22PAT%22%20and%20dataclass=%22STD%22%20and%20mol_type=%22mRNA%22%22&result=sequence_release&display=fasta";
    if ($debug) {print "$url\n";}
    $name="$output_dir/${species}_${taxid}_mRNA.txt";
  }
  elsif ($molecule =~ /STD/) {
    $url = "http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22tax_eq%28$taxid%29%20AND%20dataclass!=%22STS%22%20and%20dataclass!=%22PAT%22%20and%20dataclass=%22STD%22%20and%20mol_type=%22transcribed RNA%22%22&result=sequence_release&display=fasta";
    if ($debug) {print "$url\n";}
    $name="$output_dir/${species}_${taxid}_STD.txt";
  }
  else {
    $log->write_to("UNKNOWM molecule type\n");
    die;
  }
  
  getstore( $url, $name);
  open(FIN, "< $name ") || die "Can't open $name file: $!\n";
  my $line;
  while(<FIN>) {
    if ($_ =~ /^ID/) {
      $log->write_to("EMBL: Got it.\n");
      print "SAMPLE ENTRY $_\n" if ($verbose);
      $work = "1";
    }
    else {
      $log->write_to("No EST data retrieved\n");
      $work = "0";
    }
    last;
  }
  close(FIN);
  $log->write_to("EMBL: Finished getting $species $molecule data.\n");
  return $work;
}

#####################################################################################################################
# generate_data_from_fasta                                                                                          #
# it is assumed that we are reading data for $species                                                               #
# the molecule type and other values are set by the parameter -source e.g. 'trinity'
# no longtext file is produced                                                             #
# subroutine to parse a file containing multiple fasta entries to produce .ace .dna and .longtext files for acedb.  #
#####################################################################################################################

sub generate_data_from_fasta {
  my $inputfile = shift;
  
  my $ID;
  my $taxid = $species2taxonid{$species};
  my $full_species_name = $wormbase->full_name;
  
  
  # Output files
  $acefile = "${output_dir}/$inputfile.ace";
  if ($dna) {
    $dnafile = "${output_dir}/$inputfile.dna";
  }
  
  # Remove stale data.
  $wormbase->run_command ("rm $dnafile", $log) if ($dna && -e $dnafile);
  $wormbase->run_command ("rm $acefile", $log) if (-e $acefile);
  
  #open output file for full entries.
  open (OUT_ACE,  ">$acefile");
  print "$acefile opened.\n" if $debug;
  open (OUT_DNA,  ">$dnafile") if defined($dna);
  
  my (@description, $seq, $status, $protid, $protver, $start, $idF, $svF, $idF2, $type, $def, $sequencelength, $subspecies, $species_name);
  
# Trinity entries look like:
#>SRX004867_TR2|c0_g1_i1 len=238 path=[216:0-237] [-1, 216, -2]
#TTTTTTTTATGAAAGCTGCTTTTATATAGTTTACACCTTCTATGGTTTTCTGTTGCTTTT
#TGAACCAAACACATGAGCCACCATTGAATTCATTACAGTAGTTGATTCCAATTCCGTCGA
#AATCTACCGAAACGCTGTTATCCCATAATCTAATTGATATTTTGTTGAAGTGAGTTCCGA
#GTTTATCAATAGCGGTGCGTAAACTTCTACACACCTTTCGGCAGGTTAGTATTTCCAC

  open (NEW_SEQUENCE, "<$inputfile");
  my $done = 0;
  while (<NEW_SEQUENCE>) {
    if (/^>(\S+)/) {
      if (defined $idF) { #  print out the DNA sequence of the previous entry
	$seq = &feature_finder ($seq, $idF, 1); # want to check both orientations and return rev-comp if it appears to be reversed.
	output_seq($seq, $idF, $ID, $full_species_name, $source);
	$seq = '';
      }
      
      $idF = $1;
      if ($source eq 'trinity') {
	$ID = $idF; 
      }
      $retrieved++;
      $seq = "";
      
      # Is this entry already in the database?
      if (defined $acc_sv2sequence{$idF}) {
	print "$idF is already in the database. Not loaded again.\n";
	$idF = undef; # don't output the sequence
	next;
      }
      
      print OUT_DNA ">$idF\n" if $dna;
      
    } else { # sequence line
      chomp;
      $seq .= $_;       
    }
  }
  if (defined $idF) { 
    $seq = &feature_finder ($seq, $idF, 1); # want to check both orientations and return rev-comp if it appears to be reversed.
    output_seq($seq, $idF, $ID, $full_species_name, $source); #  print out the DNA sequence of the last entry
  }
  
  close (NEW_SEQUENCE);
  close (OUT_ACE);
  close (OUT_DNA) if $dna;
  $log->write_to("\n\nOutput Files:\n");
  $log->write_to("Ace file => $acefile\n");
  $log->write_to ("DNA file => $dnafile\n") if ($dna);
}

#####################################################################################################################
# output_seq
# print details of an entry from the fasta file to the ace files.
#####################################################################################################################

sub output_seq {
  my ($seq, $idF, $ID, $full_species_name, $source) = @_;
  
  print OUT_ACE "\nDNA \: \"$idF\"\n";
  print OUT_DNA "$seq\n" if ($dna);
  print OUT_ACE "$seq\n";
  
  print OUT_ACE "\nSequence : \"$idF\"\n";
  print OUT_ACE "DNA $idF\n";
  print OUT_ACE "Database $source AC $ID\n";
  print OUT_ACE "Species \"$full_species_name\"\n";
  
  # Properties depend on the molecule subdivision of embl that the sequence was fetched from.
  
  if ($source eq 'trinity') {
    print OUT_ACE "Properties cDNA cDNA_EST\n";
    print OUT_ACE "Method \"RNASeq_${source}\"\n";
  } 
  else {
    die "Can't determine the type from the ID line for\n";
  }
}
  
  
  
#####################################################################################################################
# generate_data_from_flat this is now the main processing sub.                                                      #
# subroutine to parse a file containing multiple EMBL flat files to produce .ace .dna and .longtext files for acedb.#
#####################################################################################################################

 sub generate_data_from_flat {

  my $inputfile = shift;
  my $molecule = shift; 
  my $taxid = shift; 
  my $full_species_name = shift; 
  my $suborganism = shift;
  my $seen = "no";
    #Output files
  if (defined $molecule) {
    $acefile = "$output_dir/new_${suborganism}_$molecule.ace";
    if ($dna) {
      $dnafile = "$output_dir/new_${suborganism}_$molecule.dna";
    }
    $Longtextfile = "$output_dir/new_${suborganism}_$molecule.txt";
  }
  else {
    $output_dir = "/nfs/wormpub/";
    $acefile = "${output_dir}/$inputfile.ace";
    if ($dna) {
      $dnafile = "${output_dir}/$inputfile.dna";
    }
    $Longtextfile = "${output_dir}/$inputfile.txt";
  }

  # Remove stale data.
  if ($dna) {
    $wormbase->run_command ("rm $dnafile", $log) if (-e $dnafile);
  }
  $wormbase->run_command ("rm $acefile", $log) if (-e $acefile);
  $wormbase->run_command ("rm $Longtextfile", $log) if (-e $Longtextfile); 
  
  #open output file for full entries.
  open (OUT_ACE,  ">$acefile");
  print "$acefile opened.\n" if $debug;
  open (OUT_DNA,  ">$dnafile") if defined($dna);
  open (OUT_LONG, ">$Longtextfile") unless (defined $nolongtext);
  
  my (@description,$seq,$status,$protid,$protver,$start,$idF,$svF,$idF2,$type,$def,$sequencelength,$subspecies,$species_name);
  
 
  #     $new_sequence = "DQ342049"; #get 1 entry DQ342049.1
  open (NEW_SEQUENCE, "<$inputfile");
 
  while (<NEW_SEQUENCE>) {
      #my $idF2;
      #Extract ID and Sequence Version.
     
     if (/^ID\s+/) {
       $retrieved++;
       $start = "1";
       @description = "";
       $seq = "";
       $start = "";
       $sequencelength = "";
     }
     unless (/^ID/){
       next unless ($seen eq "no");
     }

# EST
# ID   AA007700; SV 1; linear; mRNA; EST; INV; 115 BP.
# mRNA
# ID   AB016491; SV 1; linear; mRNA; STD; INV; 1020 BP.
# ID   AF273797; SV 1; linear; transcribed RNA; STD; INV; 2438 BP.
     if ((/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+linear\;\s+mRNA\;\s+(EST)/) or 
	 (/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+linear\;\s+mRNA\;\s+(GSS)/) or 
	 (/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+linear\;\s+mRNA\;\s+(STD)/) or
	 (/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+linear\;\s+transcribed RNA\;\s+(STD)/)
	) { # EST or mRNA
       $seen = "no";
       $idF = $1;
       $svF = $2;
       $type = $3;


       #Is this entry already in the database?
       if (defined $acc_sv2sequence{$idF}) { # If yes
	 print "PRIMARY:$acc_sv2sequence{$idF} EMBL:$svF\n" if ($verbose);
	 if ($acc_sv2sequence{$idF} == $svF) { # If the SV is the same don't do anything
	   $log->write_to("$idF Sequence Versions match PRIMARY:$acc_sv2sequence{$idF} EMBL:$svF\n") if $verbose;
	   $seen = "yes";
	   print OUT_ACE "\/\/ $idF seen before\n";
	   $retrievedseen++;
	 }
	 else {              # If the SV is different update it!!
	   $log->write_to("$idF - Sequence Versions don\'t match PRIMARY:$acc_sv2sequence{$idF} EMBL:$svF\n");
	   $log->write_to("Updated entry $idF fetch it!\n");
	 }
	 if ($feature2seq{$idF} eq "Yes") { # be careful of feature_data.
	   print "Warning $idF has associated Feature_data....this will need updating\n";
	   $log->write_to("Warning $idF has associated Feature_data....this will need updating\n");
	 }
       }
       

       print OUT_LONG "\nLongText : \"$1\"\n" unless (defined $nolongtext);
       
       if ($type eq 'GSS') {
	 $type = 'EST'
       }
     } 
     elsif (/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+/) { #Non-EST/mRNA entries
       print OUT_LONG "\nLongText : \"$1\"\n" unless (defined $nolongtext);
       $idF = $1;
       $svF = $2;
       $type = '';
     }
     # grab species info
     if (/^OS\s+(.+)/) {
       $subspecies = $1;
       if (/^OS\s+Caenorhabditis\s+(\S+)/) {
	 $species_name = $1;
       } elsif ($subspecies eq 'Brugia malayi') {
	 $species_name = 'brugia';
       } 
       elsif ($species) {
	 $species_name = $species;
       }
       else {$species_name = "unknown";}
     }
     my $DEline1;
     if (/^DE/) {
       if (/^DE.+(yk\d+\S+\d+)\s+:\s+(\d)/) { #DE   Caenorhabditis elegans cDNA clone yk181f5 : 5' end, single read.
	 $idF2 = "$1".".$2";
	 $DEline1 = "1";
       } 
       elsif (/^DE\s+(\S+)\s+Chris Martin.+/) { #DE   CEL01A1S1 Chris Martin sorted cDNA library Caenorhabditis elegans cDNA
	 $idF2 = $1;
	 $DEline1 = "1";
       } 
       elsif (/DE\s+(OST\S+)\s+/) {
	 $idF2 =  $1;
	 $DEline1 = "1";
       }
       elsif (/DE\s+.+clone\s+(\S+)\s+similar/) { #DE   cDNA clone CEESE35 similar to ATP synthase lipid binding protein P1, mRNA
	 $idF2 =  $1;
	 $DEline1 = "1";
       } 
       else {
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
       
       if ($type eq 'EST') { # EST
	 print OUT_ACE "Properties cDNA cDNA_EST\n";
	 print OUT_ACE "Method \"EST_${species_name}\"\n";
       } 
       elsif ($type eq 'STD') { # mRNA
	 print OUT_ACE "Properties RNA mRNA\n";
	 print OUT_ACE "Method \"NDB\"\n";
       } 
       else {
	 die "can't determine the type from the ID line for $idF2\n";
       }
       
       &feature_finder ($seq,$idF2);
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
  my ($subseq1, $subidF, $test_both_senses) = @_;

  $subseq1=~tr/ACGT/acgt/; #lowercase sequence

  my $feature1=Features::annot($subseq1, $subidF);
  if (!defined $test_both_senses) {
    chomp $feature1;
    print OUT_ACE "\n",$feature1;
    return $subseq1;
  }

  # we want to test both senses and use the one that appears to be the correct one, so get the other sense
  my $subseq2 = $wormbase->DNA_string_reverse($subseq1);
  my $feature2=Features::annot($subseq2, $subidF);

  if (!$feature1 && !$feature2) {
    return $subseq1;
  }

  # if only the forward sense has features or if it has a TSL
  if ($feature1 && !$feature2 || $feature1 =~ /\:TSL\"/) {
    chomp $feature1;
    print OUT_ACE "\n",$feature1;
    return $subseq1;
  } 

  # if only the reverse sense has features or if it has a TSL
  if (!$feature1 && $feature2 || $feature2 =~ /\:TSL\"/) { 
    chomp $feature2;
    print OUT_ACE "\n",$feature2;
    return $subseq2;
  }
  
  # use the sense that has a longer description, but we don't like to have to mess about masking out low complexity regions.
  if (length $feature1 >= length $feature2 || $feature2 =~ /low-complexity/) {
    chomp $feature1;
    print OUT_ACE "\n",$feature1;
    return $subseq1;
  } else {
    chomp $feature2;
    print OUT_ACE "\n",$feature2;
    return $subseq2;
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
