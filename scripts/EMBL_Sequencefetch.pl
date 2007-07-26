#!/software/bin/perl -w
#
# EMBL_Sequencefetch_species.pl
#
# Usage : EMBL_Sequencefetch.pl [-options]
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2007-07-26 17:05:35 $

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

my ($help, $debug, $test, $verbose, $store, $wormbase, $version, $organism, $output, $longtext, $dna, $database, $inc, $species, $dump_only);

GetOptions ("help"       => \$help, #
            "debug=s"    => \$debug, # debug option, turns on more printing and only email specified user.
            "test"       => \$test, # only genomic RNAs will be used for a quick test.
            "verbose"    => \$verbose, # additional printing
            "store:s"    => \$store, # wormbase storable object
	    "organism=s" => \$organism, # Specify an organism
	    "longtext"   => \$longtext, # Create longtext objects.
	    "dna"        => \$dna, # Create a seperate DNA file.
	    "database=s" => \$database, #database are you downloading sequence data for.
	    "inc"        => \$inc, #just query emblnew.
	    "dump_only"  => \$dump_only # only dumps the Transcript data from the primary databases, not EMBL connection is made.
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
		       );

if ($test) {
  @molecules  = ("mRNA");
  print "\n---------------------------------------------------------------------------------------------------------\n
WARNING: You are running in test mode, only molecule(s) [@molecules] will be returned!!!!!!!\n
---------------------------------------------------------------------------------------------------------\n";
}
else {
  @molecules  = ("genomic RNA","mRNA","other RNA","pre-RNA","rRNA","snRNA","unassigned RNA","tRNA" );
}

# incremental mode
if ($inc) {
  print "\n---------------------------------------------------------------------------------------------------------\n
WARNING: You are retrieving EMBL sequences as an incremental update, to run a full re-sync, remove the -inc.\n
---------------------------------------------------------------------------------------------------------\n";
}
# Full mode
if (!defined($inc)) { 
print "\n---------------------------------------------------------------------------------------------------------\n
WARNING: You are retrieving EMBL sequences as a full update (This will take a long time)!\nTo run an incremental update, add -inc to the command line.\n
---------------------------------------------------------------------------------------------------------\n";
}


##########################
# MAIN BODY OF SCRIPT
##########################

my %species = ($wormbase->species_accessors);
$species{$wormbase->species} = $wormbase;

foreach my $organism(keys %species) {
  $log->write_to("============================================\nProcessing: $organism\n============================================\n");
  my $sourceDB;
  if ($organism eq "elegans") {
    $sourceDB = $wormbase->database('camace');
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
    }
    # if there isn't any, just print this line in the log.
    if (scalar@entries < 1) {
      $log->write_to("No NEW EMBL data:$organism - $molecule\n");
    }
    # Load the data back to the primary database.
    $log->write_to("Loading $organism data into primary database\n");
    &load_data ($sourceDB, $molecule, $organism) if (!defined $dump_only);
    $log->write_to("Processed: $molecule\n============================================\n");
  }
  # Dump the data back to BUILD_DATA/cDNA/organism/
  $log->write_to("Dumping BLAT sequences($organism)\n");
  &dump_BLAT_data ($sourceDB, $organism, \%species);
  $log->write_to("Processed: $organism\n============================================\n\n");
}
$log->write_to ("WORK DONE------------FINISHED\n\n");
#Close log/files and exit
$log->mail();
exit(0);




###################################################
#                 SUBROUTINES                     #
###################################################


#######################
# fetch_database_info #
#######################

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

      # add data to hash(s) sequence_accession_no is key, status/method/acc_sv is value
#      %sequencedata = (
#		       SEQUENCE => {
#				    sequence2acc_sv => '$acc_sv',
#				    sequence2method => '$method',
#				    sequence2feature => '$status',
#				   }
#		      );
      $acc_sv2sequence{$sequenceacc} = $acc_sv;
      $method2sequence{$sequenceacc} = $method;
      $feature2seq{$sequenceacc} = $status;
    }
  }
  close TACE;
  $log->write_to("DB: Got it!!\n\n");
}

#################
# get_embl_data #
#################

sub get_embl_data {
  my $entries = shift;
  my $molecule = shift;
  my $species = shift;
#  my (@entries, $count, $molname);
  if ($inc) {
    open (SEQUENCES, "$mfetch -d emblnew -i \"mol:$molecule&org:Caenorhabditis $species\" |");
    print "Incremental update selected <mfetch -d emblnew -i \"mol:$molecule&org:Caenorhabditis $species\">\n";
  }
  else {
      open (SEQUENCES, "$mfetch -d embl -i \"mol:$molecule&org:Caenorhabditis $species\" |");
      print "Default Full update selected <mfetch -d embl -i \"mol:$molecule&org:Caenorhabditis $species\">\n";
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
#    $count = (@entries);
    $log->write_to(scalar@$entries." entrie(s) in the $species $molecule subclass need to be retrieved from embl\n");
    $log->write_to("EMBL: Got it!!\n");

    ###################################################################################
    # Let the real work commence - fetch the entire entry for new sequences from embl #
    ###################################################################################
    # examples of problems/exceptions
    #D36756.1   ID   D36756;   SV 1; linear; mRNA;        EST; INV; 360  BP. yk36g11.5 exception
    #CB400075.1 ID   CB400075; SV 1; linear; mRNA;        EST; INV; 559  BP. OST exception
    #M89002.1   ID   M89002;   SV 1; linear; mRNA;        EST; INV; 396  BP. (CEL12H5) Chris Martin EST exception
    #AU115113.2 ID   AU115113; SV 2; linear; mRNA;        EST; INV; 236  BP. yk726a3.3 exception
    #CB393073.1 ID   CB393073; SV 1; linear; mRNA;        EST; INV; 331  BP. OSTR113E12_1 OST exception
    #T00346.1   ID   T00346;   SV 1; linear; mRNA;        EST; INV; 318  BP. (CEESE35) Stratagene est exception
    #AF273829.1 ID   AF273829; SV 1; linear; mRNA;        STD; INV; 1183 BP. Standard mRNA NDB
    #L25083.1   ID   L25083;   SV 1; linear; genomic DNA; STD; INV; 412  BP. Standard DNA NDB
}

################
# get_new_data #
################
# Basically queries EMBL for all $subentries passed from 1st query and writes out the associated output files.
sub get_new_data {
  my $subentries = shift;
  my $submol = shift;
  my $suborganism = shift;

  #Output files
  $acefile = "$output_dir/new_${suborganism}_$submol.ace";
  $wormbase->run_command ("rm $acefile", $log) if (-e $acefile);

  if ($dna) {
    $dnafile = "$output_dir/new_${suborganism}_$submol.dna";
    $wormbase->run_command ("rm $dnafile", $log) if (-e $dnafile);
  }
  if ($longtext) {
    $Longtextfile = "$output_dir/new_${suborganism}_${submol}_longtext.txt";
    $wormbase->run_command ("rm $Longtextfile", $log) if (-e $Longtextfile);
  }

  # Remove stale data if it exists on disk.
  $wormbase->run_command ("rm $acefile", $log) if (-e $acefile);
  $wormbase->run_command ("rm $dnafile", $log) if (-e $dnafile);
  $wormbase->run_command ("rm $Longtextfile", $log) if (-e $Longtextfile);
							
  #open output file for full entries.
  open (OUT_ACE,  ">$acefile");
  print "$acefile opened!!\n" if $debug;
  open (OUT_DNA,  ">$dnafile") if defined($dna);
  open (OUT_LONG, ">$Longtextfile") if defined($longtext);

  foreach $new_sequence (@{$subentries}) {
    my ($seq,$status,$protid,$protver,@description,$idF,$svF,$idF2,$type,$def,$sequencelength);

    #     $new_sequence = "DQ342049"; #get 1 entry DQ342049.1
    open (NEW_SEQUENCE, "$mfetch -d embl -v full $new_sequence |");
    while (<NEW_SEQUENCE>) {
      #my $idF2;
      #Extract ID and Sequence Version.
      if ((/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+.+mRNA.+(EST)/) or (/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+.+mRNA.+(STD)/)) { #mRNA or EST
	#        if (/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+/) {
	print OUT_LONG "\nLongText : \"$1\"\n" if $longtext;
	$idF = $1;
	$svF = $2;
	$type = $3;
      } elsif (/^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+/) { #Non-EST/mRNA entries
	print OUT_LONG "\nLongText : \"$1\"\n" if $longtext;
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
      print OUT_LONG "$_" if $longtext; #print out the embl entry line unmodified to the LongText file
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
	print OUT_LONG "***LongTextEnd***\n" if $longtext;
	print "$status\n" if ($verbose);
	print OUT_ACE "\nSequence : \"$idF2\"\n";
	#	print OUT_ACE "DNA $idF2 $sequencelength\n";
	print OUT_ACE "DNA $idF2\n";
	print OUT_ACE "Database EMBL NDB_AC $idF\n";
	print OUT_ACE "Database EMBL NDB_ID $idF\n";
	print OUT_ACE "Database EMBL NDB_SV $idF\.$svF\n";
	print OUT_ACE "DB_annotation EMBL $idF\n";
	print OUT_ACE "Protein_id $idF $protid $protver\n" if (defined$protid && defined$protver);
	print OUT_ACE "Species \"Caenorhabditis $suborganism\"\n";
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
    }                           #close returned entry loop and on to the next new sequence
    
    # Check for features on the retrieved DNA.
    my $feature=Features::annot($seq,$idF2);
    if ($feature) {
      chomp $feature;
      print OUT_ACE "\n",$feature;
    }
  }  #close for each entries loop and on to the next entry.
  #close NEW_SEQUENCE file handle.
  close (NEW_SEQUENCE);
  #Close files
  close (OUT_ACE);
  close (OUT_DNA) if $dna;
  close (OUT_LONG) if $longtext;
  $log->write_to("\n\nOutput Files:\n");
  $log->write_to("Ace file for $suborganism $submol => $acefile\n");
  $log->write_to ("DNA file for $suborganism $submol => $dnafile\n") if ($dna);
  $log->write_to ("LongText file for $suborganism $submol => $Longtextfile\n\n") if $longtext;
}

##############################
# Load New Data Into Primary #
##############################
#Loads new data into primary database.
sub load_data {
  my $sub_sourceDB = shift;
  my $submol = shift;
  my $suborganism = shift;

  $acefile = "$output_dir/new_${suborganism}_$submol.ace";
  $wormbase->load_to_database($sub_sourceDB, $acefile, "EMBL_sequence_fetch.pl", $log) if (-e $acefile);
  $log->write_to("Loading $acefile into $sub_sourceDB\n\n") if (-e $acefile);

  if ($longtext) {
    $Longtextfile = "$output_dir/new_${suborganism}_${submol}_longtext.txt";
    $wormbase->load_to_database($sub_sourceDB, $Longtextfile, "EMBL_sequence_fetch.pl", $log) if (-e $Longtextfile);
    $log->write_to("Loading $Longtextfile into $sub_sourceDB\n\n") if (-e $Longtextfile);
  }
}

##########################
# dump data for BLATing  #
##########################

sub dump_BLAT_data {
  my $dbdir = shift;
  my $subspecies = shift;
  my $EST_dir = $species{$subspecies}{cdna_dir};
  
  $log->write_to("Dumping $subspecies from $dbdir\n\n");
  
  # Remove stale data if it exists on disk.
  $wormbase->run_command ("rm $EST_dir/mRNA", $log) if (-e $EST_dir."/mRNA");
  $log->write_to("Removed $EST_dir/mRNA\n\n")  if (-e $EST_dir."/mRNA" && $debug);
  $wormbase->run_command ("rm $EST_dir/ncRNA", $log) if (-e $EST_dir."/ncRNA");
  $log->write_to("Removed $EST_dir/ncRNA\n\n")  if (-e $EST_dir."/ncRNA" && $debug);
  $wormbase->run_command ("rm $EST_dir/EST", $log) if (-e $EST_dir."/EST");
  $log->write_to("Removed $EST_dir/EST\n\n")  if (-e $EST_dir."/EST" && $debug);
  $wormbase->run_command ("rm $EST_dir/OST", $log) if (-e $EST_dir."/OST");
  $log->write_to("Removed $EST_dir/OST\n\n")  if (-e $EST_dir."/OST" && $debug);
  $wormbase->run_command ("rm $EST_dir/tc1", $log) if (-e $EST_dir."/tc1");
  $log->write_to("Removed $EST_dir/tc1\n\n")  if (-e $EST_dir."/tc1" && $debug);

  my $command=<<END;
query find Sequence where method = NDB & RNA AND NEXT = mRNA\n
Dna -mismatch $EST_dir/mRNA\n
clear\n
query find Sequence where method = NDB & RNA AND NEXT != mRNA\n
Dna -mismatch $EST_dir/ncRNA\n
clear\n
query find Sequence where method = EST_$subspecies & !OST*\n
Dna -mismatch $EST_dir/EST\n
clear\n
query find Sequence where method = EST_$subspecies & OST*\n
Dna -mismatch $EST_dir/OST\n
clear\n
query find Sequence TC*\n
Dna -mismatch $EST_dir/tc1\n
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

This script does...blah blah blah

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
