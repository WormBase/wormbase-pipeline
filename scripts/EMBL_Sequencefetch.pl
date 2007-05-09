#!/nfs/disk100/wormpub/bin/perl -w
#
# EMBL_Sequencefetch.pl
#
# Usage : EMBL_Sequencefetch.pl [-options]
#
# Last edited by: $Author: pad $
# Last edited on: $Date: 2007-05-09 12:51:58 $

my $script_dir = $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'};

use strict;
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Modules::Features;

######################################
# variables and command-line options #
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $version, $organism, $output,$longtext,$dna,$database,$inc);

GetOptions ("help"       => \$help, #
            "debug=s"    => \$debug, # debug option, turns on more printing and only email specified user.
            "test"       => \$test, # only genomic RNAs will be used for a quick test.
            "verbose"    => \$verbose, # additional printing
            "store:s"    => \$store,
	    "version=s"  => \$version, # Specifies a build version for using common data etc.
	    "organism=s" => \$organism, # Specify an organism
	    "output=s"   => \$output, # Specify output directory.
	    "longtext"   => \$longtext, # Create longtext objects.
	    "dna"        => \$dna, # Create a seperate DNA file.
	    "database=s" => \$database, #database are you downloading sequence data for.
	    "inc"        => \$inc #just query emblnew.
	   );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			   );
}

# Display help if required
&usage("Help") if ($help);

# establish log file.
my $log = Log_files->make_build_log($wormbase);

#################################
# Set up some useful paths      #
#################################
my $logs_dir        = $wormbase->logs; # AUTOACE LOGS
my $sourceDB;

if ($database) {
  $sourceDB = $database;
} 
else {
  $sourceDB = $wormbase->database('camace');
}

my $output_dir;
if (defined $output) {
  $output_dir = $output;
} 
else {
  $output_dir = $sourceDB."/EMBL_sequence_info";
  $wormbase->run_command ("mkdir $output_dir", $log) if (!-e $output_dir);
}

# other paths
my $mfetch = "/usr/local/pubseq/scripts/mfetch"; #mfetch script


##########################
# MAIN BODY OF SCRIPT
##########################

# main stuff goes here
my (%acc_sv2sequence,%method2sequence,%feature2seq);
my %molecule2rnatype;
my $name;
my @molecules;
my @species;
my %acc2molecule;

#Species hash
#push(@species,"Caenorhabditis briggsae") if (!defined $organism);
push(@species,"Caenorhabditis elegans") if (!defined $organism);
push(@species,"$organism") if ($organism);


#Molecules types
#molecule type to Type tag data hash. key:molecule value:type
$molecule2rnatype{'genomic RNA'} = "ncRNA";
$molecule2rnatype{'mRNA'} = "ESTormRNA";
$molecule2rnatype{'other RNA'} = "ncRNA"; 
$molecule2rnatype{'pre-RNA'} = "ESTorRNA";
$molecule2rnatype{'rRNA'} = "rRNA";
$molecule2rnatype{'snRNA'} = "snRNA";
$molecule2rnatype{'unassigned RNA'} = "ncRNA";
$molecule2rnatype{'tRNA'} = "tRNA";

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

##########################################################################
# Retrieve datahash of Sequence entries from camace and associated data. #
##########################################################################
$log->write_to("Fetching sequence data from $sourceDB\n");
&fetch_database_info;


###########
# Species #
###########
foreach my $species (@species) {
  if ($species =~ (/(\S+)\s+(\S+)/)) {
    $name = "$1_$2";
    #print "$name\n";
  }
  $log->write_to("\n============================================\nProcessing $name\n============================================\n");

##### Need to add database stuff here, so that each database can be updated seperately.....when we have all of them here. #####

  #################
  # Molecule type #
  #################
  foreach my $molecule (@molecules) {
    my $molname;
    my $count;
    my @entries;
    if ($molecule =~ (/(\S+)\s+(\S+)/)) {
      $molname = "$1_$2";
    } else {
      $molname = $molecule;
    }
    $log->write_to("\n============================================\nProcessing: $molecule\n\n");

    #############################################
    # Limit the work to new or updated entries. #
    #############################################
    if ($inc) {
      open (SEQUENCES, "$mfetch -d emblnew -i \"mol:$molecule&org:$species\" |");
      print "Incremental update selected <mfetch -d emblnew -i \"mol:$molecule&org:$species\">\n";
    }
    else {
      open (SEQUENCES, "$mfetch -d embl -i \"mol:$molecule&org:$species\" |");
      print "Default Full update selected <mfetch -d embl -i \"mol:$molecule&org:$species\">\n";
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
	  print "CAMACE:$acc_sv2sequence{$accR} EMBL:$svR\n" if ($verbose);
	  if ($acc_sv2sequence{$accR} == $svR) { # If the SV is the same don't do anything
	    $log->write_to("$accR Sequence Versions match CAMACE:$acc_sv2sequence{$accR} EMBL:$svR\n") if $verbose;
	  } else {		# If the SV is different update it!!
	    push(@entries,"$accR\.$svR");
	    $log->write_to("$accR - Sequence Versions don\'t match CAMACE:$acc_sv2sequence{$accR} EMBL:$svR\n");
	    $log->write_to("Updated entry $accR fetch it!\n");
	  }
	  if ($feature2seq{$accR} eq "Yes") { # be careful of feature_data.
	    print "Warning $accR has associated Feature_data....this will need updating\n";
	    $log->write_to("Warning $accR has associated Feature_data....this will need updating\n");
	  }
	} else {
	  $log->write_to("New entry $accR fetch it!\n");
	  push(@entries,"$accR\.$svR");
	}
      }
    }				# while SEQUENCE
    close (SEQUENCES);
    $count = (@entries);
    $log->write_to("$count entrie(s) in the $species $molecule subclass need to be retrieved from embl\n");
    $log->write_to("\nProcessed: $molecule\n============================================\n");

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
    
    if ($count > 0) {
      $log->write_to("\nWorking to get new data.........\n\n");
      &get_embl_data (\@entries, $molname);
    }
  }				#close foreach molecule loop and on to the molecule type.
  $log->write_to ("============================================\nFinished Processing $species\n============================================\n");
}				#close foreach species loop and on to the next species.

$log->write_to ("WORK DONE------------FINISHED\n\n");
#Close log/files and exit
$log->mail();
exit(0);


###################
# ==Subroutines== #
###################

sub fetch_database_info {
  my $tace = $wormbase->tace; # TACE PATH
  my $def_dir = "$sourceDB/wquery";
  my $tablemaker_query =  "${def_dir}/SCRIPT:estfetch.def";
  my $command = "Table-maker -p $tablemaker_query\nquit\n";

  open (TACE, "echo '$command' | $tace $sourceDB |");
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
      }				#print "Sequence: $sequence has the feature $4\n\n"}
      if (defined$feature) {
	$status = "yes";
      } else {
	$status = "No";
      }
      # add data to hash(s) sequence_accession_no is key, status/method/acc_sv is value
      $acc_sv2sequence{$sequenceacc} = $acc_sv;
      $method2sequence{$sequenceacc} = $method;
      $feature2seq{$sequenceacc} = $status;
    }
  }
  close TACE;
  $log->write_to("=>\n=>\n=>\nGot it!!\n\n");
}

sub get_embl_data {
  my $new_sequence;
  my $subentries = shift;
  my $submol = shift;
  my ($acefile, $dnafile, $Longtextfile);
  #Output files
  $acefile = "$output_dir/new_${name}_$submol.ace";
  $wormbase->run_command ("rm $acefile", $log) if (-e $acefile);

  if ($dna) {
    $dnafile = "$output_dir/new_${name}_$submol.dna";
    $wormbase->run_command ("rm $dnafile", $log) if (-e $dnafile);
  }
  if ($longtext) {
    $Longtextfile = "$output_dir/new_${name}_${submol}_longtext.txt";
    $wormbase->run_command ("rm $Longtextfile", $log) if (-e $Longtextfile);
  }
  
  # Remove stale data if it exists on disk.
  #$wormbase->run_command ("rm $acefile", $log) if (-e $acefile);
  #$wormbase->run_command ("rm $dnafile", $log) if (-e $dnafile);
  #$wormbase->run_command ("rm $Longtextfile", $log) if (-e $Longtextfile);
							
  #open output file for full entries.
  open (OUT_ACE,  ">$acefile");
  open (OUT_DNA,  ">$dnafile") if defined($dna);
  open (OUT_LONG, ">$Longtextfile") if defined($longtext);

  foreach $new_sequence (@{$subentries}) {
    my $seq;
    my $status;
    my ($protid,$protver,@description,$idF,$svF,$idF2,$type,$def,$sequencelength);
    
    #     $new_sequence = "DQ342049"; #get 1 entry DQ342049.1
    open (NEW_SEQUENCE, "$mfetch -d embl -v full $new_sequence\" |");
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
	print OUT_ACE "Species \"Caenorhabditis elegans\"\n";
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
	    my $moleculetype = ($acc2molecule{$idF2});
	    my $rna_value = $moleculetype;
#	    my $rna_value = ($molecule2rnatype{$moleculetype});
	    print OUT_ACE "Properties RNA $rna_value\n" unless ($rna_value eq "ESTormRNA");
	    if ($rna_value eq "ESTormRNA") {
	      print OUT_ACE "Properties RNA mRNA\n";
	    }
	  }
	  # method depends on molecule type?? NDB or EST_elegans && will depend on organism.
	  if ($type eq "STD") {
	    print OUT_ACE "Method \"NDB\"\n";
	  }
	  if ($type eq "EST") {
	    print OUT_ACE "Method \"EST_elegans\"\n";
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
  $log->write_to("Ace file for $name $submol => $acefile\n");
  $log->write_to ("DNA file for $name $submol => $dnafile\n") if ($dna);
  $log->write_to ("LongText file for $name $submol => $Longtextfile\n\n") if $longtext;
}


sub usage {
  my $error = shift;
  
  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################


# Add perl documentation in POD format
# This should expad on your brief description above and
#  add details of any options that can be used with the program.
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - estfetch.pl

=head1 USAGE

=over 4

=item script_template.pl  [-options]

=back

This script does...blah blah blah

script_template.pl MANDATORY arguments:

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

=item -test, Test mode, run the script, but dont change anything.

=back

=over 4

=item -verbose, output lots of chatty test messages

=back

=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Keith Bradnam (krb@sanger.ac.uk)

=back

=cut
