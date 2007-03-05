#!/nfs/disk100/wormpub/bin/perl -w
use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Modules::Features;

######################################
# variables and command-line options #
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $version, $organism, $output,$longtext,$dna,);

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
my $output_dir;
my $sourceDB;
if (defined $output) {
  $output_dir = $output;
} 
else {
  $output_dir = "/nfs/team71/worm/pad/Scripts";
  #print "\n$output_dir\n";
}
# other paths
my $mfetch = "/usr/local/pubseq/scripts/mfetch"; #mfetch script
if ($test && $debug) {
  $sourceDB = "/nfs/disk100/wormpub/DATABASES/TEST_DBs/camace";
} 
else {
  $sourceDB = $wormbase->database('camace');
}

##########################
# MAIN BODY OF SCRIPT
##########################

# main stuff goes here
my (%acc_sv2sequence,%method2sequence,%feature2seq);
my %molecult2rnatype;
my $name;
my @molecules;
my @species;
my %acc2molecule;
my %molecule2type;
#Species hash
#push(@species,"Caenorhabditis briggsae") if (!defined $organism);
push(@species,"Caenorhabditis elegans") if (!defined $organism);
push(@species,"$organism") if ($organism);


#Molecules types
#molecule type to Type tag data hash.
$molecult2rnatype{'genomic RNA'} = "ncRNA";
$molecult2rnatype{'mRNA'} = "ESTormRNA";
$molecult2rnatype{'other RNA'} = "ncRNA"; 
$molecult2rnatype{'pre-RNA'} = "ESTorRNA";
$molecult2rnatype{'rRNA'} = "rRNA";
$molecult2rnatype{'snRNA'} = "snRNA";
$molecult2rnatype{'unassigned RNA'} = "ncRNA";
$molecult2rnatype{'tRNA'} = "tRNA";
 
if ($test) {
  @molecules  = ("mRNA");
  #@molecules  = ("pre-RNA","tRNA","unassigned RNA","rRNA","snRNA");
  #  @molecules  = ("unassigned RNA","snRNA","genomic RNA");
} else {
  @molecules  = ("genomic RNA","mRNA","other RNA","pre-RNA","rRNA","snRNA","unassigned RNA","tRNA" );
}

# in test mode(s)
if ($test) {
  print "\n---------------------------------------------------------------------------------------------------------\n
WARNING: You are running in test mode, only molecule(s) [@molecules] will be returned!!!!!!!\n
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
  #  my @entries;
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
    $log->write_to("\n\nProcessing: $molecule\n\n");

    #############################################
    # Limit the work to new or updated entries. #
    #############################################
    open (SEQUENCES, "$mfetch -d embl -f id -i \"mol:$molecule&org:$species\" |");
    while (<SEQUENCES>) {
      chomp;
      if (/no match/) {
	$log->write_to("No entries were retrieved for Species:$species in Molecule type:$molecule\n");
      }
      #eg. ID   AX254400; SV 1; linear; unassigned RNA; PAT; INV; 1161 BP.
      elsif (my($accR,$svR) = /^ID\s+(\S+)\;\s+\SV\s+(\d+)\;\s+/) {
	print "$accR.$svR\n" if ($verbose);

	# Store data about the sequence is datahash.
	# add to hash. Accession is key, molecule is value
	$acc2molecule{$accR} = $molecule2type{$molecule};

	#Is this entry already in the database?
	if (defined $acc_sv2sequence{$accR}) { # If yes
	  print "CAMACE:$acc_sv2sequence{$accR} EMBL:$svR\n" if ($verbose);
	  if ($acc_sv2sequence{$accR} == $svR) { # If the SV is the same don't do anything
	    $log->write_to("$accR Sequence Versions match CAMACE:$acc_sv2sequence{$accR} EMBL:$svR\n");
	  } else {		# If the SV is different update it!!
	    push(@entries,"$accR\.$svR");
	    $log->write_to("$accR - Sequence Versions don\'t match CAMACE:$acc_sv2sequence{$accR} EMBL:$svR\n");
	    $log->write_to("Updated entry $accR fetch it!\n");
	  }
	  if ($feature2seq{$accR} eq "Yes") { # be careful of feature_data.
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
    $log->write_to("\n\nProcessed: $molecule\n\n");

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
      $log->write_to("\n\nWorking to get new data.........\n\n");
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
  my $tace            = $wormbase->tace; # TACE PATH
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
  #Output files
  my $acefile = "$output_dir/new_${name}_$submol.ace";
  my $dnafile = "/nfs/team71/worm/pad/Scripts/new_${name}_$submol.dna" if $dna;
  my $Longtextfile = "$output_dir/new_${name}_${submol}_longtext.txt" if $longtext;
  # Is there stale data on disk??
  if (-e "$output_dir/new_${name}_$submol.ace") { #remove old output files if they exist!
    $wormbase->run_command ("rm $acefile", $log);
    $wormbase->run_command ("rm $dnafile", $log) if ($dna);
    $wormbase->run_command ("rm $Longtextfile", $log) if ($longtext);
  }
  #open output file for full entries.
  open (OUT_ACE,  ">$acefile");
  open (OUT_DNA,  ">$dnafile") if $dna;
  open (OUT_LONG, ">$Longtextfile") if $longtext;

  foreach $new_sequence (@{$subentries}) {
    my @seq;
    my $status;
    my ($protid,$protver,@description,$idF,$svF,$idF2,$type,$def,$sequencelength);
    
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
	push (@seq,"$_\n") if ($status eq "1");
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
	    my $rna_value = $molecult2rnatype{$moleculetype};
	    print OUT_ACE "Properties RNA $rna_value\n" unless ($rna_value eq "ESTorRNA");
	    if ($rna_value eq "ESTorRNA") {
	      print OUT_ACE "Properties RNA mRNA\n";
	    }
	  }
	  # method depends on molecule type?? NDB or EST_elegans
	  if ($type eq "STD") {
	    print OUT_ACE "Method \"NDB\"\n";
	  }
	  if ($type eq "EST") {
	    print OUT_ACE "Method \"EST_elegans\"\n";
	  }
	} elsif (!defined $type) {
	  print OUT_ACE "Method \"NDB\"\n";
	  print OUT_ACE "Properties RNA mRNA\n";
	}
	
	#	# method depends on molecule type?? NDB or EST_elegans
	#	if (defined$type) {
	#	  if ($type eq "STD") {
	#	    print OUT_ACE "Method \"NDB\"\n";
	#	  }
	#	  if ($type eq "EST") {
	#	    print OUT_ACE "Method \"EST_elegans\"\n";
	#	  }
	#	} 
	#	elsif (!defined $type) {
	#	  print OUT_ACE "Method \"NDB\"\n";
	#	}
      }				#end of record flag loop
      
    }                           #close returned entry loop and on to the next new sequence
    
    # Check for features on the retrieved DNA.
    my $seq = "@seq";
    print "\n\n\n***************\nSequence ID = $idF2\n>Sequence : \"$idF2\"\n$seq\n***************\n";
    my $feature=Features::annot($seq,$idF2);
    if ($feature) {
      chomp $feature;
      print OUT_ACE "\n",$feature;
    }
  }				#close for each entries loop and on to the next entry.
  close (NEW_SEQUENCE);
  #Close files
  close (OUT_ACE);
  close (OUT_DNA) if $dna;
  close (OUT_LONG) if $longtext;
  $log->write_to("\n\nOutput Files:\n\n");
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
# This should expand on your brief description above and
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
