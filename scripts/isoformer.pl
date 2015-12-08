#!/usr/local/bin/perl5.8.0 -w
#
# isoformer.pl
# 
# by Gary Williams                         
#
# This does stuff with what is in the active zone
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2015-04-29 12:30:27 $      

# Things isoformer gets confused by or misses:
# - non-canonical spliced introns where the RNASeq intron is placed on the positive strand and so is missing from reverse-strand genes
# - retained introns
# - isoforms that start in the second (or more) exon
# - isoforms that terminate prematurely
# - two or more separate sites of TSL in the same intron cause multiple identical structuers to be created
# - existing structures where the Sequence span is not the same as the span of the exons so a new structure look unique when compared to it


use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
use Sequence_extract;
use Coords_converter;
use Modules::Isoformer;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($species, $database, $gff, $notsl);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "species:s"  => \$species,
	    "database:s" => \$database, # database being curated
	    "gff:s"      => \$gff, # optional location of the GFF file if it is not in the normal place
	    "notsl"      => \$notsl, # don't try to make TSL isoforms - for debugging purposes
	    );

# always in debug mode
$debug = $ENV{USER};


if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
			     );
}

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

if (! defined $species) {$species = $wormbase->species}

my $USER = $ENV{USER};
if (! defined $database) {
  if ($species eq 'elegans') {
    $database = "/nfs/wormpub/camace_${USER}";
  } else {
    $database = "/nfs/wormpub/${species}_curation";
  }
}

##########################
# MAIN BODY OF SCRIPT
##########################

# splices are held as array of splice structures

# splice node hash structure
# int - flag to ignore this splice
# char - seq - chromosome or contig
# int - start
# int - end
# char - sense
# int - score
# int array - list of possible child introns in the structure
# int - pos (current position in list of child nodes)

my $output = "$database/tmp/isoformer_method$$";

my $coords = Coords_converter->invoke($database, undef, $wormbase);

mkdir "$database/CHROMOSOMES", 0777; # we may need to write a new set of chromosome files if they do not yet exist
my $seq_obj = Sequence_extract->invoke($database, undef, $wormbase);

my $Iso = Isoformer->new($wormbase, $log, $database, $gff, $notsl);

# name of the method and base-name of the objects
my $CDS_name = "isoformer";
my $ncRNA_name = "non_coding_transcript_isoformer";

# load the method to display the isoformer CDS structures
&load_isoformer_method();


MAIN:
while (1) {
  
  my ($chromosome, $region_start, $region_end, $sense, $biotype, $gene);
  
  do {
    print "SAVE your session > ";
    #print "Connecting to Ace\n";
    my $userinput =  <STDIN>;
    $Iso->db_close();
    $Iso->db_connect();
    chomp ($userinput);
    #print "user input: $userinput\n";
    if (!defined $userinput || $userinput eq '') {next}   # no input
    if ($userinput eq '?' || $userinput eq 'h' || $userinput eq 'help') { # help
      print "?, h, help             : this help\n";
      print "q, quit                : quit\n";
      print "cds_name               : search for structures in the region covered by the CDS\n";
      #      print "cds_name -100          : use the region starting 100 bases before the CDS\n";
      #      print "cds_name -100 +200     : use the region starting 100 bases before and 200 bases after the CDS\n";
      print "clear, clear all       : clear all isoformer objects\n";
      print "clear isoformer_8      : clear object isoformer_8\n";
      print "clear 8 9 10           : clear object isoformer_8, isoformer_9 and isoformer_10\n";
      print "clean\n                : an alias for 'clear'";
      print "what                   : reports the isoformer object that are saved in the database\n";
      print "fix isoformer_1 AC3.3c : fix isoformer_1 to CDS/Transcript, creating it if necessary\n";
      print "pseudogene non_coding_transcript_1 AC3.3 : convert the CDS to a Pseudogene, using the isoformer structure\n";
      print "check AC3.3c           : check if the specified object's structure looks OK\n";
      print "\n";
      next
    }
    if ($userinput eq 'q' || $userinput eq 'quit') {last MAIN} # quit
    
    if ($userinput =~ /^clear\b/) {
      $Iso->clear($userinput);
      next;
    }
    if ($userinput =~ /^clean\b/) {
      $Iso->clear($userinput);
      next;
    }
    if ($userinput =~ /^what\b/) {
      $Iso->what($userinput);
      next;
    }
    if ($userinput =~ /^fix\b/) {
      &fix($userinput);
      next;
    }
    if ($userinput =~ /^pseud\b/) {
      &pseud($userinput);
      next;
    }
    if ($userinput =~ /^check\b/) {
      $Iso->check($userinput);
      next;
    }
    if ($userinput =~ /^tsl\b/) {
      $Iso->check_tsls($userinput);
      next;
    }
    
    # get the region of interest from the CDS name or clone positions
    ($chromosome, $region_start, $region_end, $sense, $biotype, $gene) = $Iso->get_active_region($userinput);
    
  } while (! defined $chromosome);
  
  
  my ($confirmed, $not_confirmed, $created, $warnings) = $Iso->make_isoforms_in_region($chromosome, $region_start, $region_end, $sense, $biotype, $gene);

  open (ISOFORM, "> $output") || die "Can't open $output to write the Method\n";
  print ISOFORM $Iso->aceout();
  $Iso->aceclear();
  close (ISOFORM);

  my $return_status = system("xremote -remote 'parse $output'");
  if ( ( $return_status >> 8 ) != 0 ) {
    die ("WARNING - X11 connection appears to be lost\n");
    #      &error_warning("WARNING", "X11 connection appears to be lost");
  } else {
    print "Loaded isoforms.\n";
  }

  
  if (@{$confirmed}) {print "\n*** The following structures were confirmed: @{$confirmed}\n";}
  if (@{$not_confirmed}) {print "\n*** THE FOLLOWING STRUCTURES WERE NOT CONFIRMED: @{$not_confirmed}\n";}
  if (@{$created}) {print "\n*** The following novel structures were created: @{$created}\n";}
  if (@{$warnings}) {print "\n@{$warnings}\n\n";}
  
  # close the ACE connection
  $Iso->db_close();
}


$log->mail();
print "Finished.\n" if ($verbose);
exit(0);



##############################################################
#
# Subroutines
#
##############################################################

# load the method to display the isoformer CDS structures

sub load_isoformer_method {

  if (!-e "$database/tmp") {mkdir "$database/tmp", 0777};

  open (OUT, "> $output") || die "Can't open $output to write the Method\n";
  $Iso->load_isoformer_method();
  print OUT $Iso->aceout();
  $Iso->aceclear();
  close(OUT);

  my $return_status = system("xremote -remote 'parse $output'");
  if ( ( $return_status >> 8 ) != 0 ) {
    die ("WARNING - X11 connection appears to be lost\n");
    #      &error_warning("WARNING", "X11 connection appears to be lost");
  } else {
    print "Loaded isoformer Method\n";
  }

}


###############################################################################
# fix selected isoformer structures to a specified name
# fix isoformer_1 AC3.3b - fix specified object to CDS/Transcript, creating it if necessary

# check that the isoformer object exists
# check whether or not the target object exists
# check whether an existing target object is of the same class
# rename the isoformer object if there is no existing target of the same class
# else transfer location and tags to the existing object
# if the target object name ends with a letter, set the Isoform tag
# check that the Gene tag is populated correctly
# warn the user if the Gene tag is not populated
# warn the user if History should be made


sub fix {
  my ($userinput) = @_;

  $Iso->fix($userinput);

  open (TARGET, ">$output") or die "cant open $output\n";
  print TARGET $Iso->aceout();
  $Iso->aceclear();
  close TARGET;

  my $return_status = system("xremote -remote 'parse $output'");
  if ( ( $return_status >> 8 ) != 0 ) {
    die("WARNING - X11 connection appears to be lost\n");
  }
}

###############################################################################
# convert a CDS to an Pseudogene using the spcified non-coding isoformer structure
# pseud isoformer_1 AC3.3 - fix specified object to CDS, converting it to a Pseudogene

# check that the isoformer object exists
# check that the isoformer object is a non-coding transcript
# check that the target object exists
# check that the target object is a CDS
# rename the isoformer object to the target CDS name
# check that the Gene tag is populated correctly
# warn the user if the Gene tag is not populated
# warn the user if History should be made
# set the Last_reviewed tag

sub pseud {
  my ($userinput) = @_;

  $Iso->pseud($userinput);

  open (TARGET, ">$output") or die "cant open $output\n";
  print TARGET $Iso->aceout();
  $Iso->aceclear();
  close TARGET;

  my $return_status = system("xremote -remote 'parse $output'");
  if ( ( $return_status >> 8 ) != 0 ) {
    die("WARNING - X11 connection appears to be lost\n");
  }
}

