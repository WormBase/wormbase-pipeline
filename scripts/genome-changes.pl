#!/usr/local/bin/perl5.8.0 -w
#
# genome-changes.pl                           
# 
# by Gary Williams                        
#
# This is a script to aid making changes to the sequence of a clone.
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2012-10-26 08:56:54 $      

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

######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($database, $new_database, $infile, $nomove, $noload, $nohomol, $ignore_assembly_tags, $species, $positionsfile, $stlace, $camace, $ignore_clone_overlaps, $dump_dna_every_time);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "database:s" => \$database,	# specify database other than camace for test purposes
	    "new_database:s" => \$new_database,	# specify the path to the copy that is made that is to be changed
	    "infile:s"     => \$infile, # file specifying the changes to be made
	    "positions:s"  => \$positionsfile, # file specifying the positions of the changes to be made - an alternative to using -infile
	    "nomove"     => \$nomove, # assume that the new database has been copied already
	    "noload"     => \$noload,	# for testing only, don't load ace files
	    "nohomol"    => \$nohomol,	# for testing only, don't change superlink homol data
	    "ignore_assembly_tags" => \$ignore_assembly_tags, # ignore assembly tags, report them instead of dieing
	    "ignore_clone_overlaps" => \$ignore_clone_overlaps, # ignore clone overlaps, report them instead of dieing
	    "species=s"  => \$species,
	    "stlace"     => \$stlace, # tell makesuperlinks.pl that the stlace superlinks should be made
	    "camace"     => \$camace, # tell makesuperlinks.pl that the camace superlinks should be made
	    "dump_dna_every_time" => \$dump_dna_every_time # test every correction by doing a DNA dump of superlinks and chromosomes - slow!
	    );

#if (! defined $database || $database eq "") {
#  die "-database not specified\n";
#}

if (! defined $new_database || $new_database eq "") {
  die "-new_database not specified\n";
}

if ($database eq $new_database) {
  die "the old and the new databases are the same\n";
}

if (!defined $infile && !defined $positionsfile) {
  die "You must define the changes to be made in one of -infile or -positions\n";
}

# always in debug mode while testing
#$test = 1;
#$debug = "gw3";

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     -organism => $species,
			     );
}

$species = $wormbase->species;

# Display help if required
&usage("Help") if ($help);

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

#################################
# Set up some useful paths      #
#################################

# Set up top level base directories (these are different if in test mode)
my $basedir         = $wormbase->basedir; # BASE DIR
my $ace_dir         = $wormbase->autoace; # AUTOACE DATABASE DIR

# some database paths
my $currentdb = $wormbase->database('current');
my %clonelab;
if ($species eq 'elegans') {
  $wormbase->FetchData('clone2centre', \%clonelab, "$currentdb/COMMON_DATA/");
}

# other paths
my $tace            = $wormbase->tace; # TACE PATH



##########################
# MAIN BODY OF SCRIPT
##########################

# main stuff goes here

# See: /nfs/WWWdev/INTWEB_docs/htdocs/Projects/C_elegans/DOCS/Sequence_update.shtml
#
#- make a copy of camace database to work on
#- input details of sequence change and check they are properly formatetd and correctly identify a genomic position
#- for each genomic change:
#  - make the changes to the DNA sequence in the database
#  - shift up gene models etc. in the database
#  - load the resulting clone ace file into the copy of camace
#  - report any objects that are not easy to shift and need manual curation
#  - run 'makesuperlinks.pl' to recreate the superlink structure from its clones
#  - check we have one of each superlink
#  - load the resulting superlink ace file into the copy of camace
#  - check on updated superlink coords
#  - dump out superlinks and check for non AGCT bases
#  - dump out curated CDSs and check for '*' characters
#  - update cosmid clone sequences on file

# 
# All clone and superlink change information and ACE commands for a change
# are stored in the hash %change.
# The ACE commands for deletions are stored under a list like:
# @{$change->{'ace-delete'}{"Sequence : \"$clone\""}}, "-D $line";
# and the ACE commands for the replacement new data are stored under a list like:
# @{$change->{'ace-add'}{"Sequence : \"$clone\""}}, "$new_line";
#
# A new ace connection is made between writing changes to the new
# database and the next read frmo it so that fresh data is read in teh
# event of making several changes to the same clone.


my $errors='';

if (!defined $database) {die "-database not defined\n"}

#
# Make a working copy of camace database before making any changes we
# could then simply restart the program if there was a problem.  We
# will read data from, and write ace files to the copy
#

if ($nomove) {
  print "-nomove specified. It is therefore assumed that the new database is already in place waiting to be updated\n";
} else {
  $wormbase->run_script( "TransferDB.pl -start $database -end $new_database -name 'copy_for_genome_changes' -all", $log) && die "The new copy of the database was not transferred correctly. No changes were made.\n";
}

# open an ACE connection to parse details for mapping to genome
print "Connecting to Ace\n";
my $ace = Ace->connect (-path => $new_database,
			-program => $tace) || die "cannot connect to database at $new_database\n";


# get the coords for the chromosome to clone coord conversion simply from the currentdb
# this conversion is done before any genome changes are made
my $coords = Coords_converter->invoke($currentdb, undef, $wormbase);


if (defined $positionsfile) {
  
  # If the lines starts with '#' or does not start with 'CHROMOSOME' or 'chr' then it is ignored
  
  # The columns are then:
  
  # Chromosome : e.g. CHROMOSOME_I - the sequence to change
  # start : e.g. 232039 - the start position or the position to insert to the right of
  # end : e.g. 232039 - the end position
  # type of change : Insertion, Deletion, SNP (i.e. substitution)
  # Base(s) to change from : or '-' if insertion
  # Base(s) to change to : or '-' if deletion
  
  # any thing else on the line will be ignored.
  
  
  my %change = &read_positions_file($positionsfile);
  my %ordered = &sort_by_position(%change);
  &make_changes(%ordered);
  
  
} else {
  
  
  
  #
  # do a quick test that all of the changes map to a unique position 
  # before starting the time-consuming changes
  #
  open (IN, "< $infile") || die "Can't open the input file $infile\n";
  
  my $errors = 0;
  my $prev_type = "";
  my $prev_clone = "";
  my $prev_removed = "";
  
  while (my $input_line = <IN>) {
    chomp $input_line;
    
    if ($input_line =~ /^\s*$/) {next;} # skip blank lines
    if ($input_line =~ /^#/) {next;}    # skip comments
    
    my ($clone, $change_type, $change_sequence) = split /\s+/, $input_line; 
    
    # check that we are not trying to make a change in the mitochondrion - we don't own this sequence
    if ($clone =~ /MtDNA/) {next;}
    
    my $change = {};		# hash ref to hold details of the next change
    $change->{'clone'} = $clone;
    $change->{'change_type'} = $change_type;
    $change->{'region'} = $change_sequence; 
    
    print "\nChecking change specification: clone $change->{'clone'} $change->{'change_type'} of $change->{'region'}\n";
    
    # if the previous type is a deletion and this is an insertion and the
    # clone name is the same and the lowercase part of the sequence is the
    # same then this is a substitution and we need not worry about the
    # sequence not matching the genomic sequence because the genomic
    # sequence will be changed when we do the insertion for real.
    my ($removed, @positions) = &find_changed_part_of_region($change_sequence);
    if (lc($prev_type) eq "deletion" &&
	lc($change_type) eq "insertion" &&
	$prev_clone eq $clone &&
	$prev_removed eq $removed) {
      print "Found a substitution - OK\n";
    } else {
      #  - find the changes for the DNA sequence
      #print "change DNA sequence\n";
      if (&change_DNA_in_db($change)) {
	print "Error found - correct this line in the input file\n";
	$errors++;
      } 
    }
    
    $prev_type = $change_type;
    $prev_clone = $clone;
    $prev_removed = $removed;
  }
  close (IN);
  if ($errors) {die "There were $errors errors in the input file\n";}
  print "\n That all looks OK.\n\n";
  
  
  
  
  #
  # now do it all for real
  #
  
  open (IN, "< $infile") || die "Can't open the input file $infile\n";
  
  # next change loop
  while (my $input_line = <IN>) {
    chomp $input_line;
    
    if ($input_line =~ /^\s*$/) {next;} # skip blank lines
    if ($input_line =~ /^#/) {next;}    # skip comments
    
    my ($clone, $change_type, $change_sequence) = split /\s+/, $input_line; 
    
    my $change = {};		# hash ref to hold details of the next change
    $change->{'clone'} = $clone;
    $change->{'change_type'} = $change_type;
    $change->{'region'} = $change_sequence; 
    
    print "\nNext change: clone $change->{'clone'} $change->{'change_type'} of $change->{'region'}\n";
    
    #  - make the changes to the DNA sequence in the database
    print "change DNA sequence\n";
    if (&change_DNA_in_db($change)) {die "ERROR in change_DNA_in_db\n";} 
    
    #  - shift up gene models etc. in the database
    #  - report any objects that are not easy to shift and need manual curation
    print "shift up gene models\n";
    if (&shift_gene_models($change)) {die "ERROR in shift_gene_models\n";}
    
    &write_ace_lines($change);
    
    
    
    
    # now run 'makesuperlinks.pl' on the last clone change    
    my $super = "/tmp/makesuperlinks$$.ace";
    my $ans = "";
    while (lc($ans) ne 'n') {
      print "run 'makesuperlinks.pl' on new database. Any errors will be mailed to you.\n";
      my $stl = '';
      if ($stlace) {$stl = '-stlace'}
      my $superlink_errors = $wormbase->run_script( "makesuperlinks.pl -db $new_database -acefile $super $stl", $log);
      
      #  - check we have one of each superlink
      my @superlinks = `grep SUPERLINK $super | sort -u  | sort -n`;
      if ($camace && scalar @superlinks != 8) {
	print "\nERROR, the following superlinks names have been found:\n@superlinks\n\n";
	$superlink_errors = 1;
      }
      
      # if there were any errors, loop to run makesuperlinks again when the user is ready
      if ($superlink_errors) {
	print "\n\nErrors have been mailed to you in the report.\nCorrect them in the database $new_database now.\nThen answer 'y' to run 'makesuperlinks' again.\n";
	while (1) {
	  print "Run makesuperlinks.pl again now? [Y/N] ";
	  $ans = <STDIN>;
	  chomp $ans;
	  if (lc($ans) eq 'n' || lc($ans) eq 'y') {last;}
	}
      } else {
	$ans = 'n';			# no errors, break from while loop
      }
    }
    
    # fnished running 'makesuperlinks.pl' - load the resulting ACE file
    if ($noload) {
      print "-noload specified. No changes to the new database will be made\n";
    } else {
      # load the resulting $super ace file into newdatabase
      print "load ace file of superlinks '$super' into new database\n";
      $wormbase->load_to_database($new_database, $super, "genome_changes_pl", $log);
      
      # close then reopen an ACE connection to parse details for mapping to genome
      print "Connecting to Ace\n";
      $ace->close;
      $ace = Ace->connect (-path => $new_database,
			   -program => $tace) || die "cannot connect to database at $new_database\n";
    }
    
    
  }  # loop back to next change
  
  close (IN);
}

# --- now all changes have been made in new database 

print "\n\n\nThe following errors were noted:\n\n";
print "$errors\n\n\n";


#  - dump out superlink sequences and check for non AGCT bases
&check_superlink_sequence() || die "ERRORS were found in the superlink sequences\n";

## create ~wormpub/analysis/cosmids/$clone/YYMMDD directories to store the clone sequences
#my $clonedates = "/tmp/makeclonedates.ace";
#open (IN, "< $infile") || die "Can't open the input file $infile\n";
#while (my $input_line = <IN>) {
#  chomp $input_line;
#  if ($input_line =~ /^\s*$/) {next;} # skip blank lines
#  if ($input_line =~ /^#/) {next;}    # skip comments
#  my ($clone, $change_type, $change_sequence) = split /\s+/, $input_line; 
#  print "Writing sequence for clone $clone\n";
#  &write_clone_sequence($clone, $clonedates);
#}
#close (IN);

#if ($noload) {
#  print "-noload specified. No changes to the new database will be made\n";
#} else {
#  # load the resulting $super ace file into newdatabase
#  print "load ace file of clone dates '$clonedates' into new database\n";
#  my $command = "pparse $clonedates\nsave\nquit\n";
#  open( WRITEDAT, "| $tace -tsuser genome_changes_pl $new_database " ) || die "Couldn't open pipe to database $new_database\n";
#  print WRITEDAT $command;
#  close(WRITEDAT);
#}

# close the ACE connection
$ace->close;


$log->mail();

exit(0);






##############################################################
#
# Subroutines
#
##############################################################

##########################################
#  my %change = &read_positions_file($positionsfile);
# for the positions file input format

sub read_positions_file {
  my %change;
  
  open (IN, "$positionsfile") || die "Can't open the file $positionsfile\n";
  while (my $line = <IN>) {
    chomp $line;
    if ($line =~ /^#/) {next;}
    if ($line !~ /CHROMOSOME\S+/ && $line !~ /chr\S/) {next;}
    my ($chrom, $start, $end, $id, $type, $frombase, $tobase) = ($line =~ /(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/);

    if ($frombase eq '-' || $frombase eq '.') {$frombase = ''}
    if ($tobase eq '-'|| $tobase eq '.') {$tobase = ''}

    $change{$chrom}{$start}{$end}{'ID'} = $id;
    $change{$chrom}{$start}{$end}{'TYPE'} = $type;
    $change{$chrom}{$start}{$end}{'FROMBASE'} = $frombase;
    $change{$chrom}{$start}{$end}{'TOBASE'} = $tobase;
  }
  
  close (IN);
  return %change;
}




##########################################
#  &sort_by_position(%change);
# for the positions file input format
# sort the changes by position in the genome and give them an ordinal number
# get the clone positions now before any changes are made to the clone sequences
sub sort_by_position {
  my (%change) = @_;
  my %ordered;
  
  my $ordinal = 1;
  foreach my $chrom (sort keys %change) {
    foreach my $start (sort { $a <=> $b } keys %{$change{$chrom}}) { # numberic sort of positions
      foreach my $end (sort { $a <=> $b } keys %{$change{$chrom}{$start}}) { # numberic sort of positions
	# get the clone position
	my @clone_coords = $coords->LocateSpan($chrom, $start, $end);
	my $clone = $clone_coords[0];

	if ($clone =~ /SUPERLINK/) {
	  my $ID = $change{$chrom}{$start}{$end}{'ID'};
	  $log->write_and_die("The location $ID is on the $clone at ",$clone_coords[1],"..",$clone_coords[2]," which will cause problems - maybe move the clone boundaries to avoid it?\n");
	}

	# see if this clone is on the side of the genome we want
	my $lab =  $clonelab{$clone};

	if ($lab eq 'RW' && $camace) {next}
	if ($lab eq 'HX' && $stlace) {next}

	if ($clone =~ /SUPERLINK/) {$errors .= "ID  $change{$chrom}{$start}{$end}{'ID'} is on $clone - skipping\n"; next}
	if ($clone =~ /CHROMOSOME/) {$errors .= "ID  $change{$chrom}{$start}{$end}{'ID'} is on $clone - skipping\n"; next}
	#$change{$chrom}{$start}{$end}{'ORDINAL'} = $ordinal;
	$ordered{$ordinal}{'clone'} = $clone_coords[0];
	$ordered{$ordinal}{'start_pos'} = $clone_coords[1];
	$ordered{$ordinal}{'end_pos'} = $clone_coords[2];
	$ordered{$ordinal}{'count_bases'} = abs ((length $change{$chrom}{$start}{$end}{'TOBASE'}) - (length $change{$chrom}{$start}{$end}{'FROMBASE'}));
	$ordered{$ordinal}{'change_type'} = $change{$chrom}{$start}{$end}{'TYPE'};
	$ordered{$ordinal}{'from_base'} = $change{$chrom}{$start}{$end}{'FROMBASE'};
	$ordered{$ordinal}{'to_base'} = $change{$chrom}{$start}{$end}{'TOBASE'};
	$ordered{$ordinal}{'ID'} = $change{$chrom}{$start}{$end}{'ID'};
	$ordered{$ordinal}{'chrom'} = $chrom; # copy over the chrom positions for debugging purposes
	$ordered{$ordinal}{'start'} = $start;
	$ordered{$ordinal}{'end'} = $end;

	$ordinal++;
	
      }
    }
  }
  return %ordered;
}


##########################################
#  &make_changes(%change);
# for the positions file input format
# do all the changes here

sub make_changes {
  my (%ordered) = @_;

  foreach my $ordinal (sort { $a <=> $b } keys %ordered) { # get each change in order down the genome
    
    # make the change in the clone sequence
    &position_change_DNA_in_db($ordered{$ordinal});


    # change the positions of objects in this clone
    #print "shift up gene models\n";
    if (&shift_gene_models($ordered{$ordinal})) {die "ERROR in shift_gene_models()\n";}


    # shift up/down the positions of other changes to be done in this clone
    &shift_locations($ordinal, \%ordered);

    # if the $ID column is the Feature id, then convert it into a 'corrected_genome_error' Feature
    &change_feature_method($ordered{$ordinal});

    # write the ACE changes for this location back to the database
    &write_ace_lines($ordered{$ordinal});
      

    # close then reopen an ACE connection to parse details for mapping to genome
    print "Connecting to Ace\n";
    $ace->close;
    $ace = Ace->connect (-path => $new_database,
			   -program => $tace) || die "cannot connect to database at $new_database\n";

    # save space
    delete $ordered{$ordinal};

    # dump out superlink sequences and check for non AGCT bases
    if ($dump_dna_every_time) {
      &check_superlink_sequence() || die "ERRORS were found in the superlink sequences\n";
    }

    print "===============================================================================================\n";

  }

}

##########################################
# &shift_locations($ordinal, $ordered{$ordinal});

sub shift_locations {

  my ($ordinal, $ordered) = @_;

  my $clone = $ordered->{$ordinal}{'clone'};
  my $original_clone_name = $ordered->{$ordinal}{'original_clone_name'};
  my $count_bases = $ordered->{$ordinal}{'count_bases'};
  my $type = $ordered->{$ordinal}{'change_type'};
  $ordinal++;

  while ((exists $ordered->{$ordinal}{'clone'} && $ordered->{$ordinal}{'clone'} eq $clone) || 
	 (exists $ordered->{$ordinal}{'clone'} && $ordered->{$ordinal}{'clone'} eq $original_clone_name)) {
    if (lc $type eq 'deletion') {
      $ordered->{$ordinal}{'start_pos'} -= $count_bases;
      $ordered->{$ordinal}{'end_pos'} -= $count_bases;

    } elsif (lc $type eq 'insertion') {
      print "Insertion start_pos ",$ordered->{$ordinal}{'start_pos'}," ",$ordered->{$ordinal}{'ID'}," changed by $count_bases to ";
      $ordered->{$ordinal}{'start_pos'} += $count_bases;
      print " new start_pos ",$ordered->{$ordinal}{'start_pos'},"\n";
      $ordered->{$ordinal}{'end_pos'} += $count_bases;


    }

    $ordinal++;
  }


}



##########################################
# &position_change_DNA_in_db($ordered{$ordinal});
# edit the DNA sequence and write it back to the database
sub position_change_DNA_in_db {
  my ($change) = @_;

  #print "get ace clone object for $clone\n";
  my $clone = $change->{'clone'};
  my $clone_obj = $ace->fetch(Sequence => $clone);
  my $ID = $change->{'ID'};
  my $type = $change->{'change_type'};
  print "Making change $type $ID to clone $clone\n\n";


  #print "reading DNA sequence\n";
  my $dna = $clone_obj->asDNA();
  $dna =~ s/\>(\w+)\n//;	# remove title line
  $dna =~ s/\n//g;		# make into one line
  #print "dna = $dna\n";

  &position_make_changes_to_DNA($change, $dna);

}


##########################################
# &position_make_changes_to_DNA($change, $dna);
# insertions are made to the right of the start_postion

sub position_make_changes_to_DNA {

  my ($change, $dna) = @_;

  my $clone = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};
  my $from_base = $change->{'from_base'};
  my $to_base = $change->{'to_base'};
  my $ID = $change->{'ID'};
  my $chrom = $change->{'chrom'};
  my $start = $change->{'start'};


  ######################
  # DELETION
  ######################

  if (lc($change_type) eq 'deletion') {
    my $clone_base = substr($dna, $start_pos-1, $count_bases);
    if (lc $clone_base ne lc $from_base) {die "In $ID deletion of $chrom $start bases $from_base to $to_base ($count_bases bases changed) : there is a mismatch to the clone base found: $clone_base\n"}
    substr($dna, $start_pos-1, $count_bases) = '';
  }

  ######################
  # INSERTION
  ######################

  elsif (lc($change_type) eq 'insertion') {
    substr($dna, $start_pos, 0) = lc($to_base);
  }
  
  ######################
  # SUBSTITUTION
  ######################
  
  elsif (lc($change_type) eq 'substitution' || lc($change_type) eq 'snp') {
    my $clone_base = substr($dna, $start_pos-1, length $change->{'from_base'});
    if (lc $clone_base ne lc $from_base) {die "In $ID substitution of $chrom $start bases $from_base to $to_base ($count_bases bases changed) : there is a mismatch to the clone base found: $clone_base\n"}
    substr($dna, $start_pos-1, length $change->{'from_base'}) = lc($to_base);
  }

  else {
    die "Unknown type of change: $change_type in $ID\n";
  }
  
  
  # reformat the DNA sequence into lines of 60 characters
  $dna = &reformat($dna);
  
  # write the dna back to the database
  push @{$change->{'ace-add'}{"DNA : \"$clone\""}}, "$dna"; # create output to be written to the ace file


}

##########################################
# if the $ID column is the Feature id, then convert it into a 'corrected_genome_error' Feature
# &change_feature_method($ordered{$ordinal});

sub change_feature_method {

  my ($change) = @_;

  my $ID = $change->{'ID'};

  if ($ID =~ /^WBsf\d+/) {
    push @{$change->{'ace-add'}{"Feature : \"$ID\""}}, "Method \"Corrected_genome_sequence_error\""; # create output to be written to the ace file
    push @{$change->{'ace-add'}{"Feature : \"$ID\""}}, "Remark \"This genome error was corrected in WS235.\""; # create output to be written to the ace file
  }

}

##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################

# ALL OF THE FOLLOWING WERE ORIGINALLY WRITTEN FOR THE FLANKING-SEQUENCE STYLE OF INPUT FILE

##########################################
# write_ace_lines($change);
# write out all the changes for this location and read them into the database

sub write_ace_lines {

  my ($change) = @_;

  # output all of the ace lines stored in $change->{'ace-delete'} to the output ACE file
  my $acefile = "/tmp/genome-changes$$.ace";
  print "write ace file of changes to $acefile\n";
  open (ACEOUT, ">$acefile") || die "Can't open file $acefile\n";
  foreach my $ace_object (keys %{$change->{'ace-delete'}}) {
    print ACEOUT "\n$ace_object\n";
    foreach my $ace_line (@{$change->{'ace-delete'}{$ace_object}}) {
      print ACEOUT "$ace_line\n";
    }
    print ACEOUT "\n";
  }
  # now output all of the ace lines stored in $change->{'ace-add'} to the output ACE file
  foreach my $ace_object (keys %{$change->{'ace-add'}}) {
    print ACEOUT "\n$ace_object\n";
    foreach my $ace_line (@{$change->{'ace-add'}{$ace_object}}) {
      print ACEOUT "$ace_line\n";
    }
    print ACEOUT "\n";
  }
  close(ACEOUT);
  
  
  if ($noload) {
    print "-noload specified. No changes to the new database will be made\n";
  } else {
    # load ace file into new database
    print "load ace file of changes into new database\n";
    $wormbase->load_to_database($new_database, $acefile, "genome_changes_pl", $log);
  }

}

##########################################
# result = &change_DNA_in_db($change)
# make the specified change to the DNA of the clone in the database
# return non-zero if hit a problem
# the input hash ref is updated to hold information on the site of the changes in the clone sequence

sub change_DNA_in_db {
  my ($change) = @_;

  my $clone = $change->{'clone'};
  my $region = $change->{'region'};

  #print "get ace clone object for $clone\n";
  my $clone_obj = $ace->fetch(Sequence => $clone);
  my $canonical = $clone_obj->at('Properties.Genomic_canonical');
  if ($canonical !~ /Genomic_canonical/) {
    $log->write_to( "*** The clone $clone is not genomic_canonical\nNo changes were made to $clone\n");
    return 1;			# error returned
  }

  #print "check it reads OK: ", $clone_obj->Overlap_left, "\n";
  
  #print "reading DNA sequence\n";
  my $dna = $clone_obj->asDNA();

  $dna =~ s/\>(\w+)\n//;	# remove title line
  $dna =~ s/\n//g;		# make into one line

  #print "dna = $dna\n";

  # find the part of the region we wish to change, 
  # i.e. the uppercased part of the region sequence
  my ($removed, @positions) = &find_changed_part_of_region($region);

  # make the changes to the DNA and store the result in the $change hash_ref
  my $error = &make_changes_to_DNA($change, $dna, $removed, @positions);

  if ($error == 2) {		# the region was not found - try the reverse complement
    my $region_rev = &DNA_string_reverse($region);
    $change->{'region'} = $region_rev;
    my ($removed_rev, @positions_rev) = &find_changed_part_of_region($region_rev);
    $error = &make_changes_to_DNA($change, $dna, $removed_rev, @positions_rev);
  }

  # still no sequence found
  if ($error == 2) {
    $log->write_to( "*** The sequence $removed was not found in $clone\nNo changes were made to $clone\n");
  }

  return $error;
}

##########################################
# make the desired changes to the DNA and store the result in in the $change hash_ref
# '2' is returned if the dna region is not found in the main sequence
# 1 is returned for all other errors
# $error = &make_changes_to_DNA($change, $dna, $removed, @positions);

sub make_changes_to_DNA {
  my ($change, $dna, $removed, @positions) = @_;

  my $clone = $change->{'clone'};
  my $region = $change->{'region'};
  my $change_type = $change->{'change_type'};

  # find the region we wish to change
  my $match_pos;


  ######################
  # DELETION
  ######################

  if (lc($change_type) eq 'deletion') {
    $match_pos = index $dna, lc($region);
    if ($match_pos == -1) {
      #$log->write_to( "*** The sequence $region was not found in $clone\nNo changes were made to $clone\n");
      return 2;			# error returned
    }
    # check the region is unique
    if ((my $pos2 = index $dna, lc($region), $match_pos+1) != -1) {
      $log->write_to( "*** The sequence $region is not unique in $clone e.g. positions: $match_pos and $pos2\nNo changes were made to $clone\n");
      return 1;			# error returned      
    }

    my $count = 0;
    foreach my $pos (@positions) {
      #print "pos=$count $pos\n";
      if ($pos == 1) {
	#print "delete at pos $match_pos + $count\n";
	# store the start and end positions of the change in the hash ref
	if (! exists $change->{'start_pos'}) {$change->{'start_pos'} = $match_pos + $count + 1;}
	$change->{'end_pos'} = $match_pos + $count + 1;
	$change->{'count_bases'}--;
	substr($dna, $match_pos + $count, 1) = '';
	$match_pos--;
      }	
      $count++;
    }

  ######################
  # INSERTION
  ######################

  } elsif (lc($change_type) eq 'insertion') {
    $match_pos = index $dna, lc($removed);
    if ($match_pos == -1) {
      return 2;			# error returned - sequence not found
    }
    # check the region is unique
    if ((my $pos2 = index $dna, lc($removed), $match_pos+1) != -1) {
      print "The sequence $removed is not unique in $clone e.g. positions: $match_pos and $pos2\nNo changes were made to $clone\n";
      return 1;			# error returned      
    }

    my $count = 0;
    foreach my $pos (@positions) {
      if ($pos == 1) {
	my $chr = uc(substr($region, $count, 1));
	#print "insert $chr at pos $match_pos + $count\n";
	# store the start and end positions of the change in the hash ref
	if (! exists $change->{'start_pos'}) {$change->{'start_pos'} = $match_pos + $count + 1;}
	$change->{'end_pos'} = $match_pos + $count + 1;
	$change->{'count_bases'}++;
	substr($dna, $match_pos + $count, 0) = lc(substr($region, $count, 1));
      }				
      $count++;
    }

  ###############
  # SHIFT-OVERLAP - this is where we wish to shift the superlink boundary without making a real insertion
  ###############

  } elsif (lc($change_type) eq 'shift-overlap') {
    $match_pos = index $dna, lc($removed);
    if ($match_pos == -1) {
      return 2;			# error returned - sequence not found
    }
    # check the region is unique
    if ((my $pos2 = index $dna, lc($removed), $match_pos+1) != -1) {
      print "The sequence $removed is not unique in $clone e.g. positions: $match_pos and $pos2\nNo changes were made to $clone\n";
      return 1;			# error returned      
    }

    my $count = 0;
    foreach my $pos (@positions) {
      if ($pos == 1) {
	my $chr = uc(substr($region, $count, 1));
	#print "insert $chr at pos $match_pos + $count\n";
	# store the start and end positions of the change in the hash ref
	if (! exists $change->{'start_pos'}) {$change->{'start_pos'} = $match_pos + $count + 1;}
	$change->{'end_pos'} = $match_pos + $count + 1;
	$change->{'count_bases'}++;
	substr($dna, $match_pos + $count, 0) = lc(substr($region, $count, 1));
      }				
      $count++;
    }




  } else {
    $log->write_to( "*** The type of change '$change_type' is not recognised\nNo changes were made to $clone\n");
    return 1;			# error returned
  }

  # see if there were any changes to be made
  if ($change->{'count_bases'} == 0) {
    $log->write_to( "*** No bases were in uppercase in the input '$region' in clone $clone\nNo changes were made to $clone\n");
    return 1;			# error returned
  }

# reformat the DNA sequence into lines of 60 characters
  $dna = &reformat($dna);

# write the dna back to the database
  push @{$change->{'ace-add'}{"DNA : \"$clone\""}}, "$dna"; # create output to be written to the ace file


  return 0;
}

##########################################
# find the part of the DNA region that is in uppercase and so marked to be changed
# returns the region sequence with the uppercase section removed
# and the array of positions to be changed
#
#  my ($removed, @positions) = &find_changed_part_of_region($region);

sub find_changed_part_of_region {
  my ($region) = @_;

  my @positions;
  my $count = 0;
  my $removed = '';		# region sequence with the uppercased bit removed
  foreach my $chr (split //, $region) {
    if (lc($chr) ne $chr) {
      $positions[$count] = 1;	# change this position
      #print "change at position $count, chr=$chr\n";
    } else {
      $positions[$count] = 0;
      $removed .= $chr;
    }
    $count++;
  }

  return ($removed, @positions);

}

##########################################
# $result = &shift_gene_models($change);
# change objects on the clone to shift them down
# return non-zero if hit a problem

sub shift_gene_models {
  my ($change) = @_;

  my $clone = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};

  # open tace connection to clone sequence and slurp up the contents
  my $cmd = "find Sequence $clone\nshow -a\nquit\n";
  open (TACE, "echo '$cmd' | $tace $new_database |");
  my @slurp = <TACE>;
  close TACE;

  #print "Contents of $clone:\n";
  #print @slurp;


  ##############################
  #
  # Clone data to be changed
  #
  ##############################

  # there is no point in shift things around if we are only doing a substitution
  if (lc($change_type) eq 'deletion' || lc($change_type) eq 'insertion') {

    # get and change Gene_child - works
    # NO - there are no Gene objects in camace
    #  print "change clone Genes\n";
    #  if (&get_and_change($change, 1, "Gene_child", @slurp)) {return 1;}

    # get and change CDS_child
    print "change clone CDSs\n";
    if (&get_and_change($change, 1, "CDS_child", @slurp)) {return 1;}

    # get and change Pseudogene
    print "change clone Pseudogenes\n";
    if (&get_and_change($change, 1, "Pseudogene", @slurp)) {return 1;}

    # get and change Transposon
    print "change clone Transposons\n";
    if (&get_and_change($change, 1, "Transposon", @slurp)) {return 1;}

    # get and change Transcript (these are probable non-coding RNA objects) - works
    print "change clone Transcripts\n";
    if (&get_and_change($change, 1, "Transcript", @slurp)) {return 1;}

    # get and change Feature_object (this is automatically done by the feature mapper, but it doesn't hurt to do it here as well)
    print "change clone Feature_objects\n";
    if (&get_and_change($change, 0, "Feature_object", @slurp)) {return 1;}

    # get into and change Feature_data - not important, just useful for display purposes
    print "change clone Feature_data\n";
    if ($nohomol) {
      print "NOT DONE DURING TESTING BECAUSE IT TAKES SO LONG\n";
    } else {
      if (&change_feature_data_on_clone($change, "Feature_data", @slurp)) {return 1;}
    }

    # get into and change Homol_data - not important, just useful for display purposes
    print "change clone Homol_data\n";
    if ($nohomol) {
      print "NOT DONE DURING TESTING BECAUSE IT TAKES SO LONG\n";
    } else {
      if (&change_homol_data_on_clone($change, "Homol_data", @slurp)) {return 1;}
    }

    # get and change Confirmed_intron - definitely want to update this
    print "change clone Confirmed_intron\n";
    if (&change_confirmed_intron($change, "Confirmed_intron", @slurp)) {return 1;}

    # Assembly tags are manual annotations done by the finishers.
    # they contain a text part (ignore this, it is simply a remark by
    # the finisher, like "clone right end", or "finished right")
    # followed by the start and end display position and the name of the
    # clone that the finisher was talking about.
    # Find the assembly tags with the clone name that we are working on
    # and change the assembly tag position in line with the changes.
    
    # get and change Assembly_tags - definitely want to update this
    print "change Assembly_tags\n";
    if (&change_assembly_tags($change, "Assembly_tags", @slurp)) {return 1;}
    
    
# these are not in camace - ignore them - they are mapped elsewhere - think about informing other sites?
#  # get and change Gene_child
#  # get and change Nongenomic
#  # get and change PCR_product
#  # get and change Allele
#  # get and change Oligo_set
#  # get and change SAGE_transcript
#  # get and change Oligo


  }

  # get and change Overlap_right - definitely want to update this
  if (&change_clone_overlaps($change, @slurp)) {return 1;}

  # get and change Remark - definitely want to update this - works
  if (&add_remark($change)) {return 1;}


  ##############################
  #
  # Superlink data to be changed
  #
  ##############################
  

  # get the name of the Superlink
  my ($source_line) = (grep /Source/, @slurp);
  my ($superlink) = ($source_line =~ /Source\s+\"(\S+)\"/);

  # open tace connection to superlink sequence and slurp up the contents
  $cmd = "find Sequence \"$superlink\"\nshow -a\nquit\n";
  open (TACE, "echo '$cmd' | $tace $new_database |");
  my @super_slurp = <TACE>;
  close TACE;

  #print "Contents of $superlink:\n";
  #print @super_slurp;

  # convert to using superlink coordinates instead of clone
  # coordinates for the changed region
  # and change the clone name to be the superlink name so that the
  # subroutines can be re-used easily
  if (&use_superlink_coords($change, $superlink, @super_slurp)) {return 1;};

  # shift all objects on the superlink downstream of the changed region
  # there is no point in shifting things around if we are only doing a substitution
  if (lc($change_type) eq 'deletion' || lc($change_type) eq 'insertion') {

    # get and change Gene_child
    # NO - there are no Gene objects in camace
    #  print "change Superlink Genes\n";
    #  if (&get_and_change($change, 1, "Gene_child", @super_slurp)) {return 1;}
    
    # get and change CDS_child
    print "change Superlink CDS\n";
    if (&get_and_change($change, 1, "CDS_child", @super_slurp)) {return 1;}
    
    # get and change Pseudogene
    print "change Superlink Pseudogene\n";
    if (&get_and_change($change, 1, "Pseudogene", @super_slurp)) {return 1;}
    
    # get and change Transposon
    print "change Superlink Transposon\n";
    if (&get_and_change($change, 1, "Transposon", @super_slurp)) {return 1;}
    
    # get and change Transcript (these are probable non-coding RNA objects)
    print "change Superlink Transcript\n";
    if (&get_and_change($change, 1, "Transcript", @super_slurp)) {return 1;}
    
    # get and change Feature_object (this is automatically done by the feature mapper, but it doesn't hurt to do it here as well)
    print "change Superlink Feature_object\n";
    if (&get_and_change($change, 0, "Feature_object", @super_slurp)) {return 1;}
    
    # get into and change Feature_data - not important, just useful for display purposes
    print "change Superlink Feature_data\n";
    if ($nohomol) {
      print "NOT DONE DURING TESTING BECAUSE IT TAKES SO LONG\n";
    } else {
      if (&change_feature_data_on_superlink($change, "Feature_data", @super_slurp)) {return 1;}
    }

    # get into and change Homol_data - not important, just useful for display purposes
    print "change Superlink Homol_data\n";
    if ($nohomol) {
      print "NOT DONE DURING TESTING BECAUSE IT TAKES SO LONG\n";
    } else {
      if (&change_homol_data_on_superlink($change, "Homol_data", @super_slurp)) {return 1;}
    }

    # get and change Confirmed_intron - definitely want to update this
    print "change Superlink Confirmed_intron\n";
    if (&change_confirmed_intron($change, "Confirmed_intron", @super_slurp)) {return 1;}

    if (&change_clone_length_on_superlink($change, $clone, @super_slurp)) {return 1;}
  }


  ##############################
  #
  # Chromosome data to be changed
  #
  ##############################
  

  # get the name of the Chromosome
  ($source_line) = (grep /Source/, @super_slurp);
  my ($chrom) = ($source_line =~ /Source\s+\"(\S+)\"/);

  # open tace connection to chrom sequence and slurp up the contents
  $cmd = "find Sequence \"$chrom\"\nshow -a\nquit\n";
  open (TACE, "echo '$cmd' | $tace $new_database |");
  my @chrom_slurp = <TACE>;
  close TACE;

  #print "Contents of $chrom:\n";
  #print @chrom_slurp;

  # convert to using chrom coordinates instead of clone
  # coordinates for the changed region
  # and change the clone name to be the chrom name so that the
  # subroutines can be re-used easily
  if (&use_chrom_coords($change, $superlink, $chrom, @chrom_slurp)) {return 1;};

  # shift all objects on the chrom downstream of the changed region
  # there is no point in shifting things around if we are only doing a substitution
  if (lc($change_type) eq 'deletion' || lc($change_type) eq 'insertion') {

    # get into and change Feature_data - not important, just useful for display purposes
    print "change Chrom Feature_data\n";
    if (&change_feature_data_on_chrom($change, "Feature_data", @chrom_slurp)) {return 1;}

    # get into and change Homol_data - not important, just useful for display purposes
    print "change Chrom Homol_data\n";
    if ($nohomol) {
      print "NOT DONE DURING TESTING BECAUSE IT TAKES SO LONG\n";
    } else {
      if (&change_homol_data_on_chrom($change, "Homol_data", @chrom_slurp)) {return 1;}
    }

    # get the length of the affected superlink on the chromosome object and
    # update its lengths
    print "change Chromosome superlink lengths\n";
    if (&change_superlink_length_on_chrom($change, $superlink, @chrom_slurp)) {return 1;}
  
    print "change Chromosome length\n";
    # this doesn't actually have any effect because we don't use
    # $change->{chrom_length} after this, but we may do further things
    # in the future
    &change_chromosome_length_on_chrom($change);
  }


# these are not in camace - ignore them - they are mapped elsewhere - think about informing other sites?
#  # get and change Gene_child
#  # get and change Nongenomic
#  # get and change PCR_product
#  # get and change Allele
#  # get and change Oligo_set
#  # get and change SAGE_transcript
#  # get and change Oligo


  return 0;
}

##########################################
# $result = &get_and_change($change, $important, @lines);
# work through the lines changing the ones that occur after the changed region
# If $important is non-zero then follow up any changes which occur within
# an object with a report to the log file
# return non-zero if hit a problem

sub get_and_change {

  my ($change, $important, $type, @lines) = @_;

  my $clone = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }

  #print "start, end and count of changed bases = $start_pos, $end_pos, $count_bases\n";

  foreach my $line (grep /^$type/, @lines) {
    chomp $line;

    my ($id, $start, $end) = ($line =~ /^$type\s+(\S+)\s+(\S+)\s+(\S+)/);
    next if (! defined $id || ! defined $start || ! defined $end);

    $id =~ s/\"//g;		# strip quotes from ID

    # debug
    #print "in get_and_change, line = $line\n";

    my $change_made = 0;	# a change has been made to the object's location
    my $within_object = 0;	# the change is within the object

    my $sense;
    if ($start < $end) {	# forward sense
      $sense = 1;

      # if we have a CDS in superlink, then only change twinscan CDS, not gene CDS or genefinder CDS
      # because these other CDSs are moved by the makesuperlinks.pl script later
      # (also want to move CDSs that are on the changed clone)
# now move up everything as we now do the superlink shuffling in this program rather than in makesuperlinks.pl
#      if (exists $change->{superlink} && 
#	  $type eq "CDS_child" && 
#	  $id !~ /tw$/ &&
#	  $start > $change->{clone_end_on_superlink}
#	  ) {next;}

      if ($end <= $start_pos) { # object before the change
	# no change
      } elsif ($start >= $end_pos) { # object after the change
	$start += $count_bases;
	$end += $count_bases;
	$change_made = 1;
      } elsif ($start <= $end_pos && $end >= $start_pos) { # object overlaps the changed region
	$within_object = 1;
	if ($start_pos <= $start) { # this will not make an accurate change if the start/end overlaps a changed region of several bases
	  $start += $count_bases;
	  $change_made = 1;
	}
	if ($end_pos <= $end) {     # this will not make an accurate change if the start/end overlaps a changed region of several bases
	  $end += $count_bases;
	  $change_made = 1;
	}
      }

    } elsif ($start > $end) {	# reverse sense
      $sense = -1;

      # if we have a CDS in superlink, then only change twinscan CDS, not gene CDS or genefinder CDS
      # because these other CDSs are moved by the makesuperlinks.pl script later
      # (also want to move CDSs that are on the changed clone)
# now move up everything as we now do the superlink shuffling in this program rather than in makesuperlinks.pl
#      if (exists $change->{superlink} && 
#	  $type eq "CDS_child" && 
#	  $id !~ /tw$/ &&
#	  $end > $change->{clone_end_on_superlink}       
#	  ) {next;}

      if ($start <= $start_pos) { # object before the change
	# no change
      } elsif ($end >= $end_pos) { # object after the change
	$start += $count_bases;
	$end += $count_bases;
	$change_made = 1;
      } elsif ($end <= $end_pos && $start >= $start_pos) { # object overlaps the changed region
	$within_object = 1;
	if ($start_pos <= $end) { # this will not make an accurate change if the start/end overlaps a changed region of several bases
	  $end += $count_bases;
	  $change_made = 1;
	}
	if ($end_pos <= $start) { # this will not make an accurate change if the start/end overlaps a changed region of several bases
	  $start += $count_bases;
	  $change_made = 1;
	}
      }
    } 

    # ouput to ace file
    if ($change_made) {
      push @{$change->{'ace-delete'}{"Sequence : \"$clone\""}}, "-D $line"; # create output to be written to the ace file
      #print "$line\n";
      push @{$change->{'ace-add'}{"Sequence : \"$clone\""}}, "$type \"$id\" $start $end"; # create output to be written to the ace file
      #print "$type \"$id\" $start $end\n";
    }

    # chase up changes within an object's exons
    if ($within_object && $type ne "Feature_object" && $type ne "Transposon") {	# transposons and features don't have exons
      if ($important) {
	$log->write_to( "*** Changes were made that overlapped with the $type\nobject '$id' in the clone $clone.\nYou should check whether '$id' is OK.\n\n");
      }
      # change the exons
      &change_exons($change, $type, $id, $start, $sense);
    }

    if ($type eq "Feature_object") { # need to change flanking sequences if the change is within about 30 bases of the feature object
      if (($start < $end && $end_pos > $start-30 && $start_pos < $end+30) ||    # forward sense
	  ($start > $end && $end_pos > $end-30 && $start_pos < $start+30)   ) { # reverse sense
	$log->write_to( "*** Changes were made near to the $type\nobject '$id' in the clone $clone.\nYou should check whether '$id' is OK. You may need to change the flanking sequences.\n\n");
	# change the flanking sequence
#	&change_flanking_sequences($change, $type, $id, $start, $end);
      }

    }

  }


  return 0;

}

##########################################
# attempt to shift up the exons of a CDS_child or Pseudogene
# when there is a change in the middle of the structure

sub change_exons {
  my ($change, $type, $id, $start, $sense) = @_; # 'start' is the start of the CDS on the clone

  my $clone = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }


  if ($type eq 'CDS_child') {
    $type = 'CDS';
  }

  my $cds_obj = $ace->fetch($type => $id);
  # then get the data following a tag in that object
  # and the NEXT data following the tag
  my @exons_start = $cds_obj->Source_exons; # these exon start positions are the position within the start to end on the clone
  my @exons_end = $cds_obj->Source_exons(2);
  $cds_obj->DESTROY();
  if (! @exons_start || ! @exons_end) {
    $log->write_to("Can't fetch exons for $id - no changes made to its exon structure\n");
    return;
  }

  my $warned = 0;
  my $changed = 0;
  for (my $i=0; $i < scalar @exons_start; $i++) {
    my $changed_this_exon = 0;
    my $newstart = $exons_start[$i];
    my $newend = $exons_end[$i];

    if ($sense == 1) {		# forward sense
#
# FORWARD SENSE     
#      
# CDS start   
# exon 1
#  |   
# exon 2
#  |   
# exon 3
#  |      <- change made here ($start_pos..$end_pos)
# exon 4   increment this exon's start & end
#  |   
# exon 5   increment this exon's start & end   
#  |   
# exon 6   increment this exon's start & end   
# CDS end - increment this exon's end position
#      
#      
      if ($exons_start[$i] >= $end_pos - $start + 1) {
	$newstart += $count_bases;
	$changed++;
	$changed_this_exon++;
      }
      if ($exons_end[$i] >= $end_pos - $start + 1) {
	$newend += $count_bases;
	$changed++;
	$changed_this_exon++;
      }
    } else {			# reverse sense

      # N.B. we need to shift along the exons BEFORE (genomic) and
      # AFTER (counting the exons) the change so that they fill the
      # whole of the CDS region
#      
# REVERSE SENSE     
#      
# CDS end   
# exon 6   increment this exon's start & end
#  |   
# exon 5   increment this exon's start & end
#  |   
# exon 4   increment this exon's start & end
#  |      <- change made here ($start_pos..$end_pos)
# exon 3
#  |   
# exon 2   
#  |   
# exon 1   
# CDS start - increment this exon's start position
#      
#      
#      
      
      if ($exons_start[$i] >= $start - $end_pos + 1) {
	$newstart += $count_bases;
	$changed++;
	$changed_this_exon++;
      }
      if ($exons_end[$i] >= $start - $end_pos + 1) {
	$newend += $count_bases;
	$changed++;
	$changed_this_exon++;
      }
    }

    if ($changed_this_exon) {
      push @{$change->{'ace-delete'}{"$type : \"$id\""}}, "-D Source_exons $exons_start[$i] $exons_end[$i]";
      push @{$change->{'ace-add'}{"$type : \"$id\""}}, "Source_exons $newstart $newend";
    }

    if ($changed == 1) {
      $log->write_to( "* (Changes were made that overlapped with exon $i of '$id')\n\n");
      $warned = 1;		# we have now warned the user - don't display a message for this gene again
    } elsif ($changed == 2 && !$warned) {
      $log->write_to( "* (Changes were made in the intron before exon $i of '$id')\n\n");
      $warned = 1;		# we have now warned the user - don't display this message for this gene again
    }

  }



}
##########################################
# attempt to repair the flanking sequences of a Feature_object
# when there is a change near to that object
sub change_flanking_sequences {
  my ($change, $type, $id, $start, $end) = @_;

  my $clone = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }


  my $feature_obj = $ace->fetch(Feature => $id);
  # then get the flanking sequences
  my $flank_clone = $feature_obj->Flanking_sequences;
  my $flank1 = $feature_obj->Flanking_sequences(2);
  my $flank2 = $feature_obj->Flanking_sequences(3);

# +++ here we need to match the section of flanking sequence to the
# change and apply the change to the flanking sequence

  if ($end_pos > $start-40 && $start_pos < $start) { # change in first flank
    $log->write_to( "*** Changes were made that in the first flanking sequence '$flank1' of $type\nobject '$id' in the clone $clone.\nYou should correct this sequence.\n\n");

  } elsif ($start_pos < $end+40 && $end_pos > $end) {
    $log->write_to( "*** Changes were made that in the second flanking sequence '$flank2' of $type\nobject '$id' in the clone $clone.\nYou should correct this sequence.\n\n");

  }


}
##########################################
# $result = &change_feature_data_on_clone($change, @lines);
# change the clone's Feature_data line to reflect the new clone length
# change the positions in the Feature_data object to reflect the changes
# return non-zero if hit a problem

sub change_feature_data_on_clone {

  my ($change, $type, @lines) = @_;

  my $clone = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }


  foreach my $line (grep /^$type/, @lines) {
    chomp $line;

    my ($id, $clone_length) = ($line =~ /^$type\s+\"(\S+)\"\s+1\s+(\d+)/);
    next if (! defined $clone_length);

    #print "in change_feature_data(), line = $line\n";
    #print "id = $id\n";
    #print "clone_length = $clone_length\n";

    # update the clone length
    $clone_length += $count_bases;

    # print to ace file
    push @{$change->{'ace-delete'}{"Sequence : \"$clone\""}}, "-D $line"; # create output to be written to the ace file
    #print "$line\n";
    push @{$change->{'ace-add'}{"Sequence : \"$clone\""}}, "Feature_data $id 1 $clone_length"; # create output to be written to the ace file
    #print "Feature_data $id 1 $clone_length\n";


    # open tace connection to feature_data object and slurp up the contents
    my $cmd = "find Feature_data $id\nshow -a\nquit\n";
    open (TACE, "echo '$cmd' | $tace $new_database |");
    my @slurp = <TACE>;
    close TACE;


    # change the start and end positions if they are past the changed
    # region. Any changes will be repaired when the feature is next mapped to the genome    
    foreach my $feature (grep /^Feature/, @slurp) {
      chomp $feature;
      if ($feature =~ /^Feature_data/) {next;} # ignore the Feature_data line, we just want the 'Feature' lines
      my @split = split /\s+/, $feature;
      if ($split[2] > $end_pos || $split[3] > $end_pos) {
	$split[2] += $count_bases;
	$split[3] += $count_bases;
	my $new_feature = join " ", @split;
	push @{$change->{'ace-delete'}{"Feature_data : \"$id\""}}, "-D $feature"; # create output to be written to the ace file
	#print "$feature\n";
	push @{$change->{'ace-add'}{"Feature_data : \"$id\""}}, "$new_feature"; # create output to be written to the ace file
	#print "$new_feature\n";
      }
    }
  }

  return 0;
}

##########################################
# $result = &change_feature_data_on_sperlink($change, @lines);
# change the superlinks's Feature_data line to reflect the new clone length
# change the positions in the Feature_data object to reflect the changes
# return non-zero if hit a problem

sub change_feature_data_on_superlink {

  my ($change, $type, @lines) = @_;

  my $clone = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }


  my $have_changed_the_features = 0;

  my @grepped_lines = grep /^$type/, @lines;
  foreach my $line (@grepped_lines) {
    chomp $line;


    #print "$line";

    my ($id, $feat_start, $feat_end) = ($line =~ /^$type\s+\"(\S+)\"\s+(\d+)\s+(\d+)/);
    next if (! defined $feat_start || ! defined $feat_end);

    # if we have a superlink, then only change the last feature-data of a
    # set of virtual blocks

    #print "looking for superlink last feature-data in $line\n";
    #print "ID $id, feat_start $feat_start, feat_end $feat_end\n";
    if ($feat_end == $change->{superlink_length}) {
      #print "Found the last feature-data in the superlink for $id\n";
      # update the clone length
      $feat_end += $count_bases;
      # print to ace file
      push @{$change->{'ace-delete'}{"Sequence : \"$clone\""}}, "-D $line"; 
      #print "$line\n";
      push @{$change->{'ace-add'}{"Sequence : \"$clone\""}}, "Feature_data $id $feat_start $feat_end"; 
      #print "Feature_data $id $feat_start $feat_end\n";
      $have_changed_the_features = 1;
    }

    # change the start and end positions if they are past the changed
    # region - we need not get too excited about changes occurring in
    # the middle of the feature_data object
    #print " feat_end $feat_end  start_pos $start_pos\n";
    if ($feat_end > $start_pos) { # only work on this feature_data object if any of it is past the changed region

      # try to get the data in the Feature-data virtual block faster by using ACE rather than tace
      my $feature_obj = $ace->fetch(Feature_data => $id);
      foreach my $intron_obj ($feature_obj->at('Splices.Confirmed_intron')) {
	my ($intron_start, $intron_end, $intron_type, $intron_id) = $intron_obj->row;
	my $feat = "Confirmed_intron $intron_start $intron_end $intron_type $intron_id";
	#print "$feat\n";
	if (! defined $intron_start || ! defined $intron_end) {next;}
        # don't worry about getting a change in the middle of a feature correct
	if ($feat_start + $intron_start - 1 > $start_pos || # test if have a feature past the start pos of the change
	    $feat_start + $intron_end - 1 > $start_pos) {  
	  $intron_start += $count_bases;
	  $intron_end += $count_bases;
	  my $new_feat = "Confirmed_intron $intron_start $intron_end $intron_type $intron_id";
	  #print "$new_feat\n";
	  push @{$change->{'ace-delete'}{"Feature_data : \"$id\""}}, "-D $feat"; # create output to be written to the ace file
	  push @{$change->{'ace-add'}{"Feature_data : \"$id\""}}, "$new_feat"; # create output to be written to the ace file
	}
      }



      # open tace connection to feature_data object and slurp up the contents
#      my $cmd = "find Feature_data $id\nshow -a\nquit\n";
#      open (TACE, "echo '$cmd' | $tace $database |");
#      my @slurp = <TACE>;
#      close TACE;

#      foreach my $feat (grep /^Confirmed_intron/, @slurp) {
#	#print "next feat in slurp = $feat\n";
#	chomp $feat;
#	my @split = split /\s+/, $feat;
#	if (! defined $split[1] || ! defined $split[2]) {next;} # if the start or end position is not present in the Feature_data object, skip it
#        # don't worry about getting a change in the middle of a feature correct
#	if ($feat_start + $split[1] - 1 > $start_pos || # test if have a feature past the start pos of the change
#	    $feat_start + $split[2] - 1 > $start_pos) { 
#	  $split[1] += $count_bases;
#	  $split[2] += $count_bases;
#	  my $new_feat = join " ", @split;
#	  push @{$change->{'ace-delete'}{"Feature_data : \"$id\""}}, "-D $feat"; # create output to be written to the ace file
#	  #print "feat_start = $feat_start split[1] = $split[1] start_pos =$start_pos \n";
#	  #print "$feat\n";
#	  push @{$change->{'ace-add'}{"Feature_data : \"$id\""}}, "$new_feat"; # create output to be written to the ace file
#	  #print "$new_feat\n";
#	}
#      }
    }

  }

  # sanity check for when doing superlinks
  if (@grepped_lines && ! $have_changed_the_features) {

    $log->write_to( "*** For some reason we did not manage to change the feature_data
lengths of the superlink ${clone}'s last virtual feature-data block.
Check this. No changes were made.\n");
    return 1;
  }

  return 0;
}




##########################################
# $result = &change_feature_data_on_chrom($change, @lines);
# change the chroms's Feature_data line to reflect the new clone length
# change the positions in the Feature_data object to reflect the changes
# return non-zero if hit a problem

sub change_feature_data_on_chrom {

  my ($change, $type, @lines) = @_;

  my $clone = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }


  my $have_changed_the_features = 0;

  my @grepped_lines = grep /^$type/, @lines;
  foreach my $line (@grepped_lines) {
    chomp $line;


    #print "$line";

    my ($id, $feat_start, $feat_end) = ($line =~ /^$type\s+\"(\S+)\"\s+(\d+)\s+(\d+)/);
    next if (! defined $feat_start || ! defined $feat_end);

    # if we have a chrom, then only change the last feature-data of a
    # set of virtual blocks

    #print "looking for chrom last feature-data in $line\n";
    #print "ID $id, feat_start $feat_start, feat_end $feat_end\n";
    if ($feat_end == $change->{chrom_length}) {
      #print "Found the last feature-data in the chrom for $id\n";
      # update the clone length
      $feat_end += $count_bases;
      # print to ace file
      push @{$change->{'ace-delete'}{"Sequence : \"$clone\""}}, "-D $line"; 
      #print "$line\n";
      push @{$change->{'ace-add'}{"Sequence : \"$clone\""}}, "Feature_data $id $feat_start $feat_end"; 
      #print "Feature_data $id $feat_start $feat_end\n";
      $have_changed_the_features = 1;
    }

    # change the start and end positions if they are past the changed
    # region - we need not get too excited about changes occurring in
    # the middle of the feature_data object
    #print " feat_end $feat_end  start_pos $start_pos\n";
    if ($feat_end > $start_pos) { # only work on this feature_data object if any of it is past the changed region

      # try to get the data in the Feature-data virtual block faster by using ACE rather than tace
      my $feature_obj = $ace->fetch(Feature_data => $id);
      foreach my $intron_obj ($feature_obj->at('Splices.Confirmed_intron')) {
	my ($intron_start, $intron_end, $intron_type, $intron_id) = $intron_obj->row;
	my $feat = "Confirmed_intron $intron_start $intron_end $intron_type $intron_id";
	#print "$feat\n";
	if (! defined $intron_start || ! defined $intron_end) {next;}
        # don't worry about getting a change in the middle of a feature correct
	if ($feat_start + $intron_start - 1 > $start_pos || # test if have a feature past the start pos of the change
	    $feat_start + $intron_end - 1 > $start_pos) {  
	  $intron_start += $count_bases;
	  $intron_end += $count_bases;
	  my $new_feat = "Confirmed_intron $intron_start $intron_end $intron_type $intron_id";
	  #print "$new_feat\n";
	  push @{$change->{'ace-delete'}{"Feature_data : \"$id\""}}, "-D $feat"; # create output to be written to the ace file
	  push @{$change->{'ace-add'}{"Feature_data : \"$id\""}}, "$new_feat"; # create output to be written to the ace file
	}
      }

    }

  }

  # sanity check for when doing chroms
  if (@grepped_lines && ! $have_changed_the_features) {

    $log->write_to( "*** For some reason we did not manage to change the feature_data
lengths of the chrom ${clone}'s last virtual feature-data block.
Check this. No changes were made.\n");
    return 1;
  }

  return 0;
}




##########################################
# $result = &change_homol_data_on_clone($change, @lines);
# change the clone's Homol_data lines to reflect the new clone length
# change the positions in the Feature_data object to reflect the changes
# return non-zero if hit a problem

sub change_homol_data_on_clone {

  my ($change, $type, @lines) = @_;

  my $clone = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }


  my $have_changed_the_homols = 0;

  foreach my $line (grep /^$type/, @lines) {
    chomp $line;

    #print "$line";

    my ($id, $homol_start, $homol_end) = ($line =~ /^$type\s+\"(\S+)\"\s+(\d+)\s+(\d+)/);
    next if (! defined $homol_start || ! defined $homol_end);

    # update the clone length
    $homol_end += $count_bases;
    # print to ace file
    push @{$change->{'ace-delete'}{"Sequence : \"$clone\""}}, "-D $line"; # create output to be written to the ace file
    #print "$line\n";
    push @{$change->{'ace-add'}{"Sequence : \"$clone\""}}, "Homol_data $id $homol_start $homol_end"; # create output to be written to the ace file
    #print "Homol_data $id $homol_start $homol_end\n";

    # open tace connection to homol_data object and slurp up the contents
    my $cmd = "find Homol_data $id\nshow -a\nquit\n";
    open (TACE, "echo '$cmd' | $tace $new_database |");
    my @slurp = <TACE>;
    close TACE;

    # change the start and end positions if they are past the changed
    # region - we need not get too excited about changes occurring in
    # the middle of the homol_data object
    foreach my $homol (grep /homol/, @slurp) {
      chomp $homol;
      my @split = split /\s+/, $homol;
      if (! defined $split[4] || ! defined $split[5]) {next;} # if the start or end position is not present in the Homol_data object, skip it
      if ($split[4] > $end_pos || $split[5] > $end_pos) {
	$split[4] += $count_bases;
	$split[5] += $count_bases;
	# deal with the AlignPep bits, if present
	if (defined $split[8] && $split[8] =~ /Align/) {
	  $split[9] += $count_bases;
	}
	my $new_homol = join " ", @split;
	push @{$change->{'ace-delete'}{"Homol_data : \"$id\""}}, "-D $homol"; # create output to be written to the ace file
	push @{$change->{'ace-add'}{"Homol_data : \"$id\""}}, "$new_homol"; # create output to be written to the ace file
      }
    }
                                                                                                                                          
  }

  return 0;
}
##########################################
# $result = &change_homol_data_on_superlink($change, @lines);
# change the clone's Homol_data lines to reflect the new clone length
# change the positions in the Feature_data object to reflect the changes
# return non-zero if hit a problem

sub change_homol_data_on_superlink {

  my ($change, $type, @lines) = @_;

  my $clone = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }


  my $have_changed_the_homols = 0;

  my @grepped_lines = grep /^$type/, @lines;
  foreach my $line (@grepped_lines) {
    chomp $line;

    #print "$line";

    my ($id, $homol_start, $homol_end) = ($line =~ /^$type\s+\"(\S+)\"\s+(\d+)\s+(\d+)/);
    next if (! defined $homol_start || ! defined $homol_end);

# if we have a superlink, then only change the last homol-data of a
# set of virtual blocks
    #print "looking for superlink last homol in $line\n";
    #print "ID $id, homol_start $homol_start, homol_end $homol_end\n";
    if ($homol_end == $change->{superlink_length}) {
      #print "Found the last homol in the superlink for $id\n";
      # update the clone length
      $homol_end += $count_bases;
      # print to ace file
      push @{$change->{'ace-delete'}{"Sequence : \"$clone\""}}, "-D $line"; 
      #print "$line\n";
      push @{$change->{'ace-add'}{"Sequence : \"$clone\""}}, "Homol_data $id $homol_start $homol_end"; 
      #print "Homol_data $id $homol_start $homol_end\n";
      $have_changed_the_homols = 1;
    }

    # change the start and end positions if they are past the changed
    # region - we need not get too excited about changes occurring in
    # the middle of the homol_data object
    if ($homol_end > $start_pos) { # only work on this homol_data object if any of it is past the changed region

# this only changes the yk82a9.5 line:
# DNA_homol yk82a9.5 BLAT_EST_BEST 97.8 26298 26394 1 97
# when doing the F52D10 test data, why doesn't it find the other two?
#
#      # try to get the data in the Homol-data virtual block faster by using ACE rather than tace
#      my $homol_obj = $ace->fetch(Homol_data => $id);
#      foreach my $blat_obj ($homol_obj->at('Homol.DNA_homol')) {
#	my ($blat_id, $blat_type, $blat_score, $blat_start, $blat_end, $other_start, $other_end) = $blat_obj->row;
#	my $blat = "DNA_homol $blat_id $blat_type $blat_score $blat_start $blat_end $other_start $other_end";
#	print "$blat\n" if ($blat_id eq "yk82a9.5");
#	if (! defined $blat_start || ! defined $blat_end) {next;}
#        # don't worry about getting a change in the middle of a feature correct
#	if ($homol_start + $blat_start - 1 > $start_pos || # test if have a feature past the start pos of the change
#	    $homol_start + $blat_end - 1 > $start_pos) {  
#	  $blat_start += $count_bases;
#	  $blat_end += $count_bases;
#	  my $new_blat = "DNA_homol $blat_id $blat_type $blat_score $blat_start $blat_end $other_start $other_end";
#	  print "$new_blat\n" if ($blat_id eq "yk82a9.5");
#	  push @{$change->{'ace-delete'}{"Homol_data : \"$id\""}}, "-D $blat"; # create output to be written to the ace file
#	  push @{$change->{'ace-add'}{"Homol_data : \"$id\""}}, "$new_blat"; # create output to be written to the ace file
#	}
#      }


      # open tace connection to homol_data object and slurp up the contents
      my $cmd = "find Homol_data $id\nshow -a\nquit\n";
      open (TACE, "echo '$cmd' | $tace $new_database |");
      my @slurp = <TACE>;
      close TACE;

      foreach my $homol (grep /homol/, @slurp) {
	chomp $homol;
	my @split = split /\s+/, $homol;
	if (! defined $split[4] || ! defined $split[5]) {next;} # if the start or end position is not present in the Homol_data object, skip it
	#print "$homol\n" if ($split[1] eq "\"yk82a9.5\"");
        # don't worry about getting a change in the middle of a homol correct
	if ($homol_start + $split[4] - 1 > $start_pos || # test if have a homol past the start pos of the change
	    $homol_start + $split[5] - 1 > $start_pos) { 
	  #print "Changing this line\n" if ($split[1] eq "\"yk82a9.5\"");

	  $split[4] += $count_bases;
	  $split[5] += $count_bases;
	  # deal with the AlignPep bits, if present
	  if (defined $split[8] && $split[8] =~ /Align/) {
	    $split[9] += $count_bases;
	  }
	  my $new_homol = join " ", @split;
	  push @{$change->{'ace-delete'}{"Homol_data : \"$id\""}}, "-D $homol"; # create output to be written to the ace file
	  #print "$homol\n";
	  push @{$change->{'ace-add'}{"Homol_data : \"$id\""}}, "$new_homol"; # create output to be written to the ace file
	  #print "$new_homol\n";
	}
      }




    }
  }

  # sanity check for when doing superlinks
  if (@grepped_lines && ! $have_changed_the_homols) {

    $log->write_to( "*** For some reason we did not manage to change the homol_data
lengths of the superlink ${clone}'s last virtual block.
Check this. No changes were made.\n");
    return 1;
  }

  return 0;
}

##########################################
# $result = &change_homol_data_on_chrom($change, @lines);
# change the clone's Homol_data lines to reflect the new clone length
# change the positions in the Feature_data object to reflect the changes
# return non-zero if hit a problem

sub change_homol_data_on_chrom {

  my ($change, $type, @lines) = @_;

  my $clone = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }


  my $have_changed_the_homols = 0;

  my @grepped_lines = grep /^$type/, @lines;
  foreach my $line (@grepped_lines) {
    chomp $line;

    #print "$line";

    my ($id, $homol_start, $homol_end) = ($line =~ /^$type\s+\"(\S+)\"\s+(\d+)\s+(\d+)/);
    next if (! defined $homol_start || ! defined $homol_end);

# if we have a chrom, then only change the last homol-data of a
# set of virtual blocks
    #print "looking for chrom last homol in $line\n";
    #print "ID $id, homol_start $homol_start, homol_end $homol_end\n";
    if ($homol_end == $change->{chrom_length}) {
      #print "Found the last homol in the chrom for $id\n";
      # update the clone length
      $homol_end += $count_bases;
      # print to ace file
      push @{$change->{'ace-delete'}{"Sequence : \"$clone\""}}, "-D $line"; 
      #print "$line\n";
      push @{$change->{'ace-add'}{"Sequence : \"$clone\""}}, "Homol_data $id $homol_start $homol_end"; 
      #print "Homol_data $id $homol_start $homol_end\n";
      $have_changed_the_homols = 1;
    }

    # change the start and end positions if they are past the changed
    # region - we need not get too excited about changes occurring in
    # the middle of the homol_data object
    if ($homol_end > $start_pos) { # only work on this homol_data object if any of it is past the changed region

      # open tace connection to homol_data object and slurp up the contents
      my $cmd = "find Homol_data $id\nshow -a\nquit\n";
      open (TACE, "echo '$cmd' | $tace $new_database |");
      my @slurp = <TACE>;
      close TACE;

      foreach my $homol (grep /homol/, @slurp) {
	chomp $homol;
	my @split = split /\s+/, $homol;
	if (! defined $split[4] || ! defined $split[5]) {next;} # if the start or end position is not present in the Homol_data object, skip it
	#print "$homol\n" if ($split[1] eq "\"yk82a9.5\"");
        # don't worry about getting a change in the middle of a homol correct
	if ($homol_start + $split[4] - 1 > $start_pos || # test if have a homol past the start pos of the change
	    $homol_start + $split[5] - 1 > $start_pos) { 
	  #print "Changing this line\n" if ($split[1] eq "\"yk82a9.5\"");

	  $split[4] += $count_bases;
	  $split[5] += $count_bases;
	  # deal with the AlignPep bits, if present
	  if (defined $split[8] && $split[8] =~ /Align\s/) {
	    $split[9] += $count_bases;
	  }
	  my $new_homol = join " ", @split;
	  push @{$change->{'ace-delete'}{"Homol_data : \"$id\""}}, "-D $homol"; # create output to be written to the ace file
	  #print "$homol\n";
	  push @{$change->{'ace-add'}{"Homol_data : \"$id\""}}, "$new_homol"; # create output to be written to the ace file
	  #print "$new_homol\n";
	}
      }




    }
  }

  # sanity check for when doing chroms
  if (@grepped_lines && ! $have_changed_the_homols) {

    $log->write_to( "*** For some reason we did not manage to change the homol_data
lengths of the chrom ${clone}'s last virtual block.
Check this. No changes were made.\n");
    return 1;
  }

  return 0;
}


##########################################
# change_clone_length_on_superlink($change, $clone, @super_slurp)
# change the superlink's Subsequence clone line to reflect the new clone length
# return non-zero if hit a problem

sub change_clone_length_on_superlink {

  my ($change, $clone, @lines) = @_;


  my $superlink = $change->{'clone'};
  #my $clone = $change->{'original_clone_name'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }


  my $have_changed_the_length = 0;

  my @grepped_lines = grep /^Subsequence/, @lines;
  foreach my $line (@grepped_lines) {
    chomp $line;

    #print "$line";

    my ($id, $sl_start, $sl_end) = ($line =~ /^Subsequence\s+\"(\S+)\"\s+(\d+)\s+(\d+)/);
    next if (! defined $sl_start || ! defined $sl_end);

    if ($id eq $clone) {
      $sl_end += $count_bases;
      # print to ace file
      push @{$change->{'ace-delete'}{"Sequence : \"$superlink\""}}, "-D $line"; 
      #print "$line\n";
      push @{$change->{'ace-add'}{"Sequence : \"$superlink\""}}, "Subsequence $id $sl_start $sl_end"; 
      #print "Homol_data $id $homol_start $homol_end\n";
      $have_changed_the_length = 1;
    }

    # now change the downstream superlinks
    if ($sl_start > $start_pos) {
      $sl_start += $count_bases;
      $sl_end += $count_bases;
      # print to ace file
      push @{$change->{'ace-delete'}{"Sequence : \"$superlink\""}}, "-D $line"; 
      #print "$line\n";
      push @{$change->{'ace-add'}{"Sequence : \"$superlink\""}}, "Subsequence $id $sl_start $sl_end"; 
      #print "Homol_data $id $homol_start $homol_end\n";
      $have_changed_the_length = 1;
    }

  }

  # sanity check for when doing superlinks
  if (@grepped_lines && ! $have_changed_the_length) {

    $log->write_to( "*** For some reason we did not manage to change the superlink length of the superlink ${superlink}. Check this.\n");
    return 1;
  }

  return 0;

}

##########################################
# $result = &change_superlink_length_on_chrom($change, $superlink, @chrom_slurp)
# change the chromosome's Subsequence superlink line to reflect the new superlink length
# return non-zero if hit a problem

sub change_superlink_length_on_chrom {

  my ($change, $superlink, @lines) = @_;

  my $chrom = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }


  my $have_changed_the_length = 0;

  my @grepped_lines = grep /^Subsequence/, @lines;
  foreach my $line (@grepped_lines) {
    chomp $line;

    #print "$line";

    my ($id, $sl_start, $sl_end) = ($line =~ /^Subsequence\s+\"(\S+)\"\s+(\d+)\s+(\d+)/);
    next if (! defined $sl_start || ! defined $sl_end);

    if ($id eq $superlink) {
      $sl_end += $count_bases;
      # print to ace file
      push @{$change->{'ace-delete'}{"Sequence : \"$chrom\""}}, "-D $line"; 
      #print "$line\n";
      push @{$change->{'ace-add'}{"Sequence : \"$chrom\""}}, "Subsequence $id $sl_start $sl_end"; 
      #print "Homol_data $id $homol_start $homol_end\n";
      $have_changed_the_length = 1;
    }

    # now change the downstream superlinks
    if ($sl_start > $start_pos) {
      $sl_start += $count_bases;
      $sl_end += $count_bases;
      # print to ace file
      push @{$change->{'ace-delete'}{"Sequence : \"$chrom\""}}, "-D $line"; 
      #print "$line\n";
      push @{$change->{'ace-add'}{"Sequence : \"$chrom\""}}, "Subsequence $id $sl_start $sl_end"; 
      #print "Homol_data $id $homol_start $homol_end\n";
      $have_changed_the_length = 1;
    }

  }

  # sanity check for when doing chroms
  if (@grepped_lines && ! $have_changed_the_length) {

    $log->write_to( "*** For some reason we did not manage to change the superlink
length of the chrom ${chrom}.
Check this. No changes were made.\n");
    return 1;
  }

  return 0;
}

##########################################
#  &change_chromosome_length_on_chrom($change);
# make an explicit change to the length of the chromosome stored on
# the chromosome object here

sub change_chromosome_length_on_chrom {
  my ($change) = @_;


  my $chrom = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $count_bases = $change->{'count_bases'};

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }


  $change->{chrom_length} += $count_bases;

  my $chrom_length = $change->{chrom_length};

  ## this would make a ghost chromosome DNA object.
  #push @{$change->{'ace-add'}{"Sequence : \"$chrom\""}}, "DNA \"$chrom\" $chrom_length";

}


##########################################
# $result = &change_confirmed_intron($change, $clone, @lines);
# change the clone's Homol_data lines to reflect the new clone length
# change the positions in the Feature_data object to reflect the changes
# return non-zero if hit a problem
# This information is on the lines like: "Splices Confirmed_intron 20575 19947 UTR yk1308c03.5"
# on the clone's Sequence object.
# Only present for confirmed UTR introns e.g. in clone C12D8


sub change_confirmed_intron {

  my ($change, $type, @lines) = @_;

  my $clone = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }


  my $newline;

  foreach my $line (grep /^$type/, @lines) {
    chomp $line;

    #print "$line\n";

    # if the change occurs in the middle of the intron, we simply
    # change the length of the intron by only altering the position
    # after the change.
    my @split = split /\s+/, $line; 
    if (! defined $split[1] || ! defined $split[2]) {next;} # if the start or end position is not present, skip it
    if ($split[1] !~ /^\d+$/ || $split[2] !~ /^\d+$/ ) {next;} # if the start or end position is not numeric, skip it

    if ($split[1] > $end_pos || $split[2] > $end_pos) {
      if ($split[1] > $end_pos ) {
	$split[1] += $count_bases;
      }
      if ($split[2] > $end_pos) {
	$split[2] += $count_bases;
      }
      my $newline = join " ", @split;

      # print to ace file
      push @{$change->{'ace-delete'}{"Sequence : \"$clone\""}}, "-D $line"; # create output to be written to the ace file
      #print "$line\n";
      push @{$change->{'ace-add'}{"Sequence : \"$clone\""}}, "$newline"; # create output to be written to the ace file
      #print "$newline\n";
    }

  }

  return 0;
}

##########################################
# $result = &change_assembly_tags($change, @lines);
# Assembly tags are manual annotations done by the finishers.
# they contain a text part (ignore this, it is simply a remark by
# the finisher, like "clone right end", or "finished right")
# followed by the start and end display position and the name of the
# clone that the finisher was talking about.
# Find the assembly tags with the clone name that we are working on
# and change the assembly tag position in line with the changes.

sub change_assembly_tags {

  my ($change, $type, @lines) = @_;

  my $clone = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }


  my $newline;

  foreach my $line (grep /^$type/, @lines) {
    chomp $line;

    #print "$line\n";

    # if the change occurs in the middle of the tag, we need to raise an error
    # because changes shouldn't occur there - plenty of people will have
    # scrutinised the region and a change here is probably an error
    my @split = ($line =~ /($type)\s+(\".*?\")\s+(\d+)\s+(\d+)\s+(\".*?\")/);
    if (! defined $split[2] || ! defined $split[3]) {next;} # if the start or end position is not present, skip it
    #if (lc($split[4]) ne lc("\"$clone\"")) {next;} # we are only interested in tags on our clone
    # test for overlap (must test for assembly tag being in either sense)
    if (
	(($split[2] < $split[3]) && ($split[2] < $end_pos && $split[3] > $start_pos)) ||
	(($split[3] < $split[2]) && ($split[3] < $end_pos && $split[2] > $start_pos))
	) {
      $log->write_to( "*** A change has been attempted in a region covered by an assembly tag:\n");
      $log->write_to("$split[4]\n");
      $log->write_to( "in clone $clone, position $start_pos..$end_pos\n");
      $log->write_to( "You should investigate this change to see if it is real.\n");
      return 1 unless ($ignore_assembly_tags);
    }
    if ($split[2] > $end_pos && $split[3] > $end_pos) {
      $split[2] += $count_bases;
      $split[3] += $count_bases;

      my $newline = join " ", @split;

      # print to ace file
      push @{$change->{'ace-delete'}{"Sequence : \"$clone\""}}, "-D $line"; # create output to be written to the ace file
      #print "-D $line\n";
      push @{$change->{'ace-add'}{"Sequence : \"$clone\""}}, "$newline"; # create output to be written to the ace file
      #print "$newline\n";
    }
  }

  return 0;
}

##########################################
# $result = &change_clone_overlaps($change, @lines)
#
# @lines = contents of the current clone, not the prev_clone
#
# overlap_right = name of clone to right of this one, position in this
#                 clone that the next clone's finished sequence starts at
#
# clone_right_end = name of clone to the left of this one, position in
#                   this clone that the previous clone ends at 
#                   - but this is unfinished sequence so we don't want to use this value to find overlaps
#
#                   and/or, name of this clone, length of this clone
#
# clone_left_end = name of clone to right of this one, position in
#                  this clone that the next clone starts at
#                  - but this is unfinished sequence so we don't want to use this value to find overlaps
#
#                   and/or, name of this clone, 1
#
# overlap_left = name of clone to the left of this one
#                - look in this clone to get its length and overlap_right value
#		   so that we can see the amount of overlap with this clone
#
# ------ = finished sequence
# ...... = unfinished sequence
#                 
#                 
#                 
#                 
#------------------------|..........|                                                                         previous clone
#                 
#          |.........|--------------V--------------------------V-------------------V-------|.......|          this clone
#                                   clone                      clone               overlap
#                                   right                      left                right
#                                   end                        end
#                                                              |...................|-----------------------   next clone
#                                                             
#
# If the change occurs in the regions of overlap at the ends, then we
# need to correct both clones
#
# The unfinished sequence region does not need the other clone to be
# changed because the sequence usually doesn't exist - it has been
# clipped off
#
# The region between overlap_right and the end of the clone is the
# overlap of existing clone sequences. Changes within this
# regionrequire both clones to be changed.


sub change_clone_overlaps {

  my ($change, @lines) = @_;

  my $clone = $change->{'clone'};
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};
  my $id = $change->{'ID'};

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }


  my $newline;

  foreach my $line (@lines) {
    chomp $line;

    #print "$line\n";

    # if the change occurs in the overlap regions, we need to raise an error
    # because changes shouldn't occur there - plenty of people will have
    # scrutinised the region and a change here is probably an error
    if ($line =~ /^Overlap_right\s+\S+\s+(\d+)/) {
      my $overlap_right = $1;
      if ($end_pos >=  $overlap_right) {
	$log->write_to( "*** A $change_type change for $id has been attempted in the region $clone:$start_pos..$end_pos overlapping the clone after $clone,\n");
	$log->write_to( "(the region after $clone:$overlap_right)\n");
	$log->write_to( "This is specified by the line: $line\n");
	$log->write_to( "You should investigate this change to see if it is real.\n\n");
	$errors .= "A change for $id has been attempted in the region $clone:$start_pos..$end_pos overlapping the clone after $clone, (the region after $clone:$overlap_right)\n";
	return 1;
	# if we want to deal with this situation, want to see if we should shift up objects on the region overlap_right to the end of the clone
      }

      my @f = split /\s+/, $line;
      if ($f[2] > $end_pos) {
	$f[2] += $count_bases;

	my $newline = join " ", @f;

	# print to ace file
	push @{$change->{'ace-delete'}{"Sequence : \"$clone\""}}, "-D $line"; # create output to be written to the ace file
	#print "$line\n";
	push @{$change->{'ace-add'}{"Sequence : \"$clone\""}}, "$newline"; # create output to be written to the ace file
	#print "$newline\n";
      }

    } elsif ($line =~ /^Overlap_left\s+(\S+)/) { 
      if (lc($1) eq lc("\"$clone\"")) {next;} 		# what follows 'Overlap_left' shouldn't match "$clone" 
      my $prev_clone = $1;
      $prev_clone =~ s/\"//g;		# strip quotes
      # get the length and overlap_right value for this clone to the left
      #print "prev_clone = $prev_clone\n";
      my $prev_clone_obj = $ace->fetch(Sequence => $prev_clone);
      my $prev_overlap_right = $prev_clone_obj->Overlap_right(2);
#      my $prev_overlap_right = $prev_clone_obj->at('Structure.Overlap_right[2]');
      #print "prev_overlap_right = $prev_overlap_right\n";
      my $prev_length = $prev_clone_obj->DNA(2);
#      my $prev_length = $prev_clone_obj->at('DNA[2]');
      #print "prev_length = $prev_length\n";
      my $prev_overlap = $prev_length - $prev_overlap_right + 1;

      if ($start_pos <= $prev_overlap && $end_pos > $prev_overlap) {
	if ($start_pos <= $prev_overlap && lc $change_type eq 'deletion') {
	  $log->write_to( "*** A change for $id has been attempted in the region $clone:$start_pos..$end_pos overlapping the clone before $clone,\n");
	  $log->write_to( "*** But crucially, it also is in the region past the end of the overlap, so some of the change does not affect the previous clone - consider changing the clone boundaries.\n");
	  $log->write_to( "You should investigate this change to see if it is real.\n\n");
	  $errors .= "A change for $id has been attempted in the region $clone:$start_pos..$end_pos overlapping the clone before $clone but also out of this region\n";
	  return 1;
	} elsif ($start_pos == $prev_overlap && (lc $change_type eq 'insertion' || $change_type eq 'SNP' || $change_type eq 'substitution'))  {
	  $log->write_to( "*** A change for $id has been attempted in the region $clone:$start_pos..$end_pos overlapping the clone before $clone,\n");
	  $log->write_to( "*** But crucially, it also is in the region past the end of the overlap, so some of the change does not affect the previous clone - consider changing the clone boundaries.\n");
	  $log->write_to( "You should investigate this change to see if it is real.\n\n");
	  $log->write_to( "It is not a $change_type, the clone end is at $prev_overlap and the change starts at $start_pos\n");
	  $errors .= "A change for $id has been attempted in the region $clone:$start_pos..$end_pos overlapping the clone before $clone but also out of this region\n";	  
	  $errors .= "(It is a $change_type, the clone end is at $prev_overlap and the change starts at $start_pos so I think this is OK)\n";
	  return 0;
	} else {
	  $log->write_to( "*** A change for $id has been attempted in the region $clone:$start_pos..$end_pos overlapping the clone before $clone,\n");
	  $log->write_to( "*** But crucially, it also is in the region past the end of the overlap, so some of the change does not affect the previous clone - consider changing the clone boundaries.\n");
	  $log->write_to( "You should investigate this change to see if it is real. \n\n");
	  $log->write_to( "It is a $change_type, the clone end is at $prev_overlap and the change starts at $start_pos\n");
	  $errors .= "A change for $id has been attempted in the region $clone:$start_pos..$end_pos overlapping the clone before $clone but also out of this region\n";	  
	  return 1;
	}
      }

      if ($end_pos < $prev_overlap) {
	$log->write_to( "*** A change for $id has been attempted in the region $clone:$start_pos..$end_pos overlapping the clone before $clone,\n");
	$log->write_to( "You should investigate this change to see if it is real.\n\n");
	$errors .= "A change for $id has been attempted in the region $clone:$start_pos..$end_pos overlapping the clone before $clone - not an error, but it should be checked carefully.\n";

	if ($change_type eq "shift-overlap") { # this is not finished
	  # print the superlink shift to ace file
	  my $or = $prev_overlap_right - $count_bases;
#	  push @{$change->{'ace-delete'}{"Sequence : \"$prev_clone\""}}, "-D Overlap_right $clone $prev_overlap_right"; 
#	  push @{$change->{'ace-add'}{"Sequence : \"$prev_clone\""}}, "Overlap_right $clone $or"; 
	  print "#######################################################################\n";
	  print "You will probably have to ask St. Louis to add the following to stlace:\n\n";
	  print "Sequence : \"$prev_clone\"\n";
	  print "-D Overlap_right $clone $prev_overlap_right\n\n";
	  print "Sequence : \"$prev_clone\"\n";
	  print "Overlap_right $clone $or\n\n";
	  print "#######################################################################\n";

	} else {
	  # make change to sequence in clone to left
	  &change_clone_to_left($prev_clone, $prev_overlap_right, $prev_length, $prev_overlap, $change);
	}
      }


    } elsif ($line =~ /^Clone_right_end\s+(\S+)\s+(\d+)/) {

      my $prev_clone = $1;
      my $clone_right_end = $2;
      if (! defined $prev_clone || ! defined $clone_right_end) {next;}

      if ($prev_clone eq $clone) {
	$clone_right_end += $count_bases; # update length of this clone

	# print to ace file
	push @{$change->{'ace-delete'}{"Sequence : \"$clone\""}}, "-D $line"; 
	push @{$change->{'ace-add'}{"Sequence : \"$clone\""}}, "Clone_right_end $prev_clone $clone_right_end"; 


      } else {			# we have the clone to the left
	if ($end_pos <= $clone_right_end) {
	  $clone_right_end += $count_bases; # update the end of the previous clone ni this clone if the change region is in the overlap

	  # print to ace file
	  push @{$change->{'ace-delete'}{"Sequence : \"$clone\""}}, "-D $line"; 
	  push @{$change->{'ace-add'}{"Sequence : \"$clone\""}}, "Clone_right_end $prev_clone $clone_right_end"; 

	}
      }


    } elsif ($line =~ /^Clone_left_end\s+(\S+)\s+(\d+)/) {

      my $next_clone = $1;
      my $clone_left_end = $2;
      if (! defined $next_clone || ! defined $clone_left_end) {next;}

      if ($next_clone ne $clone) {
	if ($end_pos <= $clone_left_end) {
	  $clone_left_end += $count_bases; # update the end of the previous clone in this clone if the change region is in the overlap

	  # print to ace file
	  push @{$change->{'ace-delete'}{"Sequence : \"$clone\""}}, "-D $line"; 
	  push @{$change->{'ace-add'}{"Sequence : \"$clone\""}}, "Clone_left_end $next_clone $clone_left_end"; 

	}
      }
    }
  }

  return 0;
}
##########################################
# make change to sequence in clone to left
# this is for when the change location is in the overlap region between
# this clone and the previous one

sub change_clone_to_left {
  my ($prev_clone, $prev_overlap_right, $prev_length, $overlap, $change) = @_;


  # change the sequence of the previous clone
  my $clone_obj = $ace->fetch(Sequence => $prev_clone);
  my $ID = $change->{'ID'};
  my $type = $change->{'change_type'};
  print "Making change $type $ID to previous clone $prev_clone\n\n";


  #print "reading DNA sequence\n";
  my $dna = $clone_obj->asDNA();
  $dna =~ s/\>(\w+)\n//;	# remove title line
  $dna =~ s/\n//g;		# make into one line
  #print "dna = $dna\n";


  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};
  ## if we are using the positions-file input, then count_bases needs
  ## to be changed from the length of the changed region to the amount
  ## the clone sequence changes
  #if ($positionsfile) {
  #  if (lc $change_type eq 'deletion') {
  #    $count_bases = -$count_bases;
  #  }
  #}

  my $from_base = $change->{'from_base'};
  my $to_base = $change->{'to_base'};
  my $chrom = $change->{'chrom'};
  my $start = $change->{'start'};

  # get the position of the change in the clone to the left
  my $prev_start_pos = $prev_overlap_right + $start_pos - 1;
  my $prev_end_pos = $prev_overlap_right + $end_pos - 1;
  if ($prev_end_pos > $prev_length) {
    die "In change_clone_to_left() doing $ID we have the end of the change region past the end of the previous clone\n";
  }

  # error checks
  if ($prev_start_pos == $prev_length && lc($change_type) eq 'insertion') {
    die "$ID has Insertion right at end of previous clone $prev_clone\n";
  }
  if ($prev_start_pos + $count_bases - 1 > $prev_length && lc($change_type) eq 'deletion') {
    die "$ID has Deletion past end of previous clone $prev_clone\n";
  }
  




  ######################
  # DELETION
  ######################

  if (lc($change_type) eq 'deletion') {
    my $clone_base = substr($dna, $prev_start_pos-1, $count_bases);
    if (lc $clone_base ne lc $from_base) {die "In $ID deletion in previous overlapping clone $prev_clone of $chrom $start bases $from_base to $to_base ($count_bases bases changed) : there is a mismatch to the clone base found: $clone_base\n"}
    substr($dna, $prev_start_pos-1, $count_bases) = '';
  }

  ######################
  # INSERTION
  ######################

  elsif (lc($change_type) eq 'insertion') {
    substr($dna, $prev_start_pos, 0) = lc($to_base);
  }
  
  ######################
  # SUBSTITUTION
  ######################
  
  elsif (lc($change_type) eq 'substitution' || lc($change_type) eq 'snp') {
    my $clone_base = substr($dna, $prev_start_pos-1, length $change->{'from_base'});
    if (lc $clone_base ne lc $from_base) {die "In $ID substitution in previous overlapping clone $prev_clone of $chrom $start bases $from_base to $to_base ($count_bases bases changed) : there is a mismatch to the clone base found: $clone_base\n"}
    substr($dna, $prev_start_pos-1, length $change->{'from_base'}) = lc($to_base);
  }

  else {
    die "Unknown type of change: $change_type in $ID\n";
  }
  
  
  # reformat the DNA sequence into lines of 60 characters
  $dna = &reformat($dna);
  
  # write the dna back to the database
  push @{$change->{'ace-add'}{"DNA : \"$prev_clone\""}}, "$dna"; # create output to be written to the ace file

  # shift up objects in the clone to the left

}


##########################################
# add a Remark - definitely want to update this
sub add_remark {

  my ($change, @lines) = @_;

  my $clone = $change->{'clone'};
  my $region = $change->{'region'};
  if (!defined $region) {$region = ''}
  my $change_type = $change->{'change_type'};
  my $start_pos = $change->{'start_pos'};
  my $end_pos = $change->{'end_pos'};
  my $count_bases = $change->{'count_bases'};
  my $from_base = $change->{'from_base'};
  my $to_base = $change->{'to_base'};
  my $ID = $change->{'ID'};

  if (!defined $ID) {$ID = ''}

  # if we are using the positions-file input, then count_bases needs
  # to be changed from the length of the changed region to the amount
  # the clone sequence changes
  if ($positionsfile) {
    if (lc $change_type eq 'deletion') {
      $count_bases = -$count_bases;
    }
  }


  #get date for remark
  my ($day, $mon, $yr)  = (localtime)[3,4,5];
  my $date = sprintf("%02d%02d%02d",$yr-100, $mon+1, $day);
  my $remark = "DB_remark \"[$date] Sequence correction $ID : $change_type ";
  if ($from_base ne '' && $to_base ne '') {$remark .= "$from_base to $to_base"}
  if ($from_base eq '') {$remark .= "$to_base"}
  if ($to_base eq '') {$remark .= "$from_base"}

  if ($count_bases > 0) {
    if ($count_bases < 20) {
      $remark .= "$count_bases bases from $region @ $start_pos";
    } else {
      $remark .= "$count_bases bases @ $start_pos";
    }
  } else {
    if ($count_bases > -20) {
      $remark .= "$count_bases bases $region @ $start_pos";
    } else {
      $remark .= "$count_bases bases @ $start_pos";
    }    
  }


  # print to ace file
  push @{$change->{'ace-add'}{"Sequence : \"$clone\""}}, $remark; # create output to be written to the ace file

  return 0;
}

##########################################
# change the stored start and end positions of the changed region from
# the clone coords to the superlink coords

#  &use_superlink_coords($change, $superlink, @slurp);

sub use_superlink_coords {

  my ($change, $superlink, @lines) = @_;

  my $clone = $change->{'clone'};

  # get the position of the end of the last clone on this superlink
  # so that we can get the length of the superlink
  $change->{superlink_length} = 0;
  foreach my $line (grep /Subsequence/, @lines) {
    if ($line =~ /Subsequence\s+\S+\s+\d+\s+(\d+)/) {
      if ($1 > $change->{superlink_length} ) {
	$change->{superlink_length} = $1; 
      }
    }
  }
  
  foreach my $line (grep /Subsequence/, @lines) {
    chomp $line;

    #print "$line\n";

    my @split = split /\s+/, $line; 
    # if the start or end position is not present, skip it
    if (! defined $split[1] || ! defined $split[2] || ! defined $split[3]) {next;} 
    if ($split[1] eq "\"$clone\"") {
      # now we have the start position of the clone on the superlink in $split[2]
      # use this to change the position of the changed region to superlink coords
      $change->{'start_pos_on_clone'} = $change->{'start_pos'};
      $change->{'end_pos_on_clone'} = $change->{'end_pos'};
      $change->{'start_pos'} = $change->{'start_pos'} + $split[2] - 1;
      $change->{'end_pos'} = $change->{'end_pos'} + $split[2] - 1;

      # store the start and end position of the clone on the superlink
      $change->{'clone_start_on_superlink'} = $split[2];
      $change->{'clone_end_on_superlink'} = $split[3];

      # store the clone name in case we need it in future
      $change->{'original_clone_name'} = $change->{'clone'};

      # change the name of the clone to that of the superlink for ease of
      # using the subroutines that expect a clone name
      $change->{'clone'} = $superlink;

      # and note that we are using the superlink
      $change->{'superlink'} = 1;

      return 0;
    }
  }

  $log->write_to( "*** The position of the clone $clone was not found in the superlink $superlink.\n");
  $log->write_to( "The coordinates of objects in the superlink cannot therefore be adjusted.\n");
  return 1;
}


##########################################
# change the stored start and end positions of the changed region from
# the superlink coords to the chrom coords

#  &use_chrom_coords($change, $superlink, $chrom, @slurp);

sub use_chrom_coords {

  my ($change, $superlink, $chrom, @lines) = @_;

  my $clone = $change->{'clone'};
  my $ID = $change->{'ID'};
  print "looking for $clone with Feature $ID\n";

  # get the position of the end of the last superlink on this chrom
  # so that we can get the length of the chrom
  $change->{chrom_length} = 0;
  foreach my $line (grep /Subsequence/, @lines) {
    if ($line =~ /Subsequence\s+\S+\s+\d+\s+(\d+)/) {
      if ($1 > $change->{chrom_length} ) {
	$change->{chrom_length} = $1; 
      }
    }
  }
  
  foreach my $line (grep /Subsequence/, @lines) {
    chomp $line;

    #print "$line\n";

    my @split = split /\s+/, $line; 
    # if the start or end position is not present, skip it
    if (! defined $split[1] || ! defined $split[2] || ! defined $split[3]) {next;} 
    if ($split[1] eq "\"$clone\"") {
      # now we have the start position of the superlink on the chrom in $split[2]
      # use this to change the position of the changed region to chrom coords
      $change->{'start_pos_on_superlink'} = $change->{'start_pos'};
      $change->{'end_pos_on_superlink'} = $change->{'end_pos'};
      $change->{'start_pos'} = $change->{'start_pos'} + $split[2] - 1;
      $change->{'end_pos'} = $change->{'end_pos'} + $split[2] - 1;


      $change->{'superlink_start_on_chrom'} = $split[2];
      $change->{'superlink_end_on_chrom'} = $split[3];

      # change the name of the clone to that of the chrom for ease of
      # using the subroutines that expect a clone name
      $change->{'clone'} = $chrom;

      # and note that we are using the chrom
      $change->{'chrom'} = 1;

      return 0;
    }
  }

  $log->write_to( "*** The position of the superlink $superlink was not found in the chrom $chrom.\n");
  $log->write_to( "The coordinates of objects in the chrom cannot therefore be adjusted.\n");
  return 1;
}



##########################################
# returns the reverse complement of a DNA sequence
# preserving the case

sub DNA_string_reverse {
  my $revseq = reverse shift;
  $revseq =~ tr/a/x/;
  $revseq =~ tr/A/X/;
  $revseq =~ tr/t/a/;
  $revseq =~ tr/T/A/;
  $revseq =~ tr/x/t/;
  $revseq =~ tr/X/T/;
  $revseq =~ tr/g/x/;
  $revseq =~ tr/G/X/;
  $revseq =~ tr/c/g/;
  $revseq =~ tr/C/G/;
  $revseq =~ tr/x/c/;
  $revseq =~ tr/X/C/;
  return ($revseq);
}

##########################################
# change the data on the clone object

sub write_clone_sequence_date {
  my ($clone, $clonedates) = @_;

  my $dat = `date +%y%m%d`;
  chomp $dat;

# write the ace file to change the clone sequence date to match the clone's sequence directory date
  open (CD, ">> $clonedates") || die "cant open the clone date ace file '$clonedates'\n";
  print CD "\nSequence $clone\n";
  print CD "Date_directory $dat\n\n";
  close (CD);

}

##########################################
# get the composition of a sequence

sub composition {
  my ($dna) = @_;
  $dna =~ s/\n//g;		# make into one line

  my ($a, $c, $g, $t, $n, $length_dna);

  $a = $dna =~ tr/[aA]/A/;
  $c = $dna =~ tr/[cC]/C/;
  $g = $dna =~ tr/[gG]/G/;
  $t = $dna =~ tr/[tT]/T/;

  # the Ns are whatever is not ACGT
  $length_dna = $n = length $dna;
  $n -= $a;
  $n -= $c;
  $n -= $g;
  $n -= $t;

  return ($a, $c, $g, $t, $n, $length_dna);

}

##########################################
# get the previous date folder from the current.versions file
# update the 'current.versions' file to point to the new directory
# my $prev_dat = &get_prev_dat($clone, $dat, $noload);

sub get_prev_dat {
  my ($clone, $dat, $noload) = @_;

  my $prev_dat;
  my $file = "/lustre/scratch101/ensembl/wormpipe/wormpub/analysis/cosmids/current.versions";

  open (INPUT, "< $file") || die "Can't open file $file\n";
  open (NEW, "> $file.new") || die "Can't open file $file.new\n";
  while (my $line = <INPUT>) {
    chomp $line;
    my ($next_clone, $next_dat) = ($line =~ /(\S+)\/(\d+)/ );
    if ($next_clone eq $clone) {
      $prev_dat = $next_dat;
      print NEW "$clone/$dat\n";
    } else {
      print NEW "$line\n";
    }
  }
  close(NEW);
  close(INPUT);
 
  if ($noload) {
    print "*** -noload specified: NOT MOVING $file.new TO BE $file\n"; 
    system("mv $file.new $file.test");
  } else {
    my $time = `date +%y%m%d_%H:%m:%S`;
    chomp $time;
    my $olddir = "/lustre/scratch101/ensembl/wormpipe/wormpub/analysis/cosmids/old_current.versions/";
    system("mv $file $olddir/current.versions.$time");
    system("mv $file.new $file");
  }

  return $prev_dat;
}


##########################################
# this looks at the sequence of the superlinks and checks that their sequence is OK

sub check_superlink_sequence {
  my $result = 1;		# OK so far


  my @superlinks;
##  if ($species eq 'elegans' && $stlace) {
##    @superlinks = qw(SUPERLINK_RW1 SUPERLINK_RW1R SUPERLINK_RW2 SUPERLINK_RW2R SUPERLINK_RW3A SUPERLINK_RW3B SUPERLINK_RW4 SUPERLINK_RW5 SUPERLINK_RWXL SUPERLINK_RWXR);
##  } elsif ($species eq 'elegans' || $camace) {
##    @superlinks = qw(SUPERLINK_CB_I SUPERLINK_CB_II SUPERLINK_CB_IIIL SUPERLINK_CB_IIIR SUPERLINK_CB_IR SUPERLINK_CB_IV SUPERLINK_CB_V SUPERLINK_CB_X);
##  }

  @superlinks = qw(SUPERLINK_RW1 SUPERLINK_RW1R SUPERLINK_RW2 SUPERLINK_RW2R SUPERLINK_RW3A SUPERLINK_RW3B SUPERLINK_RW4 SUPERLINK_RW5 SUPERLINK_RWXL SUPERLINK_RWXR SUPERLINK_CB_I SUPERLINK_CB_II SUPERLINK_CB_IIIL SUPERLINK_CB_IIIR SUPERLINK_CB_IR SUPERLINK_CB_IV SUPERLINK_CB_V SUPERLINK_CB_X);
#  @superlinks = qw(SUPERLINK_RW5 SUPERLINK_CB_V);

  print "Check superlink composition:\n";
  print "Superlink\t\tA\tC\tG\tT\tN\tLength\n";

  foreach my $superlink (@superlinks) {
    my $super_obj = $ace->fetch(Sequence => $superlink);
    my $dna = $super_obj->asDNA();
    $dna =~ s/\>(\w+)\n//;	# remove title line
    $dna =~ s/\n//g;	        # remove newline characters
    my ($a, $c, $g, $t, $n, $length_dna) = &composition($dna);
    print "$superlink\t\t$a\t$c\t$g\t$t\t$n\t$length_dna\n";
    if ($n != 0) {
      print "ERROR There are non-ACGT bases in the superlink $superlink!\n";
      $result = 0;
    }
  }

  my @chromosomes = qw(CHROMOSOME_I CHROMOSOME_II CHROMOSOME_III CHROMOSOME_IV CHROMOSOME_V CHROMOSOME_X);
#  my @chromosomes = qw(CHROMOSOME_V);

  print "\nCheck chromosome composition:\n";
  print "Chromosome\t\tA\tC\tG\tT\tN\tLength\n";

  foreach my $chromosome (@chromosomes) {
    my $chr_obj = $ace->fetch(Sequence => $chromosome);
    my $dna = $chr_obj->asDNA();
    $dna =~ s/\>(\w+)\n//;	# remove title line
    $dna =~ s/\n//g;	        # remove newline characters
    my ($a, $c, $g, $t, $n, $length_dna) = &composition($dna);
    print "$chromosome\t\t$a\t$c\t$g\t$t\t$n\t$length_dna\n";
    if ($n != 0) {
      print "ERROR There are non-ACGT bases in the superlink $chromosome!\n";
      $result = 0;
    }
  }

  return $result;
}

##########################################
# reformat a DNA string into lines of 60 characters

sub reformat {
    my $in_string = shift;
    my $out_string = "";

    my $string_len = length ($in_string);
    my $lines = int ($string_len / 60) ;

    for (my $i = 0; $i <= $lines; $i++) {
        $out_string = $out_string . substr($in_string,($i*60),60) . "\n";
    }
    return ($out_string);
}


##########################################

sub usage {
  my $error = shift;

  if ($error eq "Help") {
    # Normal help menu
    system ('perldoc',$0);
    exit (0);
  }
}

##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################
##########################################






# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - genome-changes.pl

=head1 USAGE

=over 4

=item genome-changes.pl  [-options]

=back

This script reads a data file specifying changes to be made to the
genomic sequence of clones in camace and shifts up objects that are
affected by the change accordingly.

The format of the input data file is three columns.
The first column is the name of the clone.
The second column is the type of change, either 'insertion' or 'deletion' or 'shift-overlap'.
The third column is the lowercase flanking sequences of the bases to be changed with the bases to be changed in uppercase.
Comments lines start with '#'.
Blank lines are ignored.

An example of the input file is:
                                
C34B4    Insertion    ttgagtttgatggttcaactgaaatTggtcagtgtc
C34B4    Insertion    ggtcagtgtcTttcttcactttgcctgaaacttgga
F21G4    Insertion    accattgggaattcccggagGaaaagtgtgatgttttctttaaat

A subsitution is simply a deletion followed by an insertion where the
clone name and the flanking sequences are the same:

# this is a substitution of 'G' by 'ACT'
AC3      Deletion     agtcgagtcgtagtcgtGgatcggtagcgatgcgtgtgtt
AC3      Insertion    agtcgagtcgtagtcgtACTgatcggtagcgatgcgtgtgtt


A 'shift-overlap' change does not change the genomic sequence.  It
changes the boundary of a superlink by changing the position of the
start of the first clone on this superlink by changing the
overlap_right value in the last clone of the previous superlink.  This
previous superlink will probably belong to St. Louis, so you will be
told to ask them to change their clone's overlap_right value.  Because
your superlink and clone now have some extra sequence at their start,
they will need everything in them shifted up as if there had been a
genomic insertion.

An example of a 'shift-overlap' change input file is:

Y95D11A shift-overlap   TTATATATTTTTTTGGAAATTTATAACTCTTAAAAAAATTCAATTTTTTCAAATAAATAAAATTTCAGATGGCTTCTCAACCGGAGCTCATAATGGTTGACGAGCAAGTCGTCGCTTATGAAGTAGAAATTGATAGTTTTGATGTAAAATATGATGAAGAGGAACATGATGGTCAAGGGACACAAGATGAACCATTTTCTCATGGTACGGAACAGTTTTACGCTGAAAAATTCCAGAATTCCAAAAAATGAAACCTAAAATAGTGATAAAAAGGCGTTTTGAATATTAAATTGAAGAAAAAAATCAGCAAAAATTGTTCAAAATCAAGAATTTTAACGGAAAAGTGTAAAATCTTCTCCACGGGGAGTACACATGCTTCGTAAATCGACATATGGTCAATTTTAAAGTTTTGAAAATTGAAATGCCGGCAAAAAATCTTTTCTTGTTTTTTTTTCGCAAAAAATTCAATTTTCGAAAAAATAATTATAGAAAATTGCATTTTTTGACCGAAAAGTCAATAAAAATAACAGAAAAAATCGATAAACCGTTGAAAAATTTTTTTTTAATTCAAAAATTCAGAAATTCTTAAAATTCAAATTTCCAGATGAGCCAAGCACCAGCGGTTATCACCATCACTACCAATTTCCCAATGACGTGGATCCAAATGATGTTTATTTATTCGATGAGGTATCAATTATCCGAAATTTGGCGATTTTTGAGCCAAAACTACGGTACCCGGTCTCGACACGACAATTTTTGTTAAATTAAAAAAGGTGTGCGCCTTTGAAGGTTACTGTAGTTTCGAACTTTTGCTGATTTTTCATATTTTTTCGTTGAAAACAAAAGTATTTATTTGTTGAAAATCAGAAAATATTATCTTCGCGTCGAGACCTATTACCATTCTATTTTTGCCGCAAAAAACAAAATTTCCTTTAAAAAAAAGCTAATTTTTCCAAGTTTTTCCAGGAAACTGATCAAATTCATCAGCTCGACCCGAATCAACTCAAAAATAATGAAGAAATTGACGATGTCGAATATATTGATCAATCTGTGCCTTCCACGTCATCAATGATGACGTCACTGCCGTCAACGGTGGCTCCAGTTCAGCCAAATACGTATTACAGACGGAAATCTGGAGGCCCAACTGCAACTGGAAATGAAAAACCGAATTATAGGCCGTTGGCGTTCCAAACGGTTCGtaaaataaaaaaaatgtccatgtgtcgatt

(That is all on one line. The 'taaaataaaaaaaatgtccatgtgtcgatt' is the
old start of your clone with the sequence before it being the sequence
on the previous clone on the previous superlink that we wish to
include as part of your superlink and clone.)



An alternative input file for when you have lots of changes to be made
which may have overlapping flanking sequenecs is the -positions file
input.

This has the format:

# chromosome    start   end     feature_id      insertion_deletion      from_base       to_base mcgrath weber   hobart  RNASEQ_REF      RNASEQ_NONREF   NONREF_PERCENT  AUTO-ACCEPT/REJECT
      CONSEQUENCES
CHROMOSOME_I   5042625 5042627 WBsf898113      Deletion        GGT     -       0       1       0       -       -       -       Reject  C46H11.7:inframe_deletion
CHROMOSOME_I   10113280        10113280        WBsf268456      Deletion        G       -       1       1       0       -       -       -       Accept  Y106G6D.3:coding_sequence_variant
CHROMOSOME_I   10849084        10849087        WBsf898091      Deletion        CGAT    -       0       1       0       -       -       -       Reject  C35E7.1a:frameshift_variant
CHROMOSOME_I    13101700        13101700        WBsf268374      Deletion        A       -       1       1       0       -       -       -       Accept  Y26D4A.21:frameshift_variant
CHROMOSOME_I    13101705        13101705        WBsf268374      Deletion        T       -       1       1       0       -       -       -       Accept  Y26D4A.21:frameshift_variant
CHROMOSOME_I    13101710        13101710        WBsf268374      Deletion        A       -       1       1       0       -       -       -       Accept  Y26D4A.21:frameshift_variant
CHROMOSOME_I   11145687        11145688        WBsf267982      Insertion       -       G       1       1       0       -       -       -       Accept  ZK39.9:splice_donor_variant,frameshift_variant
CHROMOSOME_I   13617731        13617732        WBsf898399      Insertion       -       T       0       1       0       -       -       -       Reject  Y6B3A.1c:frameshift_variant
CHROMOSOME_I   14349945        14349946        WBsf267937      Insertion       -       C       1       1       0       -       -       -       Accept  Y105E8A.2:splice_donor_variant
CHROMOSOME_I   9457732 9457732 WBsf899197      SNP     G       A       1       1       1       0       8       100.0   Accept  
CHROMOSOME_I   9737907 9737907 WBsf898062      SNP     C       A       1       1       1       0       0       -       Accept  
CHROMOSOME_I   9940581 9940581 WBsf898914      SNP     A       T       1       1       1       0       4       100.0   Accept  

If the lines starts with '#' or does not start with 'CHROMOSOME' or 'chr' then it is ignored

The columns are then:

Chromosome : e.g. CHROMOSOME_I - the sequence to change
start : e.g. 232039 - the start position or the position to nisert to the right of
end : e.g. 232039 - the end position
ID : usually the Feature that marks the site
type of change : Insertion, Deletion, SNP (i.e. substitution)
Base(s) to change from : or '-' if insertion
Base(s) to change to : or '-' if deletion

any thing else on the line will be ignored.





script_template.pl MANDATORY arguments:

=over 4

=item -database, path to camace - this is the database that will be copied

=back

=over 4

=item -new_database, path to where camace should be copied - this copy of the database will be changed

=back

=over 4

=item -infile, path to file specifying the changes. Three columns: clone name, type of change ('insertion' or 'deletion'), flanking sequence of change with changed bases in uppercase.

=back

=over 4

=item -position, path to file specifying the changes. This is an alternative input file to -infile.

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -h, Help

=back

=over 4
 
=item -debug, Debug mode, set this to the username who should receive the emailed log messages. The default is that everyone in the group receives them.
 
=back

=over 4

=item -test, Test mode, run the script, but don't change anything.

=back

=over 4
    
=item -verbose, output lots of chatty test messages

=back


=over 4
    
=item -nomove, do not move the -database to the -new_database, assume it is in place already

=back


=over 4
    
=item -noload, do not load the ACE files into the copy of the database - for testing purposes

=back

=over 4
    
=item -ignore_assembly_tags - don't abort changes made in the region of an assembly tag, just report these changes

=back

=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Gary Williams (gw3@sanger.ac.uk)

=back

=cut
