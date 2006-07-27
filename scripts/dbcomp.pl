#!/usr/local/bin/perl5.8.0 -w
#
# dbcomp.pl
#
# Counts the number of objects in an ACEDB database for each Class stated in the config file
# Compares this number to those from a second database.
#
# Last updated by: $Author: mh6 $
# Last updated on: $Date: 2006-07-27 08:59:25 $


use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use IO::Handle;


$|=1;


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($database, $database2, $all, $midway, $wee);
our ($errfile,$outfile); 
our ($db_1, $db_2, $dbname_1, $dbname_2);

GetOptions (
	    "help"          => \$help,
            "debug=s"       => \$debug,
	    "test"          => \$test,
	    "verbose"       => \$verbose,
	    "store:s"         => \$store,
	    "database=s"    => \$database,
	    "database2=s"   => \$database2,
	    "all"           => \$all,
	    "wee"           => \$wee,
	    "midway"        => \$midway,
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

# in test mode?
if ($test) {
  print "In test mode\n" if ($verbose);

}

# establish log file.
my $log = Log_files->make_build_log($wormbase);


############################################
# Get paths and versions
############################################

my $WS_current  = $wormbase->get_wormbase_version;
my $WS_previous = $WS_current - 1;
my $exec        = $wormbase->tace;

#################################################################
# Compare autoace to previous build unless -database specified
#################################################################

$dbname_1    = "WS${WS_previous}";
$db_1        = $wormbase->database("WS${WS_previous}");

$dbname_2    = "WS${WS_current}";
$db_2        = $wormbase->autoace;

# First alternative database specified?
if ($database) {
    $dbname_1  = "$database";
    $db_1      = "$database"; 
}

# Second alternative database specified?
if ($database2) {
    $dbname_2  = "$database2";
    $db_2      = "$database2"; 
}

# create_log_files;

$log->write_to("\n");
$log->write_to("Previous db      : $dbname_1 '$db_1'\n");
$log->write_to("Current db       : $dbname_2 '$db_2'\n");
$log->write_to("\n\n");
				
# open two main output files to store results
$errfile = $wormbase->compare."/WS${WS_previous}-WS${WS_current}.out";
$outfile = $wormbase->compare."/WS${WS_previous}-WS${WS_current}.dbcomp";

open (OUT, ">$outfile") || die "Couldn't write to out file\n";
open (ERR, ">$errfile") || die "Couldn't write to err file\n";
  
print OUT  " +--------------------------------------------+\n";
print OUT  " | Class                  |   ACEDB database  |\n";
print OUT  " |                        +---------+---------+---------+---------+---------+\n";
printf OUT " |                        | %7s | %7s |    +    |    -    |   Net   |\n", $dbname_1,$dbname_2;
print OUT  " +------------------------+---------+---------+---------+---------+---------+\n";

#########################################################################
# Read list of classes
# The complete set of classes to dump is between the DATA and END tokens
#########################################################################

my @TotalClasses;


#@TotalClasses = &full_run if ($all);
#@TotalClasses = &mid_run if ($midway);


# run midway list if asked too, else default to the full list
if ($midway) {
    @TotalClasses = &mid_run;
}
elsif ($wee) {
    @TotalClasses = &wee_run;
}
else {
    @TotalClasses = &full_run;
}



#######################
# Main loop
#######################

my $counter = 1; # for indexing each class to be counted

foreach my $query (@TotalClasses) {
  next if ($query eq "");

  $log->write_to(" Counting '$query'\n");
  printf OUT " | %22s |", $query;

  ##################################################
  # Get class counts from both databases           #
  ################################################## 
  
  my ($class_count_1,$class_count_2) = &count_class($query,$counter);

  $log->write_to(" Counting $dbname_1 : $class_count_1\n");  
  $log->write_to(" Counting $dbname_2 : $class_count_2\n\n");
  
  printf OUT " %7s ",$class_count_1;
  print OUT "|";  
  printf OUT " %7s ",$class_count_2;


  ##################################################
  # Calculate difference between databases         #
  ##################################################
  
  my ($added, $removed) = &diff($counter);
  my $diff = $added - $removed;

#  printf OUT "| %6s |\n",$diff;

  printf OUT "| %7s ", $added;
  printf OUT "| %7s ", $removed;
  printf OUT "| %7s |\n",$diff;

  $counter++;
}


print OUT  " +------------------------+---------+---------+---------+---------+---------+\n";





# create symbolic links to current.out and current.dbcomp
$log->write_to("\nCreating 'current.out' and 'current.dbcomp'\n",$log);
$wormbase->run_command("rm -f ".$wormbase->compare."/current.out",$log);
$wormbase->run_command("rm -f ".$wormbase->compare."/current.dbcomp",$log);
$wormbase->run_command("ln -s $errfile ".$wormbase->compare."/current.out",$log);
$wormbase->run_command("ln -s $outfile ".$wormbase->compare."/current.dbcomp",$log);


close (OUT);
close (ERR);

# write to the release letter - subroutine in Wormbase.pm
$wormbase->release_databases;

# Email log file
$log->mail();

print "Finished.\n" if ($verbose);
exit (0);




##############################################################
#
# Subroutines
#
##############################################################



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

sub count_class{

    my $query  = shift;
    my $counter = shift;
    my $out;
    my $class_count1;
    my $class_count2;
    
    # Formulate query
    my $command = "query find '$query'\nlist -a\nquit\n";
    
    ####################################
    # Count objects in first database
    ####################################
    
    # open temp output file
    $out = "/tmp/dbcomp_A_${counter}";
    open (COUNT, ">$out") || die "Couldn't write to tmp file: $out\n";
    
    # open tace connection and count how many objects in that class
    open (TACE, "echo '$command' | $exec $db_1 | ");
    while (<TACE>) {
	($class_count1 = $1) if (/^\/\/ (\d+) Active Objects/);
	(print COUNT "$_")   if (/\:/); # Add list of object to temp file
    }
    close (TACE);
    close (COUNT);
    
    
    ####################################
    # Count objects in second database
    ####################################
    
    $out = "/tmp/dbcomp_B_${counter}";
    open (COUNT, ">$out") || die "Couldn't write to tmp file: $out\n";
    
    # open tace connection and count how many objects in that class
    open (TACE, "echo '$command' | $exec $db_2 | ");
    while (<TACE>) {
	($class_count2 = $1) if (/^\/\/ (\d+) Active Objects/);
	(print COUNT "$_")   if (/\:/); # Add list of object to temp file
    }
    close (TACE);
    close (COUNT);
    
    return ($class_count1, $class_count2);
    
}


##################################################

sub diff {

  # look at the differences between the two sets of tmp files (one pair of files
  # for each class being queried)
  my $counter = shift;

  my $added   = 0; 
  my $removed = 0;

  system ("cat /tmp/dbcomp_A_${counter} | sort > /tmp/look-1");
  system ("cat /tmp/dbcomp_B_${counter} | sort > /tmp/look-2");
  open (COMM, "comm -3 /tmp/look-1 /tmp/look-2 |");
  while (<COMM>) {
      if (/^(\S+.+)/){
	  print ERR " <- $dbname_1 $1\n";
	  $removed++;
      }
      if (/^\s+(\S+.+)/){
	  print ERR " -> $dbname_2 $1\n";
	  $added++;
      }
  }   
  close (COMM);

  # Tidy up after yourself
  system ("rm -f /tmp/look-1");
  system ("rm -f /tmp/look-2");

  system ("rm -f /tmp/dbcomp_A_${counter}");
  system ("rm -f /tmp/dbcomp_B_${counter}");
 
  # Add break symbol to output file to separate classes
  print ERR "\/\/ -----------------------------\n";
  return($added, $removed);
}

###############################################

sub full_run {

    my @classes = (
		   "2_point_data",
		   "Accession_number",
		   "Anatomy_name",
		   "Anatomy_term",
		   "Antibody",
		   "Author",
		   "briggsae_CDS",
		   "briggsae_genomic",
		   "cDNA_sequence",
		   "CDS",
		   "Cell",
		   "Cell_group",
		   "Class",
		   "Clone",
		   "Coding_transcripts",
		   "Comment",
		   "Condition",
		   "Contig",
		   "Database",
		   "Display",
		   "DNA",
		   "elegans_CDS",
		   "elegans_pseudogenes",
		   "elegans_RNA_genes",
		   "Expression_cluster",
		   "Expr_pattern",
		   "Expr_profile",
		   "Feature",
		   "Feature_data",
		   "Gene",
		   "Gene_class",
		   "Gene_name",
		   "Gene_regulation",
		   "Genome_Sequence",
		   "GO_code",
		   "GO_term",
		   "Homology_group",
		   "Homol_data",
		   "Interaction",
		   "Journal",
		   "Keyword",
		   "Laboratory",
		   "Life_stage",
		   "Lineage",
		   "Locus",
		   "LongText",
		   "Map",
		   "Method",
		   "Microarray_experiment",
		   "Microarray_results",
		   "Model",
		   "Motif",
		   "Movie",
		   "Multi_pt_data",
		   "NDB_Sequence",
		   "nematode_ESTs",
		   "Oligo",
		   "Oligo_set",
		   "Operon",
		   "Paper",
		   "Paper_name",
		   "PCR_product",
		   "Peptide",
		   "Person",
		   "Person_name",
		   "Phenotype",
		   "Picture",
		   "Pos_neg_data",
		   "Protein",
		   "Pseudogene",
		   "Rearrangement",
		   "RNAi",
		   "SAGE_tag",
		   "SAGE_transcript",
		   "SAGE_experiment",
		   "Sequence",
		   "SK_map",
		   "SO_term",
		   "Species",
		   "Strain",
		   "Structure_data",
		   "Table",
		   "Transcript",
		   "Transgene",
		   "Transposon",
		   "Transposon_CDS",
		   "Transposon_family",
		   "Variation",
		   "Y2H"
		   );

    return (@classes);
    
}

sub mid_run {

    my @classes = (
		   "2_point_data",
		   "3d_data",
		   "Anatomy_name",
		   "Anatomy_term",
		   "Antibody",
		   "Author",
		   "briggsae_CDS",
		   "briggsae_genomic",
		   "cDNA_sequence",
		   "CDS",
		   "Cell",
		   "Cell_group",
		   "Class",
		   "Clone",
		   "Cluster",
		   "Coding_transcripts",
		   "Comment",
		   "Condition",
		   "Contig",
		   "Database",
		   "Display",
		   "DNA",
		   "elegans_CDS",
		   "elegans_pseudogenes",
		   "elegans_RNA_genes",
		   "Expr_pattern",
		   "Expr_profile",
		   "Feature",
		   "Feature_data",
		   "Gene",
		   "Gene_class",
		   "Gene_name",
		   "Gene_regulation",
		   "Genome_Sequence",
		   "GO_term",
		   "Homology_group",
		   "Homol_data",
		   "Interaction",
		   "Journal",
		   "Keyword",
		   "Laboratory",
		   "Life_stage",
		   "Lineage",
		   "Locus",
		   "LongText",
		   "Map",
		   "Method",
		   "Microarray_experiment",
		   "Microarray_results",
		   "Model",
		   "Motif",
		   "Movie",
		   "Multi_pt_data",
		   "NDB_Sequence",
		   "nematode_ESTs",
		   "Oligo",
		   "Oligo_set",
		   "Paper",
		   "Paper_name",
		   "PCR_product",
		   "Person",
		   "Person_name",
		   "Phenotype",
		   "Picture",
		   "Pos_neg_data",
		   "Pseudogene",
		   "Rearrangement",
		   "RNAi",
		   "SAGE_tag",
		   "SAGE_transcript",
		   "SAGE_experiment",
		   "Sequence",
		   "SK_map",
		   "SO_term",
		   "Strain",
		   "Table",
		   "Transcript",
		   "Transgene",
		   "Transposon",
		   "Variation",
		   "Y2H"
		   );

    return (@classes);
    
}

sub wee_run {

    my @classes = (
		   "elegans_CDS",
		   "elegans_pseudogenes",
		   "elegans_RNA_genes",
		   "Transposon",
		   );

    return (@classes);
    
}



=pod

=head2   NAME - dbcomp.pl

=head1 USAGE

=over 4

=item dbcomp.pl [-options]

=back

This script will (by default) perform a class by class comparison of autoace and 
the previous WormBase release database.  It will calculate what database this
will be but you can force the choice of the second database by using the -database
option.  Classes to be compared are specified within the script and should be amended
as necessary when classes are added/removed to the database.


=over 4

=item MANDATORY arguments: none

=back

=over 4

=item OPTIONAL arguments: -debug, -help, -database


-debug and -help are standard Wormbase script options.

-database allows you to specify the path of a database which will
then be compared to the BUILD/autoace directory.

-database2 allows you to specify a second database to compare to that specified
by -database (rather than compare with previous build)


=back

=head1 AUTHOR - Dan Lawson (but completely rewritten by Keith Bradnam)


Email krb@sanger.ac.uk



=cut

