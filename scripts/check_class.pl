#!/usr/local/bin/perl5.8.0 -w
#
# dbcomp.pl
#
# Counts the number of objects in an ACEDB database for each Class stated in the config file
# Compares this number to those from a second database.
#
# Last updated by: $Author: gw3 $
# Last updated on: $Date: 2009-01-07 10:34:33 $


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
my ($database, $database1, $database2, $classes, $species, $stage);
my ($dbname_0, $dbname_1, $dbname_2);
my ($stlace, $camace, $genace, $csh, $caltech, $misc_static, $brigace, $incomplete, $data_sets);

GetOptions (
	    "help"          => \$help,
            "debug=s"       => \$debug,
	    "test"          => \$test,
	    "verbose"       => \$verbose,
	    "store:s"       => \$store,
	    "database=s"    => \$database,
	    "classes=s"     => \$classes,
	    "stlace"        => \$stlace,
	    "camace"        => \$camace,
	    "genace"        => \$genace,
	    "csh"           => \$csh,
	    "caltech"       => \$caltech,
	    "misc_static"   => \$misc_static,
	    "brigace"       => \$brigace,
	    "incomplete"    => \$incomplete,
	    "data_sets"     => \$data_sets,
	    "species:s"     => \$species,
	    "stage:s"       => \$stage,
	    );

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


############################################
# Get paths and versions
############################################

$species = $wormbase->species;
my $version = $wormbase->get_wormbase_version;
my $prev_version = $version-1;
my $prev_prev_version = $version-2;

my $exec = $wormbase->tace;

$database = $wormbase->autoace;
$dbname_0    = "WS${prev_prev_version}";
$dbname_1    = "WS${prev_version}";
$dbname_2    = "WS${version}";

my $file = $wormbase->build_data . "/COMPARE/class_count.dat"; # file holding pparse details from previous Builds



#########################################################################
# Main part of script
#########################################################################

my @classes = ();
@classes = split /,/, $classes if (defined $classes);

$stage="unknown" if (!defined $stage);

# get the collected sets of classes
 @classes = (@classes, &set_classes('stlace')) if ($stlace);
 @classes = (@classes, &set_classes('camace')) if ($camace);
 @classes = (@classes, &set_classes('genace')) if ($genace);
 @classes = (@classes, &set_classes('csh')) if ($csh);
 @classes = (@classes, &set_classes('caltech')) if ($caltech);
 @classes = (@classes, &set_classes('misc_static')) if ($misc_static);
 @classes = (@classes, &set_classes('brigace')) if ($brigace);
 @classes = (@classes, &set_classes('incomplete')) if ($incomplete);
 @classes = (@classes, &set_classes('data_sets')) if ($data_sets);

$log->write_to("Checking $dbname_1 vs $dbname_2 for classes:\n@classes\n\n");

my ($class_count) = &count_classes($database, @classes);

$log->write_to(sprintf("%-22s %7s %7s %7s %7s\n", "CLASS","($dbname_0)",$dbname_1,$dbname_2,"Difference"));

# don't want to report duplicate classes
my %seen;

my $count = 0;
foreach my $class (@classes) {

  ##################################################
  # Calculate difference between databases         #
  ##################################################
  
  my $count0 = &get_prev_count($species, $prev_prev_version, $class, $stage);
  my $count1 = &get_prev_count($species, $prev_version, $class, $stage);
  my $count2 = $$class_count[$count];
  &store_count($species, $version, $class, $stage, $count2);
  my $diff = $count2 - $count1;
  my $err = "";
  if ($count2 == 0) {
    $err = "***** POSSIBLE ERROR *****";
    $log->error;
  } elsif (	# we expect the 'incomplete' classes to be less than in currentdb
      ($count2 < $count1 * 0.9 || 
      $count2 > $count1 * 1.1)) {
    $err = "***** POSSIBLE ERROR *****";
    $log->error;
  }
  $count++;

  # don't want to report duplicate classes
  if ($seen{$class}) {next;}
  $seen{$class} = 1;

  $log->write_to(sprintf("%-22s %7d %7d %7d %7d %s\n", $class,$count0,$count1,$count2,$diff,$err));
}


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

sub count_classes {

  my ($database, @classes)  = @_;

  my @class_count;

  my $count=0;

  # Formulate query
  my $command;
  foreach my $class (@classes) {
    $command .= "query find $class\n";
  }
  $command .= "quit\n";


  ####################################
  # Count objects in database
  ####################################
  #print "Command = $command\n";
  # open tace connection and count how many objects in each class
  open (TACE, "echo '$command' | $exec $database | ");
  while (<TACE>) {
    (push @class_count, $1) if (/^\/\/ (\d+) Active Objects/);
    $count++ if (/^\/\/ (\d+) Active Objects/);
  }
  close (TACE);
    
    
  if ($count != @classes) {die "There are different numbers of classes and results (".scalar @classes." $count)\npossibly some results are not being returned?\n";}

  return (\@class_count);
    
}


###############################################
# get the count of this class found in the previous Build
sub get_prev_count {
  my ($species, $version, $class, $stage) = @_;

  my $last_count = 0;
  if (open (CLASS_COUNT, "< $file")) {
    while (my $line = <CLASS_COUNT>) {
      chomp $line;
      my ($cc_version, $cc_class, $cc_species, $cc_count, $cc_stage) = split /\t+/, $line;
      # we don't want to get the count from any details that may have
      # been stored by an earlier run of this script in this Build,
      # but we want the most recent version's count apart from that,
      # so get the most recent result that isn't from this Build.
      if ($cc_version == $version && $cc_class eq $class && $cc_species eq $species && $cc_stage eq $stage) {
	# store to get the last one in the previous build
	$last_count = $cc_count;
      } 
    }
    close (CLASS_COUNT);
  }
  return $last_count;
}

###############################################
# now store the details for this Build
sub store_count {

  my ($species, $version, $class, $stage, $count) = @_;


  if (open (CLASS_COUNT, ">> $file")) {
    if ($version && $class && $species && $count && $stage) {
      print CLASS_COUNT "$version\t$class\t$species\t$count\t$stage\n";
    } else {
      $log->write_to("*** POSSIBLE ERROR: Couldn't write to $file because some of the following is blank\nversion=$version, class=$class, species=$species, count=$count\n\n");
      $log->error;
    }
    close (CLASS_COUNT);
  } else {
    $log->write_to("WARNING: Couldn't write to $file\n\n");
  }
}


###############################################
# set the classes to check

sub set_classes {

  my ($mode) = @_;

  my @classes;

#
# The following are useful sets of classes that are loaded from the
# primary databases they will be missing some objects that are loaded
# later, so some classes are commented out to avoid alarming the
# Builder.
#

  if ($mode eq "stlace") {
    @classes = (
#		  "CDS", 
#		  "DNA",
		  "Transposon",
#		  "Transcript", 
		  "Genetic_code", 
		  "Pseudogene", 
		  "Variation", 
		  "Oligo", 
		  "PCR_product"
		  );
    
  } elsif ($mode eq "camace") {
    @classes = (
#		  "Sequence", 
#		  "CDS", 
#		  "DNA", 
#		  "Motif", 
		  "Transposon", 
#		  "Transcript", 
#		  "Feature_data", 
#		  "Sequence", 
		  "Genetic_code", 
		  "Pseudogene", 
		  "Feature", 
		  "LongText"
		  );


  } elsif ($mode eq "genace") {
    @classes = (
		  "2_point_data",
#		  "Clone",
		  "Contig",
		  "Gene_class",
		  "Grid",
		  "Laboratory",
		  "Locus",
		  "Gene",
		  "Map",
		  "Multi_pt_data",
		  "Oligo",
		  "PCR_product",
		  "Picture",
		  "Pos_neg_data",
		  "Rearrangement",
		  "Strain",
		  "Variation",
		  "View",
#		  "Sequence",
		  "Operon",
		  "Database"
		  );


  } elsif ($mode eq "csh") {
    @classes = ( 
#		 "Sequence", 
		 "LongText", 
		 "Movie", 
		 "Oligo", 
		 "PCR_product", 
		 "Structure_data", 
#		 "Peptide", 
#		 "Protein", 
		 );


  } elsif ($mode eq "caltech") {
    @classes = (
		"Transgene",
		"Expr_pattern",
		"Life_stage",
		"Cell",
		"Cell_group",
		"Paper",
		"Paper_name",
		"Journal",
		"Author",
		"Person",
		"LongText",
		"Oligo",
		"PCR_product",
#		"Sequence",
#	        "CDS",
		"Phenotype",
		"Expr_profile",
		"SK_map",
		"Tree",
		"TreeNode",
		"Microarray",
		"Microarray_results",
		"Microarray_experiment",
		"Condition",
		"Expression_cluster",
		"Oligo_set",
		"Anatomy_name",
		"Anatomy_term",
		"Homology_group",
		"Antibody",
		"RNAi",
#	        "Homol_data",
		"SAGE_tag",
		"SAGE_experiment",
		"Gene_regulation",
		"Gene",
		"GO_term",
#		"Transcript",
		"Pseudogene",
		"Interaction",
		"Variation",
		"Database",
		"YH",
		);


  } elsif ($mode eq "misc_static") {
    @classes = (
		  "GO_code",
		  "Method",
		  "Sequence MTCE*",
		  "Oligo",
		  );


  } elsif ($mode eq "brigace") {
    @classes = (
#		  "Sequence",
#		  "DNA",
#		  "Transcript",
#		  "Feature_data",
#		  "CDS",
#		  "Protein",
#		  "Peptide",
		  "Variation",
		  );


# these classes have are partially loaded in autoace_builder.pl -build
# but there will be more instances of them added later in the Build - they
# should not be empty, but they will not have as many members as the
# corresponding currentdb class.
  } elsif ($mode eq "incomplete") { 
    @classes = (
		"CDS", 
		"Clone",
		"DNA",
		"Sequence", 
		"Motif", 
		"Feature_data", 
		"Peptide", 
		"Protein", 
		"Homol_data",
		"Transcript",
		);

# these classes should all be fully loaded after running autoace_builder -data_sets
  } elsif ($mode eq "data_sets") { 
    @classes = (
		"CDS", 
		"Clone",
		"DNA",
		"Sequence", 
		"Motif", 
		"Feature_data", 
		"Peptide", 
		"Protein", 
		"Homol_data",
		"Transcript",
		"Mass_spec_peptide",
		"Mass_Spec_experiment",
		"Gene",
		"Feature",
		"Motif",
		"Method",
		);  
  }


  return @classes;


}

###############################################

=pod

=head2   NAME - check_class.pl

=head1 USAGE

=over 4

=item check_class.pl [-options]

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

-camace
-stlace
-genace
-csh
-caltech
-misc_static
-brigace test the set of classes that are loaded from the primary
databases. (Some classes from the primary databases are omitted
because further objects of those classes are loaded later on and a
test of the numbers immediately after loading the primary databases
objects would give a apparent error.)


-classes are a comma-delimited list of explicit classes to check

-stage is a name given to this stage of the Build so you can compare
 to the count of classes in the previous Build at this same stage.

-species specifies the species BUILD database to use.
=back

=head1 AUTHOR - Gary Williams



=cut

