#!/usr/bin/env perl
#
# datomic_check_file.pl
#
# Outputs file of IDs in all classes for use in checking the Datomic database
#

use strict;
use Getopt::Long;
use Storable;

use lib $ENV{CVS_DIR};
use Wormbase;
use Log_files;




my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($database, $allout);

GetOptions (
  "help"          => \$help,
  "debug=s"       => \$debug,
  "test"          => \$test,
  "verbose"       => \$verbose,
  "store:s"       => \$store,
  "database=s"    => \$database,
  "outfile=s"     => \$allout,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}

my $log = Log_files->make_build_log($wormbase);

my $tace   = $wormbase->tace;
$database = $wormbase->autoace if not defined $database;
$allout = $wormbase->reports."/all_classes" if not defined $allout;

$log->write_to("\nGenerating object list for $database and writing to $allout\n\n");

				
my @classes;
while(<DATA>) {
  /^(\S+)/ and push @classes, $1;
}

&classes(\@classes);

# Email log file
$log->mail();
exit (0);


##########################################

sub classes {
  my ($class_list) = @_;
  

  open (my $fh, ">$allout") or $log->log_and_die("Can't open $allout\n");

  foreach my $class (@$class_list) {
    $log->write_to("Fetching names for $class...\n");

    my $command = "query find '$class'\nlist -a\nquit\n";

    # open tace connection and count how many objects in that class
    my $listen = 0;
    open (my $tace_fh, "echo '$command' | $tace $database | ") or $log->log_and_die("Could not open tace command\n");
    while (<$tace_fh>) {

      /^KeySet/ and do {
        $listen = 1;
        next;
      };

      if ($listen and /^(\S+) : (\S.+)$/) {
        my ($class, $value) = ($1, $2);

        $value =~ s/^\"(.+)\"/$1/; 
	$value =~ s/\\//g; # tace adds a backslash before unusual characters

        print $fh "\"$class\", \"$value\"\n";
      }
    }
    close ($tace_fh);
    
  }
  close($fh) or $log->log_and_die("Could not close $allout after writing\n");
  
}

##################################################

__DATA__
2_point_data
Ace2SO
Anatomy_function
Analysis
Anatomy_term
Antibody
AO_code
Author
Balancer
CDS
Clone
Condition
Construct
Contig
Database
DO_term
Expression_cluster
Expr_pattern
Expr_profile
Feature
Gene
Genetic_code
Gene_class
Gene_cluster
Grid
Grid_data
Grid_row
GO_annotation
GO_term
GO_code
Homology_group
Interaction
Laboratory
Library
Life_stage
Locus
Map
Mass_spec_experiment
Mass_spec_peptide
Method
Microarray
Microarray_results
Microarray_experiment
Mixed_grid_row
Molecule
Motif
Movie
Multi_pt_data
Oligo
Oligo_set
Operon
Paper
PATO_term
Person
PCR_product
Phenotype
Picture
Position_Matrix
Pos_neg_data
Protein
Pseudogene
Rearrangement
Reconstruction
Reference
RNAi
SAGE_experiment
SAGE_tag
Sequence
Sequence_collection
SK_map
SO_term
Species
Structure_data
Strain
Transposon
Transposon_family
Transcript
Transgene
Transcription_factor
Tree
TreeNode
Variation
WBProcess
