#!/software/bin/perl -w
#
# datomic_check_file.pl
#
# Outputs file of IDs in all classes for use in checking the Datomic database
#
# Last updated by: $Author: pad $
# Last updated on: $Date: 2015-07-02 13:46:46 $


use strict;
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use IO::Handle;


$|=1;

my @classes = (
	       "Map",
	       "Gene_class",
	       "Locus",
	       "Gene",
	       "Gene_cluster",
	       "Interaction",
	       "Balancer",
	       "Rearrangement",
	       "Strain",
	       "Clone",
	       "Grid",
	       "Grid_data",
	       "Grid_row",
	       "Mixed_grid_row",
	       "Contig",
	       "Species",
	       "2_point_data",
	       "Pos_neg_data",
	       "Multi_pt_data",
	       "GO_annotation",
	       "GO_term",
	       "GO_code",
	       "SO_term",
	       "Ace2SO",
	       "DO_term",
	       "PATO_term",
	       "AO_code",
	       "Reconstruction",
	       "Sequence",
	       "Sequence_collection",
	       "CDS",
	       "Transcript",
	       "Pseudogene",
	       "Transposon",
	       "Transposon_family",
	       "Library",
	       "Homology_group",
	       "Feature",
	       "Transcription_factor",
	       "Position_Matrix",
	       "Oligo",
	       "Genetic_code",
	       "Protein",
	       "Motif",
	       "Database",
	       "Method",
	       "Reference",
	       "Expr_pattern",
	       "Expr_profile",
	       "SK_map",
	       "Antibody",
	       "Picture",
	       "Anatomy_term",
	       "Anatomy_function",
	       "Life_stage",
	       "Tree",
	       "TreeNode",
	       "Transgene",
	       "Construct",
	       "RNAi",
	       "PCR_product",
	       "Phenotype",
	       "Operon",
	       "Movie",
	       "Microarray",
	       "Microarray_results",
	       "Microarray_experiment",
	       "Oligo_set",
	       "Expression_cluster",
	       "SAGE_tag",
	       "SAGE_experiment",
	       "Person",
	       "Author",
	       "Paper",
	       "Laboratory",
	       "Variation",
	       "Structure_data",
	       "Mass_spec_experiment",
	       "Mass_spec_peptide",
	       "Analysis",
	       "Condition",
	       "Molecule",
	       "WBProcess",
	      );


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($database);
my ($errfile,$outfile); 
my ($dbname);

GetOptions (
	    "help"          => \$help,
            "debug=s"       => \$debug,
	    "test"          => \$test,
	    "verbose"       => \$verbose,
	    "store:s"       => \$store,
	    "database=s"    => \$database,
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
my $exec        = $wormbase->tace;

$dbname       = "WS${WS_current}";
my $db        = $wormbase->autoace;

# alternative database specified?
if ($database) {
    $dbname  = "$database";
    $db      = "$database"; 
}

# create_log_files;

$log->write_to("\n");
$log->write_to("Current db       : $dbname '$db'\n");
$log->write_to("\n\n");
				
# open files for class and ID for use in QA when developing Datomic
my $allout = $wormbase->compare."/WS${WS_current}_all_classes.out";

&classes();

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

sub classes {

  open (OUT, ">>$allout") || $log->log_and_die("Can't open $allout\n");

  foreach my $query (@classes) {
    next if ($query eq "");
    print "$query\n";

    my $command = "query find '$query'\nlist -a\nquit\n";

    # open tace connection and count how many objects in that class
    open (TACE, "echo '$command' | $exec $db | ");
    while (my $line = <TACE>) {
      chomp $line;
      my @line = split /\s:\s/, $line;
      if (!defined $line[1]) {next}
      if ($line[0] eq 'KeySet') {next}
      if ($line[1] =~ /^"(\S+)"$/) {$line[1] = $1};
      $line[1] =~ s#\\##g;
      $line[1] = '"'.$line[1].'"';
      print OUT (join " : ", @line)."\n";
    }
    close (TACE);
    
  }
  close(OUT);
  
}

##################################################


###############################################

=pod

=head2   NAME - datomic_check_classes.pl

=head1 USAGE

=over 4

=item datomic_check_classes.pl [-options]

=back

Outputs all IDs of selected classes for use in checking the Datomic Build.

=over 4

=item MANDATORY arguments: none

=back

=over 4

=item OPTIONAL arguments: -debug, -help, -database


-debug and -help are standard Wormbase script options.

-database allows you to specify the path of a database which will
then be compared to the BUILD/autoace directory.

=back

=head1 AUTHOR - Gary Williams

=cut

