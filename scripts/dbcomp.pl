#!/software/bin/perl -w
#
# dbcomp.pl
#
# Counts the number of objects in an ACEDB database for each Class stated in the config file
# Compares this number to those from a second database.
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


######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($database, $database2, $mail_paulsternberg);
our ($errfile,$outfile); 
our ($db_1, $db_2, $dbname_1, $dbname_2);

GetOptions (
	    "help"          => \$help,
            "debug=s"       => \$debug,
	    "test"          => \$test,
	    "verbose"       => \$verbose,
	    "store:s"       => \$store,
	    "database=s"    => \$database,
	    "database2=s"   => \$database2,
	    "mail_paulsternberg" => \$mail_paulsternberg,
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
# Check numbers of of classes
#########################################################################

&classes("Sanger and St. Louis", &sanger_stlouis);
&classes("Sanger", &sanger);
&classes("CalTech", &caltech);
&classes("Cold Spring Harbor", &csh);

#########################################################################
# Get numbers of genes curated
#########################################################################
my %clonelab;

get_curation_stats('elegans', $wormbase);
my %accs = $wormbase->species_accessors;
foreach my $spec (sort keys %accs) {
  get_curation_stats($spec, $accs{$spec});
}

close (OUT);
close (ERR);


# create symbolic links to current.out and current.dbcomp
$log->write_to("\nCreating 'current.out' and 'current.dbcomp'\n",$log);
$wormbase->run_command("rm -f ".$wormbase->compare."/current.out",$log);
$wormbase->run_command("rm -f ".$wormbase->compare."/current.dbcomp",$log);
$wormbase->run_command("ln -s $errfile ".$wormbase->compare."/current.out",$log);
$wormbase->run_command("ln -s $outfile ".$wormbase->compare."/current.dbcomp",$log);

$log->write_to("\nClass and Gene-curation Statistics are in: $outfile\n\n");

# mail the results file to Paul Sternberg
$wormbase->mail_maintainer("Class and Gene-curation Statistics for $dbname_2", 'All', $outfile);
$wormbase->mail_maintainer("Class and Gene-curation Statistics for $dbname_2", 'pws@caltech.edu', $outfile) if ($mail_paulsternberg);

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

  my ($name, @classes) = @_;

  print OUT  " +------------------------+\n";
  printf OUT  " | %20s   |\n", $name;
  print OUT  " |------------------------|---------+---------+---------+---------+---------+\n";

  my $counter = 1; # for indexing each class to be counted

  foreach my $query (@classes) {
    next if ($query eq "");
    
    $log->write_to(" Counting '$query'\n");
    printf OUT " | %22s |", $query;
    
    ##################################################
    # Get class counts from both databases           #
    ################################################## 
    
    my ($class_count_1,$class_count_2) = &count_class($query,$counter);
    
    $log->write_to(" Counting $dbname_1 : $class_count_1\n");  
    $log->write_to(" Counting $dbname_2 : $class_count_2\n");
    
    printf OUT " %7s ",$class_count_1;
    print OUT "|";  
    printf OUT " %7s ",$class_count_2;
    
    
    ##################################################
    # Calculate difference between databases         #
    ##################################################
    
    my ($added, $removed) = &diff($counter);
    my $diff = $added - $removed;
    
    #  printf OUT "| %6s |\n",$diff;
    $log->write_to(" Diff $diff\n\n");

    printf OUT "| %7s ", $added;
    printf OUT "| %7s ", $removed;
    printf OUT "| %7s |\n",$diff;
    
    
    $counter++;
  }
  
  
  print OUT  " +------------------------+---------+---------+---------+---------+---------+\n";

}

##################################################

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

sub sanger_stlouis {
  my @classes = (
		 "Sequence",
		 "CDS",
		 "Transposon",
		 "Transcript",
		 "Pseudogene",
		 "PCR_product",
		 "Transposon_CDS",
		 "cDNA_sequence",
		 "elegans_CDS",
		 "elegans_pseudogenes",
		 "elegans_RNA_genes",
		 );
}

sub sanger {
  my @classes = (
		 "Class",
		 "Model",
		 "Method",
		 "Clone",
		 "Coding_transcripts",
		 "Accession_number",
		 "Variation",
		 "Motif",
		 "Feature",
		 "Feature_data",
		 "Laboratory",
		 "Locus",
		 "Gene",
		 "Gene_class",
		 "Gene_name",
		 "Gene_regulation",
		 "Genome_Sequence",
		 "Map",
		 "Multi_pt_data",
		 "Picture",
		 "Pos_neg_data",
		 "Rearrangement",
		 "2_point_data",
		 "Strain",
		 "Operon",
		 "Analysis",
		 "Condition",
		 "GO_code",
		 "Peptide",
		 "Protein",
		 "Species",
		 "Transposon_family",
		 "Comment",
		 "Contig",
		 "Database",
		 "Display",
		 "DNA",
		 "Homol_data",
		 "NDB_Sequence",
		 "nematode_ESTs",
		 "SO_term",
		 "Table",
		 "Mass_spec_data",
		 "Mass_spec_experiment",
		 "Mass_spec_peptide",
		 );
}

sub caltech {
  my @classes = (
		 "Transgene",
		 "Expr_pattern",
		 "Expr_profile",
		 "Life_stage",
		 "Lineage",
		 "Cell",
		 "Cell_group",
		 "Paper",
		 "Paper_name",
		 "Author",
		 "Person",
		 "Person_name",
		 "LongText",
		 "Keyword",
		 "Oligo",
		 "Oligo_set",
		 "Phenotype",
		 "Phenotype_name",
		 "SK_map",
		 "Tree",
		 "TreeNode",
		 "Microarray",
		 "Microarray_results",
		 "Microarray_experiment",
		 "Expression_cluster",
		 "Anatomy_name",
		 "Anatomy_term",
		 "Anatomy_function",
		 "Homology_group",
		 "Antibody",
		 "RNAi",
		 "SAGE_tag",
		 "SAGE_experiment",
		 "Gene_regulation",
		 "GO_term",
		 "Interaction",
		 "YH",
		 "Position_matrix",
		 );
}

sub csh {
  my @classes = (
		 "LongText",
		 "Movie",
		 "Structure_data",
		 );
}



###############################################

sub get_curation_stats {

  my ($species, $wb) = @_;

  # get the centre that sequenced the clones
  %clonelab = $wormbase->FetchData("clone2centre", undef, "$db_2/COMMON_DATA");

  my $prefix = $wb->pepdir_prefix;
  my $release = $WS_current;

  my $file = $wb->wormpep . "/${prefix}pep.diff${release}";

  if (!open (IN, "<$file")) {
    print OUT "\n\nNo curation data available for $species.\n\n";
  } else {
    
    my %lab_count;
    $lab_count{HX} = $lab_count{RW} = 0;
    $lab_count{HX_lost} = $lab_count{RW_lost} = 0;
    $lab_count{HX_changed} = $lab_count{RW_changed} = 0;
    $lab_count{HX_new_gene} = $lab_count{RW_new_gene} = 0;
    $lab_count{HX_new_isoform} = $lab_count{RW_new_isoform} = 0;
    
    while (my $line = <IN>) {
      my @f = split /\s+/, $line;
      my $type = $f[0];
      $type =~ s/://;
      my ($clone) = ($f[1] =~ /^(\S+)\.\S+/);
      #print $clone, "\n";
      my $lab;
      if ($species eq 'elegans') {
	$lab = &get_lab($clone);
      } else {
	$lab = 'RW';
      }
      #print "$clone $lab\n";
      $lab_count{$lab}++;
      
      if ($type eq 'lost') {
	$lab_count{"${lab}_lost"}++;
	print "lost $f[1] $lab\n";

      } elsif ($type eq 'changed') {
	$lab_count{"${lab}_changed"}++;
	print "changed $f[1] $lab\n";
	
      } elsif ($type eq 'new') {
	if ($f[1] =~ /[a-z]$/) { # new isoform
	  $lab_count{"${lab}_new_isoform"}++;
	  print "new isoform $f[1] $lab\n";
	  
	  # if we have a new first isoform, we should ignore the old gene that is now lost
	  if ($f[1] =~ /a$/) {
	    $lab_count{"${lab}_lost"}--;
	    $lab_count{$lab}--;
	  }
	  
	} else { # new gene (usually a result of splitting old gene)
	  $lab_count{"${lab}_new_gene"}++;
	  print "new gene $f[1] $lab\n";

	}
      } elsif ($type eq 'reappeared') { # ignore for now
      } else {
	die "Unrecognised type: $type\n";
      }
      
    }
  

    print OUT "\n\n-----------------------------------------\n";
    print OUT "Gene model curation changes for $species\n";
    print OUT "-----------------------------------------\n\n";

    my $RW = $lab_count{RW_lost} + $lab_count{RW_changed} + $lab_count{RW_new_gene} + $lab_count{RW_new_isoform};
    my $HX = $lab_count{HX_lost} + $lab_count{HX_changed} + $lab_count{HX_new_gene} + $lab_count{HX_new_isoform};
    my $total = $HX + $RW;
    my $total_lost = $lab_count{HX_lost} + $lab_count{RW_lost};
    my $total_changed = $lab_count{RW_changed} + $lab_count{HX_changed};
    my $total_new_gene = $lab_count{RW_new_gene} + $lab_count{HX_new_gene};
    my $total_iso = $lab_count{RW_new_isoform} + $lab_count{HX_new_isoform};

    if ($species eq 'elegans') {
      print OUT "Total $total\nlost: $total_lost ($lab_count{HX_lost}/$lab_count{RW_lost})\nchanged: $total_changed ($lab_count{HX_changed}/$lab_count{RW_changed})\nnew gene: $total_new_gene ($lab_count{HX_new_gene}/$lab_count{RW_new_gene})\nnew isoform: $total_iso ($lab_count{HX_new_isoform}/$lab_count{RW_new_isoform})\n";
#      print OUT "Total ", $total,"\nHinxton: $HX\nSt. Louis: $RW\n\nlost (Hinxton): $lab_count{HX_lost}\nlost (St. Louis): $lab_count{RW_lost}\n\nchanged (Hinxton): $lab_count{HX_changed}\nchanged (St. Louis): $lab_count{RW_changed}\n\nnew gene (Hinxton): $lab_count{HX_new_gene}\nnew gene (St. Louis): $lab_count{RW_new_gene}\n\nnew isoform (Hinxton): $lab_count{HX_new_isoform}\nnew isoform (St. Louis): $lab_count{RW_new_isoform}\n";
    } else {
      print OUT "Total $RW\nlost: $lab_count{RW_lost}\nchanged: $lab_count{RW_changed}\nnew gene: $lab_count{RW_new_gene}\nnew isoform: $lab_count{RW_new_isoform}\n";
    }
  }

}

###############################################

sub get_lab {

  my ($clone_name) = @_;

  my $centre = $clonelab{$clone_name};          # get the lab that sequenced this clone
  if (!defined $centre) {return 'RW'}
  return $centre;
}

###############################################

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

=head1 AUTHOR - Dan Lawson (but completely rewritten by Keith Bradnam, and rewritten again by Gary Williams)


Email krb@sanger.ac.uk



=cut

