#!/usr/bin/env perl
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
my ($database, $database2, $mail_paulsternberg, $post_gff, $pre_rel, $hinxton, $caltech, $hxcurated);
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
	    "post_gff"      => \$post_gff,
	    "pre_rel"       => \$pre_rel,
            "hinxton"       => \$hinxton,
            "hxcurated"     => \$hxcurated,
            "caltech"       => \$caltech,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
			     );
}


&usage("Help") if ($help);

my $log = Log_files->make_build_log($wormbase);


############################################
# Get paths and versions
############################################

my $WS_current  = $wormbase->get_wormbase_version;
my $WS_previous = $WS_current - 1;
my $exe        = $wormbase->tace;
my $basedir         = $wormbase->basedir;

#################################################################
# Compare autoace to previous build unless -database specified
#################################################################
unless ((defined $pre_rel) || (defined$database) || (defined$post_gff)) {
  $log->log_and_die("One of the work options needs to be specified -database/-post_gff/-pre_rel\n")
} 

if ($pre_rel) {
  $dbname_1    = "WS${WS_previous}";
  $db_1        = $wormbase->database("WS${WS_previous}");
}

#Now have an elegans backup from the pre-merge
if ($post_gff){
  my $sp = "elegans";
  $dbname_1    = "BUILD/elegans(WS${WS_previous})";
  $db_1        = "$basedir/elegans";
}

$dbname_2    = "BUILD/autoace(WS${WS_current})";
$db_2        = $wormbase->autoace;

# First alternative database specified?
if ($database) {
    $dbname_1  = "$database";
    $db_1      = "$database";
    unless (-e "$database" ) {$log->log_and_die("\nERROR: database appears to be invalid\n")} 
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
print OUT  " |                        |   ${WS_previous}   |   ${WS_current}   |    +    |    -    |   Net   |\n";
print OUT  " +------------------------+---------+---------+---------+---------+---------+\n";

#########################################################################
# Check numbers of of classes
#########################################################################

my @all_classes;
if ($hinxton) {
    push @all_classes, &hinxton();
}
elsif ($caltech) {
    push @all_classes, &caltech();
}
elsif ($hxcurated) {
    push @all_classes, &hxcurated();
}
else {
    push @all_classes, &hinxton();
    push @all_classes, &caltech();
}
&class_summary( sort @all_classes );

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

sub class_summary {

  my (@classes) = @_;

  my $counter = 1; # for indexing each class to be counted

  foreach my $query (@classes) {
    next if ($query eq "");
    
    $log->write_to(" Counting '$query'\n");
    printf OUT " | %22s |", $query;
    
    ##################################################
    # Get class counts from both databases           #
    ################################################## 
    
    my ($class_count_1,$class_count_2) = &count_class($query,$counter);
    my $err = "CLEAR";
    if ($class_count_2 == 0) {
      $err = "***** POSSIBLE ERROR *****";
      $log->error;
    } 
    elsif (
	   ($class_count_2 < $class_count_1 * 0.9 || 
	    $class_count_2 > $class_count_1 * 1.1)) {
      $err = "***** POSSIBLE ERROR *****";
      $log->error;
    }
    
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
    $log->write_to(" Diff $diff\n $err\n\n");
    
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
    open (TACE, "echo '$command' | $exe $db_1 | ");
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
    open (TACE, "echo '$command' | $exe $db_2 | ");
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

sub hxcurated {
    my @classes = (
	"elegans_CDS",
	"elegans_pseudogenes",
	"elegans_RNA_genes",
	"Transposon",
	"Transposon_CDS",
	"Transposon_family",
	"Variation",
	"Feature",
	"Gene",
	"Gene_class",
	"Strain",
	"Operon",
	);
return @classes;	
}

sub hinxton {
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
    "Class",
    "Model",
    "Method",
    "Clone",
    "Coding_transcripts",
    "Variation",
    "Motif",
    "Feature",
    "Feature_data",
    "Laboratory",
    "Locus",
    "Gene",
    "Gene_class",
    "Gene_name",
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
    "Mass_spec_experiment",
    "Mass_spec_peptide",
    "Structure_data",
    "LongText"
      );

  return @classes;
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
    "Author",
    "Person",
    "Person_name",
    "LongText",
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
    "GO_term",
    "Interaction",
    "Position_matrix",
      );
  return @classes;
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

    my $total = 0;
    my $lost = 0;
    my $changed = 0;
    my $new_gene = 0;
    my $new_isoform = 0;
    
    while (my $line = <IN>) {
      my @f = split /\s+/, $line;
      my $type = $f[0];
      $type =~ s/://;

      if ($type eq 'lost') {
        $lost++;
	print "lost $f[1]\n";

      } elsif ($type eq 'changed') {
        $changed++;
	print "changed $f[1]\n";
	
      } elsif ($type eq 'new') {
	if ($f[1] =~ /[a-z]$/) { # new isoform
          $new_isoform++;
	  print "new isoform $f[1]\n";
	  
	  # if we have a new first isoform, we should ignore the old gene that is now lost
	  if ($f[1] =~ /a$/) {
            $lost--;
            $total--;
	  }
	  
	} else { # new gene (usually a result of splitting old gene)
          $new_gene++;
	  print "new gene $f[1]\n";

	}
      } elsif ($type eq 'reappeared') { 
        # ignore for now
      } else {
	$log->log_and_die("Unrecognised type: $type\n");
      }
    }
  

    print OUT "\n\n-----------------------------------------\n";
    print OUT "Gene model curation changes for $species\n";
    print OUT "-----------------------------------------\n\n";
    print OUT "Total $total\n";
    print OUT "lost: $lost\n";
    print OUT "changed: $changed\n";
    print OUT "new gene: $new_gene\n";
    print OUT "new isoform: $new_isoform\n";
  }
}
