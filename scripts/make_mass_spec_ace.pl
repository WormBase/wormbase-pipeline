#!/software/bin/perl -w
#
# make_mass_spec_ace.pl
# 
# by Gary Williams                         
#
# This parses a file of mass spec peptide data and creates an ace file
# in ace
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2013-12-09 11:31:14 $      

use strict;                                      
use lib $ENV{'CVS_DIR'};
use Wormbase;
use Getopt::Long;
use Carp;
use Log_files;
use Storable;
use Ace;
#use Sequence_extract;
use Coords_converter;



######################################
# variables and command-line options # 
######################################

my ($help, $debug, $test, $verbose, $store, $wormbase, $species);
my ($input, $output, $database, $directory, $load);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "input:s"    => \$input,
	    "output:s"   => \$output,
	    "database:s" => \$database,
	    "directory:s"  => \$directory,
	    "load"       => \$load,
            "species:s"  => \$species, # the default is elegans
	   );

$species = 'elegans' unless $species;

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

# check parameters
if (! defined $database) {$database = $wormbase->autoace}
if (! defined $input && ! defined $directory) {$directory = $wormbase->misc_static . "/MASS-SPEC/$species/";} # /nfs/wormpub/BUILD_DATA/MISC_STATIC/ or /nfs/wormpub/TEST/BUILD_DATA/MISC_STATIC/
if (! defined $output) {$output = $wormbase->acefiles."/mass-spec-data.ace";}

my $wormpep_prefix = $wormbase->wormpep_prefix;  # colon prefix to indicate the database e.g. 'WP:' in elegans

my $history_suffix = lc($wormpep_prefix);
# japonica uses 'jp' not 'ja' to name the history objects - blame Phil for this!
if ($wormpep_prefix eq 'ja') {$wormpep_prefix = 'jp'}
if ($wormpep_prefix eq 'cn') {$wormpep_prefix = 'np'}

my $pep_prefix = $wormbase->pep_prefix; # start of protein ID e.g. 'CE' in elegans

my $pepdir_prefix = $wormbase->pepdir_prefix; # e.g. 'worm' or 'brug'

my $cds_regex = $wormbase->cds_regex;

##########################
# MAIN BODY OF SCRIPT
##########################

#
# 7 June 2006 - need to add in the genome mapping data
# This consists of:
# 
# a single Homol_data line to be added to the clone, eg:
# 
# Sequence : clone_name
# Homol_data homol_data_object_name 1 clone_length
# 
# e.g:
# 
# Sequence : AC3
# Homol_data AC3:Mass-spec 1 19177
# 
# In the Homol_data object, add a MSPeptide_homol line for each peptide
# aligning to that clone:
# 
# Homol_data : homol_data_object_name
# MSPeptide_homol peptide_id mass_spec_genome peptide_score clone_start clone_end 1 peptide_length AlignPepDNA
# (repeated for each alignment)
# 
# e.g.:
# 
# Homol_data : AC3:Mass-spec
# MSPeptide_homol MILGCDLKYV mass_spec_genome 100.0 30100 30144 1 15 AlignPepDNA
# (repeated for each alignment)
# 

# Also need to add the method (mass_spec_genome in the stuff above) to
# the methods file in MISC_STATIC. Something like:
# 
# Method : mass_spec_genome
# Colour BROWN
# Right_priority 1.5
# (check with other methods for the correct values)
# 
# 

my $database_version = $wormbase->get_wormbase_version;
if ($database_version == 666) {
  print "In test - setting database version to 257 not 666\n";
  $database_version = 257;
}

# open an ACE connection to parse details for mapping to genome
my $tace            = $wormbase->tace;        # TACE PATH
print "Connecting to Ace\n";
my $db = Ace->connect (-path => $database,
                       -program => $tace) || die "cannot connect to database at $database\n";


# open the output ACE file
open (OUT, "> $output") || die "Can't open output file $output\n";
print OUT "\n";			# blank line at start just in case we concatenate this to another ace file

# get the wormpep sequences
print "Reading wormpep sequences\n";
my $wormpep = &get_wormpep_by_wp_id;

print "Reading wormpep history file\n"; 
# get the wormpep protein IDs used by previous version of the gene
my $wormpep_history = &get_wormpep_history;

# get the hash of wormpep sequence names keyed by old wormpep protein IDs
my $wormpep_history_proteins = &get_wormpep_history_proteins($wormpep_history);

###############################
# parse the data files
###############################

# this gives a hash with key=experiment ID value=list of details
my ($experiment_hashref, $peptide_hashref) = &parse_data();

my %experiment = %{$experiment_hashref};
my %peptide = %{$peptide_hashref};
my $peptide_id_count=1;
my %peptides_used;		# quick record of the peptides that have been used, for debugging to produce small test outputs
my %unique_peptides;		# non-redundant set of peptides used in all the experiments
my %unique_clones;  		# non-redundant set of clones mapped to in all the experiments

# get the peptides used in each protein
my %proteins;
my $experiment_id;

my %peptide_count = ();		# hash to hold count of times a peptide maps to proteins

# get the set of unique peptides from all experiments
# %unique_peptides is lists of CDS matching the key peptide
foreach my $experiment_id (keys %experiment) {
  my @peptides = @{ $experiment{$experiment_id}{'PEPTIDES'} };
  # get the unique non-redundant set of peptides in all experiments
  foreach my $pep (@peptides) {
    $unique_peptides{$pep} = 1;
  }
}

# get the proteins and their peptides
# %proteins is lists of peptides matching the key CDS
foreach my $peptide_sequence (keys %unique_peptides) {
  if (! defined $peptide{$peptide_sequence} ) { # sometimes we did't identify the CDSs for a peptide sequence
    #print "No CDS names known for peptide $peptide_sequence in experiment $experiment_id\n";
    next
  } 
  my @CDS_names = @{ $peptide{$peptide_sequence} };
  foreach my $CDS (@CDS_names) {
    push @{ $proteins{$CDS} }, $peptide_sequence;
  }
}

my $final_count_ok = 0;
my $final_count_not_ok = 0;

print "Mapping peptides to genes\n";
# go through the proteins mapping peptides to them
foreach my $CDS_name (keys %proteins) {

  my $not_found_peptides;         # href of peptides that we couldn't map to the gene			
  #print "Processing CDS $CDS_name\n";

  my $cds_processed_ok;

  # do the normal CDSs first
  if ($CDS_name !~ /^clone:/) {  

    # strip off any version number and we will search all versions
    $CDS_name =~ s/:${history_suffix}\d+//;  # e.g. strip off ':wp256'


    my $CDS_name_isoform = $CDS_name; # set the CDS_name to the original name
    ($cds_processed_ok, $not_found_peptides) = &process_cds($CDS_name, $CDS_name, $wormpep_history, \%unique_clones, %proteins);

    # see if there is an isoform that we could investigate also
    if (!$cds_processed_ok) {
      $CDS_name_isoform = "${CDS_name}a";
      if (&has_current_history($wormpep_history, $CDS_name_isoform)) {
	($cds_processed_ok, $not_found_peptides) = &process_cds($CDS_name, $CDS_name_isoform, $wormpep_history, \%unique_clones, %proteins);
      }
    }
    
  }

  # do the clones and the peptides that map to their ORFs separately  
  # this way we can add peptides that don't map to a CDS as (possibly)
  # mapping to an ORF in the clone

  if ($CDS_name =~ /^clone:/) {  
    ($cds_processed_ok, $not_found_peptides) = &process_clones($CDS_name, \%unique_clones, %proteins);
  }


  if ($cds_processed_ok) {
    $final_count_ok++;
#    $log->write_to("We have found all the peptides in CDS: $CDS_name_isoform\n");
  } else {
    $final_count_not_ok++;
    my @not_found_peptides = (keys %{$not_found_peptides});
    print "********* NOT MAPPED $CDS_name with @not_found_peptides\n";
#    $log->write_to("We have not found all the peptides in CDS: $CDS_name_isoform\n");
  }
}

# go through the peptides writing out their sequences
# whether the peptide is mapped uniquely for the Mass_spec_peptide in this experiment
# and whether or not they are flagged as natural
foreach my $ms_peptide (keys %unique_peptides) {
  if (exists $peptides_used{$ms_peptide}) { # see if used this peptide - for making smaller output when debugging
    # output the MS_peptide object
    print OUT "\n";	
    print OUT "Mass_spec_peptide : \"MSP:$ms_peptide\"\n";
    print OUT "Peptide \"MSP:$ms_peptide\"\n";
    print OUT "Protein_seq \"MSP:$ms_peptide\"\n";
    foreach my $experiment_id (keys %experiment) {	  
      print OUT "Petide_is_natural\n" if (exists $experiment{$experiment_id}{natural}); # note spelling mistake
      if (exists $peptide_count{$ms_peptide}{$experiment_id} && $peptide_count{$ms_peptide}{$experiment_id} == 1) {
	print OUT "Mass_spec_experiments \"$experiment_id\" Matches_database_uniquely\n";
      }
    }
    # and output the Peptide sequence object 
    print OUT "\n";	
    print OUT "Peptide : \"MSP:$ms_peptide\"\n";
    print OUT &no_underscore($ms_peptide) . "\n"; 
  }
}

# write out the Homol_data lines for the clones and superlinks used
my $coords      = Coords_converter->invoke($database, 0, $wormbase);

foreach my $clone (keys %unique_clones) {
  my $size = $coords->Superlink_length($clone);
  print OUT "\n";	
  print OUT "Sequence : $clone\n";
  print OUT "Homol_data $clone:Mass-spec 1 $size\n";
}

# output the experiment data objects
foreach my $experiment_id (keys %experiment) {
  print OUT "\n";
  print OUT "Mass_spec_experiment : \"$experiment_id\"\n";
  print OUT "Species \"$experiment{$experiment_id}{species}\"\n" if (exists $experiment{$experiment_id}{species});
  print OUT "Strain \"$experiment{$experiment_id}{strain}\" \n" if (exists $experiment{$experiment_id}{strain});
  print OUT "Life_stage \"$experiment{$experiment_id}{life_stage}\"\n" if (exists $experiment{$experiment_id}{life_stage});
  print OUT "Sub_cellular_localization \"$experiment{$experiment_id}{sub_cellular_localization}\"\n" if (exists $experiment{$experiment_id}{sub_cellular_localization});
  print OUT "Digestion \"$experiment{$experiment_id}{digestion}\"\n" if (exists $experiment{$experiment_id}{digestion});
  print OUT "Instrumentation \"$experiment{$experiment_id}{instrumentation}\"\n" if (exists $experiment{$experiment_id}{instrumentation});
  print OUT "Ionisation_source \"$experiment{$experiment_id}{ionisation_source}\"\n" if (exists $experiment{$experiment_id}{ionisation_source});
  print OUT "Multiple_ambiguous_IDs_allowed\n" if (exists $experiment{$experiment_id}{multiple_ambiguous_IDs_allowed});
  print OUT "Remark \"$experiment{$experiment_id}{remark}\"\n" if (exists $experiment{$experiment_id}{remark});
  print OUT "Person \"$experiment{$experiment_id}{person}\" \n" if (exists $experiment{$experiment_id}{person});
  print OUT "Author \"$experiment{$experiment_id}{author}\" \n" if (exists $experiment{$experiment_id}{author});
  print OUT "Reference \"$experiment{$experiment_id}{reference}\" \n" if (exists $experiment{$experiment_id}{reference});
  print OUT "Database \"$experiment{$experiment_id}{database}\"\n" if (exists $experiment{$experiment_id}{database});
  print OUT "Program \"$experiment{$experiment_id}{program}\"\n" if (exists $experiment{$experiment_id}{program});
  print OUT "Laboratory \"$experiment{$experiment_id}{laboratory}\" \n" if (exists $experiment{$experiment_id}{laboratory});
  print OUT "Genotype \"$experiment{$experiment_id}{genotype}\" \n" if (exists $experiment{$experiment_id}{genotype});
  print OUT "Anatomy_term \"$experiment{$experiment_id}{anatomy_term}\"\n" if (exists $experiment{$experiment_id}{anatomy_term});
  print OUT "Cell_type \"$experiment{$experiment_id}{cell_type}\" \n" if (exists $experiment{$experiment_id}{cell_type});
  print OUT "Ionisation_source \"$experiment{$experiment_id}{ionisation_source}\" \n" if (exists $experiment{$experiment_id}{ionisation_source});
  print OUT "Minimum_ion_proportion $experiment{$experiment_id}{minimum_ion_proportion}\n" if (exists $experiment{$experiment_id}{minimum_ion_proportion});
  print OUT "False_discovery_rate $experiment{$experiment_id}{false_discovery_rate}\n" if (exists $experiment{$experiment_id}{false_discovery_rate});
}

close (OUT);

# close the ACE connection
$db->close;

# load the ace file
if ($load) {
  $wormbase->load_to_database($database, $output, 'make_mass_spec_ace.pl', $log);
}


# Close log files and exit
$log->write_to("\n\nStatistics\n");
$log->write_to("----------\n\n");
$log->write_to("Count of genes with peptides mapped OK:     $final_count_ok\n");
$log->write_to("Count of genes with peptides NOT mapped OK: $final_count_not_ok (due to genes becoming pseudogenes etc. - can ignore this)\n");

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################

# example peptide records
# Mass_spec_experiment = BB_1
# Species = Caenorhabditis elegans
# Strain = N2
# Life_stage = all stages
# Sub_cellular_localization = whole organism
# Digestion = Trypsin
# Instrumentation = QTOF
# Multiple_ambiguous_IDs_allowed
# Person = WBPerson8676
# #Reference = WBPaper0001234
# Database = wormpep
# Program = Sequest
# #Natural
# Posttranslation_modification_type = Phosphorylation site
# Posttranslation_modification_database = PhosphoPep
# Posttranslation_modification_URL = http://www.sbeams.org/dev1/sbeams/cgi/Glycopeptide/peptideSearch.cgi?action=Show_hits_form&search_type=gene_symbol&search_term=<PROTEINID>
# 
# # Phosphorylation is indicated by a '_' after the modified residue
# 
# PEPTIDES
# 
# 2RSSE.2	AAAFAS_NPSHSLDYQEVGASNPR	.	. BB_1
# 2RSSE.2	AAAFASNPS_HSLDYQEVGASNPR	.	. BB_1
# AC3.5	GPAKDDEMESS_EEQE	.	. BB_1
# AC3.5	IWTRPEVK	.	. BB_1
# AC3.5	MKGPAKDDEMES_SEEQE	.	. BB_1
# AC3.5	MKGPAKDDEMES_S_EEQE	.	. BB_1
# AC3.5	MKGPAKDDEMESS_EEQE	.	. BB_1
# AC3.5	TQLIDFVYANGNGS_ASNLNNR	.	. BB_1
# AC3.5	TQLIDFVYANGNGSASNLNNR	.	. BB_1
# AC7.2a	SKS_PGGIVGR	.	. BB_1
# AH10.3	RAS_AAPNVENR	.	. BB_1
# AH9.3	EVQEAVKTEEAGSS_PK	.	. BB_1


##########################################
# this gives a hash with key=experiment ID value=list of details
# %experiment = parse_data();

sub parse_data {

  my %experiment;
  my %peptide;

  my %wormpep2cds = $wormbase->FetchData("wormpep2cds", undef, $database . "/COMMON_DATA");
  my %cgc_names = $wormbase->FetchData("cds2cgc", undef, $database . "/COMMON_DATA");
  my %cds_cgc_names;
  foreach my $cgc (keys %cgc_names) { # invert the cgc_names hash
    my $cds = $cgc_names{$cgc};
    push @{$cds_cgc_names{$cds}}, $cgc;
  }

  my %uniprot2wormpep = get_uniprot();

  # if we have a directory specified, look there for files ending *.ms-dat
  if (defined $directory) {
    opendir (DIR, $directory) or $log->log_and_die("cant open directory $directory\n");
    $log->write_to("Reading input files from $directory\n");
    my @files = readdir DIR;
    foreach my $file ( @files ) {
      next if( $file eq '.' or $file eq '..');
      if( (-T $directory."/".$file) and substr($file,-7,7 ) eq ".ms-dat" ) {
	&parse_file($directory."/".$file, \%wormpep2cds, \%cgc_names, \%cds_cgc_names, \%uniprot2wormpep, \%experiment, \%peptide);
      }
    }

    close DIR;

  } else {			# else read one input file
    $log->write_to("Reading input file\n");
    &parse_file($input, \%wormpep2cds, \%cgc_names, \%cds_cgc_names, \%uniprot2wormpep, \%experiment, \%peptide);

  }

  return (\%experiment, \%peptide);

}

##########################################
# get the data out of a file
sub parse_file {
  my ($input_file, $wormpep2cds, $cgc_names, $cds_cgc_names, $uniprot2wormpep, $experiment, $peptide) = @_;

  my $peptide_sequence;
  my $peptide_probability;
  my $protein_probability;
  my $strain_name;
  my $stage_name;
  my $sub_cellular_fraction;
  my $purification_type;
  my $quantification_type;
  my $CDS_name;
  my @CDS_names;

  my $non_existent = 0;
  my $not_in_wormbase = 0;

  my $experiment_id = undef;

  my $got_peptides = 0;		# flag for got to the PEPTIDES section of the input data
  
  $log->write_to("Reading input file $input_file\n");

  open (DATA, "< $input_file") || die "Can't open $input_file\n";
  while (my $line = <DATA>) {
    chomp $line;
    if ($line =~ /^\#/) {next;} 
    if ($line =~ /^\s*$/) {next;}

    my $cgc_name;

    if (!$got_peptides) {
      if ($line =~ /^PEPTIDES\s*$/) { # got to the PEPTIDES section of the input data
	$got_peptides = 1;

      } elsif ($line =~ /Mass_spec_experiment\s*=\s*(.+)/) {
	# usually formed from the initials of the author, plus a number, e.g. 'SH_1'
	$experiment_id = $1;
	$experiment->{$experiment_id}{experiment_id} = $1;

      } elsif ($line =~ /Species\s*=\s*(.+)/) {
	$experiment->{$experiment_id}{species} = $1;

      } elsif ($line =~ /Strain\s*=\s*(.+)/) {
	$experiment->{$experiment_id}{strain} = $1;

      } elsif ($line =~ /Life_stage\s*=\s*(.+)/) {
	$experiment->{$experiment_id}{life_stage} = $1;

      } elsif ($line =~ /Sub_cellular_localization\s*=\s*(.+)/) {
 	$experiment->{$experiment_id}{sub_cellular_localization} = $1;

      } elsif ($line =~ /Digestion\s*=\s*(.+)/) {
 	$experiment->{$experiment_id}{digestion} = $1;

      } elsif ($line =~ /Instrumentation\s*=\s*(.+)/) {
 	$experiment->{$experiment_id}{instrumentation} = $1;

      } elsif ($line =~ /Ionisation_source\s*=\s*(.+)/) {
 	$experiment->{$experiment_id}{ionisation_source} = $1;

      } elsif ($line =~ /False_discovery_rate\s*=\s*(.+)/) {
 	$experiment->{$experiment_id}{false_discovery_rate} = $1;

      } elsif ($line =~ /Multiple_ambiguous_IDs_allowed/) {
 	$experiment->{$experiment_id}{multiple_ambiguous_IDs_allowed} = 1;

      } elsif ($line =~ /Person\s*=\s*(.+)/) {
 	$experiment->{$experiment_id}{person} = $1;

      } elsif ($line =~ /Reference\s*=\s*(.+)/) {
 	$experiment->{$experiment_id}{reference} = $1;

      } elsif ($line =~ /Database\s*=\s*(.+)/) {
 	$experiment->{$experiment_id}{database} = $1;

      } elsif ($line =~ /Remark\s*=\s*(.+)/) {
 	$experiment->{$experiment_id}{remark} = $1;

      } elsif ($line =~ /Program\s*=\s*(.+)/) {
 	$experiment->{$experiment_id}{program} = $1;

      } elsif ($line =~ /Natural/) {
	# specifies that the peptides in this data set are naturally occuring.
 	$experiment->{$experiment_id}{natural} = 1;

      } elsif ($line =~ /Posttranslation_modification_type\s*=\s*(.+)/) {
 	$experiment->{$experiment_id}{posttranslation_modification_type} = $1;

      } elsif ($line =~ /Posttranslation_modification_database\s*=\s*(.+)/) {
 	$experiment->{$experiment_id}{posttranslation_modification_database} = $1;

      } elsif ($line =~ /Posttranslation_modification_URL\s*=\s*(.+)/) {
 	$experiment->{$experiment_id}{posttranslation_modification_URL} = $1;

      } else {
	die "Input line not recognised:\n$line\n";
      }

    } else {

      if ($line =~ /\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
	my @f = split /\s+/, $line;
	$cgc_name = shift @f;
	$peptide_sequence = shift @f;
	while ($protein_probability = shift @f) {
	  if ($protein_probability eq '.') {$protein_probability = undef}
	  $peptide_probability = shift @f;
	  if ($peptide_probability eq '.') {$peptide_probability = undef}
	  $experiment_id = shift @f;
	  if (!defined $experiment_id || $experiment_id eq '.') {
	    print "line $line has a non-existent experiment ID\n";
	  }

	  # save the peptide data in all its experiments
	  push @{ $experiment->{$experiment_id}{'PEPTIDES'} }, $peptide_sequence;	# list of peptides in this experiment
	  $experiment->{$experiment_id}{HAS_PEPTIDE}{$peptide_sequence} = 1; # unique hash of peptides in this experiment
	  $experiment->{$experiment_id}{'PROTEIN_PROBABILITY'}{$peptide_sequence} = $protein_probability;
	  $experiment->{$experiment_id}{'PEPTIDE_PROBABILITY'}{$peptide_sequence} = $peptide_probability;
	}
	@CDS_names = ();
      }

      if ($cgc_name =~ /\S+\-\d+/ && $cgc_name !~ /^UNIPROT/) { # standard CGC-style name e.g. daf-1
	if (!exists $cds_cgc_names->{$cgc_name}) {
	  print "line $line has a non-existent cds_cgc name\n";
	  next;
	}
	# start getting the protein data
	push @CDS_names, @{$cds_cgc_names->{$cgc_name}}; # get the list of CDS names for this CGC name
  
      } elsif ($cgc_name =~ /^${pep_prefix}\d{5}$/) { # wormpep ID starting with 'CE'  (in elegans)
	if (!exists $wormpep2cds->{$cgc_name}) {
	  if (! exists $wormpep_history_proteins->{$cgc_name}) {
	    print "line $line has a non-existent wormpep ID\n";
	    next;
	  } else {
	    push @CDS_names, $wormpep_history_proteins->{$cgc_name}; # the CDS names for this wormpep ID	
	  }
	} else {
	  push @CDS_names, (split /\s+/, $wormpep2cds->{$cgc_name}); # the CDS names for this wormpep ID
	}

      } elsif ($cgc_name =~ /$cds_regex/) {	# standard Sequence style CDS name e.g. C25A7.2 (in elegans)
	push @CDS_names, $cgc_name; # use CDS name as-is

      } elsif ($cgc_name =~ /clone:\S+/) { # this peptide doesn't match a CDS or gene - it is in an ORF somewhere on this clone
	push @CDS_names, $cgc_name; # store clone name with 'clone:' prefix - this will be dealt with specially

      } elsif ($cgc_name =~ /UNIPROT:(\S+)/) { # this is a UniProt ID
	my $uniprot = $1;
	if (exists $uniprot2wormpep->{$uniprot}) {
	  my $wormpep = $uniprot2wormpep->{$uniprot};
	  if (!exists $wormpep2cds->{$wormpep}) {
	    if (! exists $wormpep_history_proteins->{$wormpep}) {
	      print "line $line with wormpep ID '$wormpep' has a non-existent wormpep ID\n";
	      $non_existent++;
	      next;
	    } else {
	      push @CDS_names, $wormpep_history_proteins->{$wormpep}; # the CDS names for this UniProt
	    }
	  } else {
	    push @CDS_names, (split /\s+/, $wormpep2cds->{$wormpep}); # the CDS names for this UniProt
	  }
	} else {
	  print "line $line has a UniProt ID '$uniprot' that is not in Wormbase\n";
	  $not_in_wormbase++;
	}
	
      } else {
	die($input_file ,":gene name '$cgc_name' not recognised\n");
      }

      # construct the peptide data - these peptides map to these CDSs
      push @{$peptide->{$peptide_sequence}}, @CDS_names;

    }
  }
  close (DATA);

  print "UniProt entries that have a non existent WormPep ID: $non_existent\n";
  print "UniProt entries that are not in WormBase: $not_in_wormbase\n";

}


##########################################
# maps all the peptides that should map to this CDS
sub process_cds {
  # CDS_name is the original CDS name (without the isoform 'a' added)
  # used as the key to look for the peptides that should map to this
  # CDS 
  # 
  # CDS_name_isoform is the same as the CDS_name unless it has been
  # changed to an isoform and we want to try mapping to it in which
  # case it has a 'a' added to the name to make the first isoform CDS
  # name.
  my ($CDS_name, $CDS_name_isoform, $wormpep_history, $unique_clones, %proteins) = @_;

  my $not_found_peptides;  # href of peptides we couldn't map

  my $all_maps_to_genome_ok = 0;		# flag for mapping all peptides to the CDS genome position OK

  if (! defined $proteins{$CDS_name}) {
    print "process_cds() : Can't find any proteins for CDS $CDS_name\n";
    return ($all_maps_to_genome_ok, $not_found_peptides)
  }

  print "$CDS_name_isoform " . @{ $proteins{$CDS_name} } . "\n" if ($verbose);

  my @peptides_in_protein = @{ $proteins{$CDS_name} };

  # go through each of the historical wormpep proteins for this gene
  # looking for the most recent one that all the peptides map to
  # and get the version of that protein so we can then get the history CDS if it is not the latest version

  my ($wormpep_ids_aref, $versions_aref) = &get_previous_wormpep_ids($CDS_name_isoform, $wormpep_history);
  my @wormpep_ids = @{$wormpep_ids_aref};
  my @versions = @{$versions_aref};

  # try mapping the peptides to each version of the wormpep protein of
  # this CDS, starting with the most recent version
  my $latest_cds = 0;		# count of the previous revisions of genes that we have to go back through before finding a match
  foreach my $wormpep_id (@wormpep_ids) {

    my $wormpep_seq = $wormpep->{$wormpep_id};
    my ($all_peptides_mapped_ok, %protein_positions);
    # check to see if the protein exists!!!
    if (defined $wormpep_seq) {
      # see if the peptides all map to this version of the wormpep protein sequence
      ($all_peptides_mapped_ok, $not_found_peptides, %protein_positions) = &map_peptides_to_protein($wormpep_seq, @peptides_in_protein);
    } else {
      $all_peptides_mapped_ok = 0;
    }


    if ($all_peptides_mapped_ok) {
      print "We have found all the peptides in CDS: $CDS_name_isoform wormpep_id: $wormpep_id version=$latest_cds\n" if ($verbose);

      # result flag is true so far if we can map this to the genome
      $all_maps_to_genome_ok = 1;

      # see if we are using the latest version of this protein or not
      # if not, then construct the history ID name of the CDS to get

#      if (defined $versions[$latest_cds] && $latest_cds == 0) {
#	$log->write_to("*** $CDS_name_isoform no longer exists - it doesn't have a wormpep history, so can't map the peptides to the genome\n");
#	$all_maps_to_genome_ok = 0;
#      } elsif (defined $versions[$latest_cds] && $latest_cds != 0) {

#      if (defined $versions[$latest_cds]) {
#	$CDS_name_isoform = $CDS_name_isoform . ":${history_suffix}" . ($versions[$latest_cds] - 1);
#	print "$CDS_name_isoform is using history version $versions[$latest_cds]\n" if ($verbose);
#      }

      # now get the name of the CDS that made the protein $wormpep_id
      $CDS_name_isoform = &get_protein_details($wormpep_id);
      if (! defined $CDS_name_isoform) {
	# print "$CDS_name : the protein $wormpep_id has no Corresponding_CDS\n";
	next;
      }

      my ($clone, $cds_start, $cds_end, $exons_start_ref, $exons_end_ref);
      if ($all_maps_to_genome_ok) {
	# get the details for mapping to the genome
	# we can't do this outside of this $wormpep_id loop because we don't know until now whether we must use a history CDS
	($clone, $cds_start, $cds_end, $exons_start_ref, $exons_end_ref) = &get_cds_details($CDS_name_isoform);
	# if we couldn't get the details of the CDS, then we just have to skip writing the genome mapping data for this protein
	if (! defined $clone) {
	  #$log->write_to("Couldn't get the mapping details of the CDS $CDS_name_isoform\n");
	  $all_maps_to_genome_ok = 0;
	}
      }


      # only write this genome mapping stuff if we could find a CDS for the protein the peptides matched
      if ($all_maps_to_genome_ok) {
	# note which clones were used for genome mapping
	$unique_clones->{$clone} = 1;
      }

      # output the details for each of the peptides for this protein
      foreach my $ms_peptide (@peptides_in_protein) {
	  
	# only write this genome mapping stuff if we could find a CDS for the protein this peptide matched
	if ($all_maps_to_genome_ok) {
	  # get the details for mapping to the genome
	  my (@homol_lines) = &get_genome_mapping($ms_peptide, $CDS_name_isoform, $cds_start, $cds_end, $exons_start_ref, $exons_end_ref, \%protein_positions); # get part of the MSPeptide_homol line like "clone_start clone_end 1 peptide_length AlignPepDNA"
	  # write the Homol_data for the genome mapping
	  print OUT "\n";
	  print OUT "// Normal output for peptide $ms_peptide matching CDS: $CDS_name_isoform\n";
	  print OUT "Homol_data : \"$clone:Mass-spec\"\n";
	  foreach my $homol_line (@homol_lines) {
	    print OUT "MSPeptide_homol \"MSP:$ms_peptide\" mass_spec_genome 100.0 $homol_line\n";
	  }
	}

	# output the MS_peptide_results hash in the Mass_spec_peptide object
	print OUT "\n";
	print OUT "Mass_spec_peptide : \"MSP:$ms_peptide\"\n";
	foreach my $experiment_id (keys %experiment) {
	  print OUT "Mass_spec_experiments \"$experiment_id\" Protein_probability ", $experiment{$experiment_id}{'PROTEIN_PROBABILITY'}{$ms_peptide} ,"\n" if (defined $experiment{$experiment_id}{'PROTEIN_PROBABILITY'}{$ms_peptide});
	  print OUT "Mass_spec_experiments \"$experiment_id\" Peptide_probability ", $experiment{$experiment_id}{'PEPTIDE_PROBABILITY'}{$ms_peptide} ,"\n" if (defined $experiment{$experiment_id}{'PEPTIDE_PROBABILITY'}{$ms_peptide});
	  # check to see if this peptide has a match in this experiment
	  if (exists $experiment{$experiment_id}{HAS_PEPTIDE}{$ms_peptide}) {
	    # the protein that this peptide in this experiment matches
	    print OUT "Mass_spec_experiments \"$experiment_id\" Protein \"$wormpep_id\"\n"; 
	    # count the number of times this peptide maps to a protein in this experiment
	    $peptide_count{$ms_peptide}{$experiment_id}++;
	  }
	}

	# output the protein object for this peptide
	print OUT "\n";
	print OUT "Protein : \"MSP:$ms_peptide\"\n";
	print OUT "Peptide \"MSP:$ms_peptide\"\n";

	# output the homology of the wormpep protein to the peptide
	my $motif_out = "";	# holds the motif output, if any

	print OUT "\n";
	print OUT "Protein : \"$$wormpep_id\"\n";                                          
	my $pos = $protein_positions{$ms_peptide}; # position in protein
	my $len = length(&no_underscore($ms_peptide)); # length of peptide sequence
	my $end = $pos+$len-1;
	print OUT "Pep_homol \"MSP:$ms_peptide\" mass-spec 1 $pos $end 1 $len\n"; # note where the peptide maps to a wormpep protein
	# output the posttranslation modifications (if any) for this protein
	foreach my $experiment_id (keys %experiment) {
	  if (exists $experiment{$experiment_id}{HAS_PEPTIDE}{$ms_peptide}) {
	    if (exists $experiment{$experiment_id}{'posttranslation_modification_database'}) {
	      my $pt_db = $experiment{$experiment_id}{'posttranslation_modification_database'};
	      my $got_a_post_trans_modification_in_this_protein = 0;
	      my @pt_positions = &pos_underscore($ms_peptide);
	      foreach my $pt_pos (@pt_positions) {
		my $pt_pos_in_protein = $pt_pos + $pos - 1;
		print OUT "Motif_homol \"$pt_db:$CDS_name_isoform\" $pt_db 0 $pt_pos_in_protein $pt_pos_in_protein 1 1\n";
		$got_a_post_trans_modification_in_this_protein = 1;
	      }
	      
	      if ($got_a_post_trans_modification_in_this_protein) {
		# and write some details for the Motif Object
		$motif_out .= "\n";
		$motif_out .= "\n";
		$motif_out .= "Motif : \"$pt_db:$CDS_name_isoform\"\n";
		$motif_out .= "Title \"$experiment{$experiment_id}{'posttranslation_modification_type'}\"\n";
		$motif_out .= "Database \"$pt_db\" \"${pt_db}_ID\" \"$CDS_name_isoform\"\n";
		$motif_out .= "Pep_homol \"$wormpep_id\"\n";
		$motif_out .= "\n";
		$motif_out .= "\n";
	      }
	    
	    } elsif (exists $experiment{$experiment_id}{'posttranslation_modification_type'}) {
	      #print "*** Posttranslation_modification_type specified, but no Posttranslation_modification_database found\nShould add in output for the postranslational data to be a Feature in the Protein? A bit like this maybe?:\n";
	      my $pt_type = $experiment{$experiment_id}{'posttranslation_modification_type'};
	      my @pt_positions = &pos_underscore($ms_peptide);
	      foreach my $pt_pos (@pt_positions) {
		my $pt_pos_in_protein = $pt_pos + $pos - 1;
		print OUT "Feature \"$pt_type\" $pt_pos_in_protein $pt_pos_in_protein 0\n";
		print OUT "Remark \"$pt_type defined by mass-spec peptide MSP:$ms_peptide from mass-spec experiment $experiment_id\"\n";
	      }
	    }
	  }
	}

	# now output the motif objects, if any
	print OUT $motif_out;
	$motif_out = "";

	# note that this peptide has been used so that we need to output it when only doing a small debugging output
	$peptides_used{$ms_peptide} = 1;
      }     # foreach my $ms_peptide 
      # don't bother to look at any more wormpep versions, we have mapped it OK
      last;
    }				# if ($all_peptides_mapped_ok) 

    # count the version back we have to go
    $latest_cds++;

  }

  return ($all_maps_to_genome_ok, $not_found_peptides);
}

##########################################
# For when we have peptides to map to the ORFs of a clone
# We can map to novel introns in the clone by specifying the positions
# around the intron as 'clone:AC3<end_exon1,start_exon2>'

sub process_clones {

  # CDS_name is the name of the clone with 'clone:' prefixed
  # $proteins{$CDS_name} contains the peptides that map to an ORF on this clone

  my ($CDS_name, $unique_clones, %proteins) = @_;

  my $not_found_peptides;   # href of peptides we couldn't map

  my @peptides_in_protein = @{ $proteins{$CDS_name} };

  my $all_maps_to_genome_ok = 0;		# flag for mapping all peptides to the CDS genome position OK

  # try mapping the peptides to the clone ORFs
  my ($clone) = ($CDS_name =~ /^clone:(\S+)/);

  my $clone_seq;
  my $have_splices = 0;

  # construct the (imaginary) exon sources
  my @exons_start;
  my @exons_end;
  my $length_orig_clone_seq;
  
  # does this clone name have splice sites?
  # if so, then make the spliced clone sequence
  my $clone_desc = $clone;
  if ($clone =~ /(\S+?)<([\d,]+)>/) {
    $clone = $1;
    my $splices = $2;
    my @splices = split /,/, $splices;
    $have_splices = 1;
    $clone_seq = &get_clone_sequence($clone); # get the clone sequence
    $length_orig_clone_seq = length $clone_seq;
    # splice the clone sequence
    my $spliced_clone_seq="";
    my $next_end = length $clone_seq;

    @exons_start = (1); # imaginary start of first exon, assume peptide matches frame 1 for now
    while (@splices) {
      # get the array values from the end of @splices two at a time
      my ($one, $two) = splice(@splices, -2);
      # check that the splice uses GT..AG splices or GC..AG splices
      my $one_spliceseq = uc substr($clone_seq, $one, 2);
      my $two_spliceseq = uc substr($clone_seq, $two-3, 2);
      if (! ((($one_spliceseq eq 'GT' || $one_spliceseq eq 'GC') && $two_spliceseq eq 'AG') || 
	     ($one_spliceseq eq 'CT' && ($two_spliceseq eq 'AC' || $two_spliceseq eq 'GC'))) ) {
#      if (! (($one_spliceseq eq 'GT' && $two_spliceseq eq 'AG') || 
#	     ($one_spliceseq eq 'CT' && $two_spliceseq eq 'AC')) ) {
	$log->write_to("This doesn't look like a canonical splice site\n");
	$log->write_to("$clone_desc ($one $one_spliceseq..$two_spliceseq $two)\n");
      }
      $spliced_clone_seq = substr($clone_seq, $two-1, $next_end-$two+1) . $spliced_clone_seq;
      $next_end = $one;

      # make the imaginary exon sources
      push @exons_start, $two;
      push @exons_end, $one;

    }
    if ($next_end > -1) {$spliced_clone_seq = substr($clone_seq, 0, $next_end) . $spliced_clone_seq;}
    push @exons_end, length($clone_seq); # imaginary end of last exon
    $clone_seq = $spliced_clone_seq;

    # sort the imaginary exons numerically
    @exons_start = sort {$a <=> $b} @exons_start;
    @exons_end = sort {$a <=> $b} @exons_end;

  } else {
    # normal clone sequence with no splicing
    $clone_seq = &get_clone_sequence($clone); # get the clone sequence
    $length_orig_clone_seq = length $clone_seq;
    @exons_start = (1);
    @exons_end = (length($clone_seq));
  }



  print "\n\nnext name name=$CDS_name" if ($verbose);

  my @wormpep_seq = &translate_clone($clone_seq); # translate the clone sequence in all six frames
  #print "clone translation = @wormpep_seq" if ($verbose);

  # note which clones were used for genome mapping
  $unique_clones->{$clone} = 1;
      

  my ($this_peptide_mapped_ok, %protein_positions);
  my $mapped_in_frame;

  # output the details for each of the peptides for this protein
  foreach my $ms_peptide (@peptides_in_protein) {

    # check to see if the protein hits an ORF in the translated clone sequence - check all 6 frames
    foreach my $frame (1..6) {
      my $wormpep_seq = $wormpep_seq[$frame];
      # see if this peptide maps to this frame of the clone translation
      ($this_peptide_mapped_ok, $not_found_peptides, %protein_positions) = &map_peptides_to_protein($wormpep_seq, ($ms_peptide));
      if ($this_peptide_mapped_ok) {
	$mapped_in_frame = $frame;
	$log->write_to("The peptide $ms_peptide matched in frame $frame OK in $clone_desc\n") if ($verbose);
	last;
      }
    }

    # only write this genome mapping stuff if it mapped to a clone translation frame
    if ($this_peptide_mapped_ok) {
      
      # get the details for mapping to the genome
      # get part of the MSPeptide_homol line like "clone_start clone_end 1 peptide_length AlignPepDNA"
      my @homol_lines;

      print "start exons @exons_start\n" if ($verbose);
      print "end exons @exons_end\n" if ($verbose);
      if ($mapped_in_frame < 4) { # forward sense translation of the clone
	# adjust the exon_starts and _ends
	foreach my $exon_start (@exons_start) { $exon_start = $exon_start - $mapped_in_frame + 1 }
	foreach my $exon_end (@exons_end) { $exon_end = $exon_end - $mapped_in_frame + 1 }
	$exons_start[0] = 1;

	@homol_lines = &get_genome_mapping($ms_peptide, $CDS_name, $mapped_in_frame, length($clone_seq), \@exons_start, \@exons_end, \%protein_positions); 

      } else { # reverse sense translation of the clone

	# adjust the exon_starts and _ends
	my $cds_start = $length_orig_clone_seq - ($mapped_in_frame - 3) + 1;
	my $cds_end = 1;
	# swap over the exon_start and exons_end arrays and use the offset from the end of the clone
	my @exons_start_new;
	my @exons_end_new;
	while (@exons_start) {
	  my $start = shift @exons_start;
	  my $end = shift @exons_end;
	  push @exons_start_new, $cds_start - $end + 1;
	  push @exons_end_new,  $cds_start - $start + 1;
	}
	@exons_start = sort {$a <=> $b} @exons_start_new;
	@exons_end = sort {$a <=> $b} @exons_end_new;
	$exons_start[0] = 1; # just to be certain

	print "after fixing to look like a CDS start exons @exons_start\n" if ($verbose);
	print "after fixing to look like a CDS exon end exons @exons_end\n" if ($verbose);

	@homol_lines = &get_genome_mapping($ms_peptide, $CDS_name, $cds_start, $cds_end, \@exons_start, \@exons_end, \%protein_positions); 
      }
      # write the Homol_data for the genome mapping
      print OUT "\n";
      print OUT "// Clone-positioned output for peptide $ms_peptide matching Clone: $clone_desc\n";
      print OUT "Homol_data : \"$clone:Mass-spec\"\n";
      foreach my $homol_line (@homol_lines) {
	print OUT "MSPeptide_homol \"MSP:$ms_peptide\" mass_spec_genome 100.0 $homol_line\n";
      }
    
      # output the MS_peptide_results hash in the Mass_spec_peptide object
      print OUT "\n";
      print OUT "// ms-peptide stuff for $CDS_name\n";
      print OUT "Mass_spec_peptide : \"MSP:$ms_peptide\"\n";
      foreach my $experiment_id (keys %experiment) {
	print OUT "Mass_spec_experiments \"$experiment_id\" Protein_probability ", $experiment{$experiment_id}{'PROTEIN_PROBABILITY'}{$ms_peptide} ,"\n" if (defined $experiment{$experiment_id}{'PROTEIN_PROBABILITY'}{$ms_peptide});
	print OUT "Mass_spec_experiments \"$experiment_id\" Peptide_probability ", $experiment{$experiment_id}{'PEPTIDE_PROBABILITY'}{$ms_peptide} ,"\n" if (defined $experiment{$experiment_id}{'PEPTIDE_PROBABILITY'}{$ms_peptide});
	# check to see if this peptide has a match in this experiment
	if (exists $experiment{$experiment_id}{HAS_PEPTIDE}{$ms_peptide}) {
	  print OUT "Mass_spec_experiments \"$experiment_id\"\n";
          # count the number of times this peptide maps to a protein in this experiment
	  $peptide_count{$ms_peptide}{$experiment_id}++;
	}
      }
      
      # note that this peptide has been used so that we need to output it when only doing a small debugging output
      $peptides_used{$ms_peptide} = 1;
    } else {
      $log->write_to("The peptide $ms_peptide didn't match in $clone_desc\n");
    }
  }     # foreach my $ms_peptide 

  return (1, $not_found_peptides);

}

##########################################
# take a protein and a set of peptides and map the peptides onto the protein
# returns a hash key=peptide sequence, value=position in protein

sub map_peptides_to_protein {
  my ($protein, @peptides) = @_;
  my %position;
  my $ok = 1;			# true if all the peptides map to this protein
  my %not_found_peptides;
  
  foreach my $peptide (@peptides) {
    my $no_underscore_peptide = &no_underscore($peptide); # remove the post-translation modification marker
    #print "mapping $no_underscore_peptide\nto $protein\n";
    print "mapping $no_underscore_peptide\n" if ($verbose);
    if (!defined $no_underscore_peptide) {print "no_underscore_peptide not defined for $peptide";}
    if (!defined $protein) {print "protein not defined for $peptide";}
    my $pos = index($protein, $no_underscore_peptide);
    print "mapping result: $pos\n" if ($verbose);
    $position{$peptide} = $pos+1; # convert from offset=0 to offset=1;
    if ($pos == -1) {
      $ok = 0;	# set the flag if a peptide is not found in the protein
      $not_found_peptides{$no_underscore_peptide} = 1;
    }
    
  }

  return ($ok, \%not_found_peptides, %position);
}


##########################################
# get the wormpep protein IDs used by previous version of the gene
#my %wormpep_history = &get_wormpep_history;

sub get_wormpep_history {
  # get database version file
  my $data_file = $wormbase->basedir . "/WORMPEP/${pepdir_prefix}pep${database_version}/${pepdir_prefix}pep.history$database_version";

  my %wormpep_history;

  open (HIST, "< $data_file") || die "Can't open $data_file\n";
  while (my $line = <HIST>) {
    chomp $line;
    my @f = split /\s+/, $line;
    # $cds_id, $wormpep_id, $version1, $version2 <- ($version2 is undef if current version)
    my $cds_id = shift @f;
    push @{$wormpep_history{$cds_id}}, [@f]; # I didn't mean to push the array of array, but we are stuck with this now, so leave it
  }
  close(HIST);


  return \%wormpep_history;

}
##########################################
# get the hash of wormpep sequence names keyed by old wormpep protein IDs

sub get_wormpep_history_proteins {
  my ($wormpep_history) = @_;

  my %wormpep_history_proteins;

  foreach my $cds_id (keys %{$wormpep_history}) {
    my @f = $wormpep_history->{$cds_id};
    my @g = @{$f[0][0]};
#    if (defined $g[2] && $g[2] ne "") { # if this is an old version of the CDS, give the 'history' sequence name
#      $cds_id .= ":${history_suffix}" . ($g[2] - 1); 
#    }
    $wormpep_history_proteins{$g[0]} = $cds_id; # key = wormpep protein ID, value = cds_id
  }


  return \%wormpep_history_proteins;

}


##########################################
# returns true if the CDS has a current history
# returns false if the CDS has been made into an isoform (or pseudogene, etc.)
sub has_current_history {
  my ($history, $CDS_name) = @_;

  # if the CDS_name is not a valid key then this CDS_name has been
  # made into a isoform (or pseudogene, etc.) and no details are known
  # of previous wormpep IDs
  if (!exists $history->{$CDS_name}) {return 0;}

  # if the last line in the history file for this CDS_name has only
  # one value, then there is still a current CDS of this name
  my @cds_history = @{$history->{$CDS_name}};
  my $last_result = pop @cds_history;
  print "@{$last_result}\n" if ($verbose);
  if (scalar @{$last_result} == 2) {return 1;} 

  # so there isn't a current CDS of this name
  return 0;
}

##########################################
# get all previous wormpep protein IDs for a given CDS
# returned in the order in which they occur in the data file
# which seems to be in historical order (oldest first)
# @worpep_ids = get_previous_wormpep_ids($gene_name);

sub get_previous_wormpep_ids {

  my ($CDS_name, $history) = @_;

  my @wormpep_ids;
  my @versions;

  if (!defined $CDS_name) {$log->log_and_die("\nNot defined CDS_name\n")}
  if (!exists $history->{$CDS_name}) {
    $log->write_to("\nIn get_previous_wormpep_ids(): history->{$CDS_name} does not exist\n");
    return \@wormpep_ids, \@versions;
  }

  my @cds_history = @{$history->{$CDS_name}};
  foreach my $version (@cds_history) {
    my ($wormpep_id, $version1, $version2) = @{$version};
    push @wormpep_ids, $wormpep_id; # the wormpep ID
    if (defined $version2) {
      push @versions, $version2; # the second (ending) Build release version, 
    } else {
      push @versions, undef;	# undef if this is still the current version
    }
  }

  # get the most recent versions first
  @wormpep_ids = reverse @wormpep_ids;
  @versions = reverse @versions;

  return \@wormpep_ids, \@versions;
}

##########################################
# get all the wormpep proteins from the wormpep fasta file by wormpep ID
# 

sub get_wormpep_by_wp_id {

  my $data_file = $wormbase->basedir . "/WORMPEP/${pepdir_prefix}pep${database_version}/${pepdir_prefix}pep.fasta${database_version}";
  my $seq="";
  my $id="";
  my %wormpep;

  open (SEQ, "<$data_file") || die "Can't open $data_file\n";
  while (my $line = <SEQ>) {
    chomp $line;
    if ($line =~ />(\S+)/) {
      if ($seq ne "") {$wormpep{$id} = $seq;} # store sequence
      $id = $1;
      $seq = "";
    } else {
      $seq .= $line;
    }
  }
  close (SEQ);
  $wormpep{$id} = $seq;		# store last sequence

  return \%wormpep;
}

##########################################
# my ($clone, $cds_start, $exons_start_ref, $exons_end_ref) = &get_cds_details($CDS_name);
# get the clone, clone start position and exons list from a CDS
#

sub get_cds_details {
  my ($CDS_name) = @_;

  print "Ace fetch->" if ($verbose);
  print "$CDS_name\n" if ($verbose);
  my $cds_obj = $db->fetch("CDS" => $CDS_name);
  # debug
  if (! defined $cds_obj) {
    print "Can't fetch CDS object for $CDS_name\n";
    return (undef, undef, undef, undef, undef);
  }

  my $clone = $cds_obj->Sequence;
  if (! defined $clone) {
    print "Can't fetch clone for $CDS_name\n";
    return (undef, undef, undef, undef, undef);
  }

  # get the named object
  # then get the data following a tag in that object
  # and the NEXT data following the tag
  my @exons_start = $cds_obj->Source_exons;
  my @exons_end = $cds_obj->Source_exons(2);
  $cds_obj->DESTROY();
  if (! @exons_start || ! @exons_end) {
    print "Can't fetch exons for $CDS_name\n";
    return (undef, undef, undef, undef, undef);
  }

  # get the start position of the CDS in the clone's SMap
  print "Ace fetch->clone = $clone\n" if ($verbose);
  my $clone_obj = $db->fetch("Sequence" => $clone);
  if (! defined $clone_obj) {
    print "Can't fetch clone object for $clone\n";
    return (undef, undef, undef, undef, undef);
  }

  my $quoted_obj = $CDS_name;
  $quoted_obj =~ s/\./\\./g;            # replace . with \.
  my ($cds, $cds_start, $cds_end) = $clone_obj->at("SMap.S_child.CDS_child.$quoted_obj")->row;
  my $foundit = 1;

  if (! defined $cds_start || ! defined $cds_end) {
    print "Can't fetch start/end positions for $CDS_name in $clone\n";
    return (undef, undef, undef, undef, undef);
  }

  print "Found $CDS_name in $clone:$cds_start..$cds_end\n" if ($verbose);

  # if the CDS is on the reverse sense, then $cds_start > $cds_end
  return ($clone, $cds_start, $cds_end, \@exons_start, \@exons_end);
}

##########################################
# get the name of the CDS that made the protein we have found

sub get_protein_details {

  my ($protein_id) = @_;

  print "Ace protein fetch->" if ($verbose);
  print "$protein_id\n" if ($verbose);

  my $protein_obj = $db->fetch("Protein" => $protein_id);

  if (! defined $protein_obj) {
    print "Can't fetch Protein object for $protein_id\n";
    return undef;
  }

  my $CDS = $protein_obj->Corresponding_CDS;

  if (defined $CDS) {
    return $CDS->name;
  } else {
    return undef;
  }

}

##########################################
# get the mapping of UniProt names to protein names
#  my %uniprot2wormpep = get_uniprot();

sub get_uniprot {

  print "Reading UniProt data\n";

  my %uniprot;
#  my $database = $wormbase->database('current');
  my $database = $wormbase->autoace;

  my $table_def = &write_uniprot_def;
  my $table_query = $wormbase->table_maker_query($database, $table_def);
  while(<$table_query>) {
    chomp;
    s/\"//g;  #remove "
    next if (/acedb/ or /\/\//);
    my @data = split(/\t/,$_);
    my ($wormpep, $uniprot) = @data;
    
    if (!defined $uniprot) {next;}
    
    $uniprot{$uniprot} = $wormpep;
    #print "$uniprot\t$wormpep\n";
  }
  
  return %uniprot;
}
############################################################################
# this will write out an acedb tablemaker defn to a temp file
############################################################################
#  my %uniprot = get_uniprot($database);

sub write_uniprot_def {
  my $def = "/tmp/Uniprot_$$.def";
  open TMP,">$def" or $log->log_and_die("cant write $def: $!\n");

  my $txt = <<END;

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Protein 
From 1 
Condition IS "${pep_prefix}*"

Colonne 2 
Width 12 
Mandatory
Hidden 
Class 
Class Database 
From 1 
Tag Database 
Condition UniProt
 
Colonne 3 
Width 12 
Mandatory
Hidden 
Class 
Class Database_field 
Right_of 2 
Tag  HERE  
Condition UniProtAcc
 
Colonne 4 
Width 30 
Mandatory
Visible 
Class 
Class Text 
Right_of 3 
Tag  HERE  
 
 
END

  print TMP $txt;
  close TMP;
  return $def;
}


############################################################################

##########################################
# my @homol_lines = &get_genome_mapping($ms_peptide, $CDS_name, $cds_start, $exons_start_ref, $exons_end_ref, \%protein_positions); 
# get part of the MSPeptide_homol line like "clone_start clone_end 1 peptide_length AlignPepDNA"
#
# take the information describing the position of the CDS's exons in
# the clone and the peptides position in the CDS and construct
# required part of the MSPeptide_homol data line for output.
#
# Check that the peptide maps to the genome correctly by translating the mapped part

sub get_genome_mapping {
  my ($ms_peptide, $CDS_name, $cds_start, $cds_end, $exons_start_ref, $exons_end_ref, $protein_positions_href) = @_;
  print "In get_genome_mapping(): $ms_peptide, $CDS_name, $cds_start, $cds_end, @$exons_start_ref, @$exons_end_ref, ",%{$protein_positions_href},"\n" if ($verbose);
  my %protein_positions = %{$protein_positions_href};

  # the positions of the exons on the clone relative to the start of the CDS
  my @exons_start = @{$exons_start_ref};
  my @exons_end = @{$exons_end_ref};

  my @output;			# returned list of output lines

  # the position in the CDS of the start and end of exons
  my @exons_internal_start = ();
  my @exons_internal_end = ();

  # get the peptide length
  my $peptide_length = length(&no_underscore($ms_peptide));

  print "ms_peptide = $ms_peptide CDS_name = $CDS_name\n" if ($verbose);

  # get the start and end position of the peptide in the protein
  my $peptide_cds_start = $protein_positions{$ms_peptide}; 
  print "position of peptide in protein = $peptide_cds_start\n" if ($verbose);
  my $peptide_cds_end = $peptide_cds_start + $peptide_length - 1;
  # convert protein positions to CDS coding positions
  $peptide_cds_start = ($peptide_cds_start * 3) - 2; # start at the beginning of the codon, not the end :-)
  $peptide_cds_end *= 3;
  print "peptide_cds_start = $peptide_cds_start\n" if ($verbose);

  # find the positions in the CDS of the start and end of the exons (i.e. the splice sites on the CDS)
  my $exon_count = 0;
  my $prev_end = 0;
  print "number of exons = ",scalar @exons_start," and ",scalar @exons_end,"\n" if ($verbose);
  foreach my $exon_start (@exons_start) {
    print "exon_start $exon_start exon_end $exons_end[$exon_count]\n" if ($verbose);
    my $exon_length = $exons_end[$exon_count] - $exon_start + 1;
    print "exon_length $exon_length\n" if ($verbose);
    my $start = $prev_end + 1;
    my $end = $start + $exon_length - 1;
    print "exon internal positions on the CDS mRNA sequence :start $start end $end\n" if ($verbose);
    push @exons_internal_start, $start;
    push @exons_internal_end, $end;
    $prev_end = $end;
    $exon_count++;
  }
  print "number of internal CDS mRNA coords exons = ",scalar @exons_internal_start," and ",scalar @exons_internal_end,"\n" if ($verbose);


  # @exons_start/end = positions of the exons start/end in the transcript, from 1 to length of transcript (i.e. length of CDS plus length of introns)
  # @exons_internal_start/end = positions in the cDNA of the start and end of the exons (i.e. the splice sites on the cDNA)
  # $peptide_cds_start/end = start/end positions of the peptide in the CDS
  # $peptide_length = length of peptide in aa's

  # go through the exons seeing if the peptide is on this exon 
  my $exon_clone_start;		# start position of exon on clone
  my $exon_clone_end;
  my $peptide_clone_start;	# start position of peptide on clone
  my $peptide_clone_end;
  my $peptide_start;		# start of peptide in aa's in this exon
  my $peptide_end;

  foreach $exon_count (0..(scalar @exons_start) - 1) {

    #print "exon_count = $exon_count\n" if ($verbose);
    #print "peptide_cds_start = $peptide_cds_start\n" if ($verbose);
    #print "peptide_cds_end = $peptide_cds_end\n" if ($verbose);
    #print "exon CDS boundaries = $exons_internal_start[$exon_count] $exons_internal_end[$exon_count]\n" if ($verbose);

    #
    # see if start of peptide in this exon
    #
    if ($peptide_cds_start >= $exons_internal_start[$exon_count] && $peptide_cds_start <= $exons_internal_end[$exon_count]) {
      print "start of peptide\n" if ($verbose);

      if ($cds_start < $cds_end) {
	$exon_clone_start = $cds_start + $exons_start[$exon_count] - 1; # get the clone position of the start of the exon
	$exon_clone_end   = $cds_start + $exons_end[$exon_count] - 1; # get the clone position of the end of the exon
      } else {			# CDS is in reverse sense
	$exon_clone_start = $cds_start - ($exons_start[$exon_count] - 1); # get the clone position of the start of the exon in reverse sense
	$exon_clone_end   = $cds_start - ($exons_end[$exon_count] - 1); # get the clone position of the end of the exon in reverse sense
      }

      # get the default start and end of the peptide in aa's in this exon
      $peptide_start = 1;
      $peptide_end = $peptide_length;

      # get the default start and end of the peptide on the clone
      my $bases_from_start_of_exon_to_peptide_start = $peptide_cds_start - $exons_internal_start[$exon_count];
      my $bases_from_start_of_exon_to_peptide_end = $peptide_cds_end - $exons_internal_start[$exon_count];
      if ($cds_start < $cds_end) {
	$peptide_clone_start = $exon_clone_start + $bases_from_start_of_exon_to_peptide_start;
	$peptide_clone_end = $exon_clone_start + $bases_from_start_of_exon_to_peptide_end;
      } else {			# CDS is in reverse sense
	$peptide_clone_start = $exon_clone_start - $bases_from_start_of_exon_to_peptide_start;
	$peptide_clone_end = $exon_clone_start - $bases_from_start_of_exon_to_peptide_end;
      }

      # see if the end of the peptide is after this exon
      if ($peptide_cds_end > $exons_internal_end[$exon_count]) {
	$peptide_clone_end = $exon_clone_end; # set the end of this bit of the peptide on the clone as being the position of the end of the exon
	$peptide_end = int (($exons_internal_end[$exon_count] - $peptide_cds_start) / 3) + 1;
      }

      push @output, "$peptide_clone_start $peptide_clone_end $peptide_start $peptide_end";
      print "output = $peptide_clone_start $peptide_clone_end $peptide_start $peptide_end\n" if ($verbose);

      if ($peptide_cds_end <= $exons_internal_end[$exon_count]) {last;}		# the peptide only covers one exon so we can get out of the loop

    #  
    # see if end of peptide in this exon
    #
    } elsif ($peptide_cds_end >= $exons_internal_start[$exon_count] && $peptide_cds_end <= $exons_internal_end[$exon_count]) {
      print "end of peptide\n" if ($verbose);

      if ($cds_start < $cds_end) {
	$exon_clone_start = $cds_start + $exons_start[$exon_count] - 1; # get the clone position of the start of the exon
	$exon_clone_end   = $cds_start + $exons_end[$exon_count] - 1; # get the clone position of the end of the exon
      } else {			# CDS is in reverse sense
	$exon_clone_start = $cds_start - ($exons_start[$exon_count] - 1); # get the clone position of the start of the exon in reverse sense
	$exon_clone_end   = $cds_start - ($exons_end[$exon_count] - 1); # get the clone position of the end of the exon in reverse sense
      }

      # get the start and end of the peptide in aa's in this exon
      $peptide_start = int (($exons_internal_start[$exon_count] - $peptide_cds_start + 2) / 3);
      $peptide_end = $peptide_length;
      $peptide_clone_start = $exon_clone_start;

      my $bases_from_start_of_exon_to_peptide_end = $peptide_cds_end - $exons_internal_start[$exon_count];
      if ($cds_start < $cds_end) {
	$peptide_clone_end = $exon_clone_start + $bases_from_start_of_exon_to_peptide_end;
      } else {
	$peptide_clone_end = $exon_clone_start - $bases_from_start_of_exon_to_peptide_end;
      }

      push @output, "$peptide_clone_start $peptide_clone_end $peptide_start $peptide_end";
      print "output = $peptide_clone_start $peptide_clone_end $peptide_start $peptide_end\n" if ($verbose);

      last;		# we have finished, so get out of the loop



    #
    # see if internal part of peptide in this exon
    #
    } elsif ($peptide_cds_start <= $exons_internal_start[$exon_count]) {
      print "internal to peptide\n" if ($verbose);
      
      if ($cds_start < $cds_end) {
	$exon_clone_start = $cds_start + $exons_start[$exon_count] - 1; # get the clone position of the start of the exon
	$exon_clone_end   = $cds_start + $exons_end[$exon_count] - 1; # get the clone position of the end of the exon
      } else {			# CDS is in reverse sense
	$exon_clone_start = $cds_start - ($exons_start[$exon_count] - 1); # get the clone position of the start of the exon in reverse sense
	$exon_clone_end   = $cds_start - ($exons_end[$exon_count] - 1); # get the clone position of the end of the exon in reverse sense
      }

      # get the start and end of the peptide in aa's in this exon
      $peptide_start = int (($exons_internal_start[$exon_count] - $peptide_cds_start + 2) / 3);
      $peptide_end = int (($exons_internal_end[$exon_count] - $peptide_cds_start) / 3) + 1;
      $peptide_clone_start = $exon_clone_start;
      $peptide_clone_end = $exon_clone_end;

      push @output, "$peptide_clone_start $peptide_clone_end $peptide_start $peptide_end";
      print "output = $peptide_clone_start $peptide_clone_end $peptide_start $peptide_end\n" if ($verbose);

    }

  }



  return @output;
}
##########################################
# returns the DNA sequence of a clone

sub get_clone_sequence {

  my ($clone) = @_;
  my $clone_obj = $db->fetch(Sequence => $clone);
  if (!defined $clone_obj) {$log->log_and_die("ERROR - have an undefined clone_obj for clone '$clone' in get_clone_sequence()\n");}
  my $dna = $clone_obj->asDNA();
  $dna =~ s/\>(\w+)\n//;	# remove title line
  $dna =~ s/\n//g;	        # remove newline chars
  return $dna;

}
##########################################


# returns the clone sequence translated in all six frames
sub translate_clone {
  my ($clone_seq) = @_;
  my @clone_pep;
  my %trans = (
	       'ttt' => 'F',
	       'ttc' => 'F',
	       'tta' => 'L',
	       'ttg' => 'L',

	       'tct' => 'S',
	       'tcc' => 'S',
	       'tca' => 'S',
	       'tcg' => 'S',

	       'tat' => 'Y',
	       'tac' => 'Y',
	       'taa' => '*',
	       'tag' => '*',

	       'tgt' => 'C',
	       'tgc' => 'C',
	       'tga' => '*',
	       'tgg' => 'W',
#
	       'ctt' => 'L',
	       'ctc' => 'L',
	       'cta' => 'L',
	       'ctg' => 'L',

	       'cct' => 'P',
	       'ccc' => 'P',
	       'cca' => 'P',
	       'ccg' => 'P',

	       'cat' => 'H',
	       'cac' => 'H',
	       'caa' => 'Q',
	       'cag' => 'Q',

	       'cgt' => 'R',
	       'cgc' => 'R',
	       'cga' => 'R',
	       'cgg' => 'R',
#
	       'att' => 'I',
	       'atc' => 'I',
	       'ata' => 'I',
	       'atg' => 'M',

	       'act' => 'T',
	       'acc' => 'T',
	       'aca' => 'T',
	       'acg' => 'T',

	       'aat' => 'N',
	       'aac' => 'N',
	       'aaa' => 'K',
	       'aag' => 'K',

	       'agt' => 'S',
	       'agc' => 'S',
	       'aga' => 'R',
	       'agg' => 'R',
#
	       'gtt' => 'V',
	       'gtc' => 'V',
	       'gta' => 'V',
	       'gtg' => 'V',

	       'gct' => 'A',
	       'gcc' => 'A',
	       'gca' => 'A',
	       'gcg' => 'A',

	       'gat' => 'D',
	       'gac' => 'D',
	       'gaa' => 'E',
	       'gag' => 'E',

	       'ggt' => 'G',
	       'ggc' => 'G',
	       'gga' => 'G',
	       'ggg' => 'G',
	      );



  foreach my $frame (1,2,3) {
    my $translation=();
    for (my $i=$frame-1; $i<length($clone_seq)-2; $i+=3) {
      my $c = substr($clone_seq, $i, 3);
      $translation .= $trans{$c};
    }
    $clone_pep[$frame] = $translation;
  }
  $clone_seq = &DNA_revcomp($clone_seq);
  foreach my $frame (1,2,3) {
    my $translation=();
    for (my $i=$frame-1; $i<length($clone_seq)-2; $i+=3) {
      my $c = substr($clone_seq, $i, 3);
      $translation .= $trans{$c};
    }
    $clone_pep[$frame+3] = $translation;
  }

  return @clone_pep;

}
##########################################

=head2

  Title   :   DNA_revcomp
  Usage   :   my $revcomp = $seq_obj->($seq)
  Function:   revcomp DNA seq
  Returns :   DNA sequence as string
  Args    :   DNA sequence as string

=cut


sub DNA_revcomp
  {
    my $revseq = reverse shift;
    $revseq =~ tr/a/x/;
    $revseq =~ tr/t/a/;
    $revseq =~ tr/x/t/;
    $revseq =~ tr/g/x/;
    $revseq =~ tr/c/g/;
    $revseq =~ tr/x/c/;
    return ($revseq);
  }


##########################################
# returns the peptide sequence with any underscore characters
# (post-translational modification) stripped out

sub no_underscore {
  my $seq = shift;
  $seq =~ s/_//g;
  return $seq;
}

##########################################
# returns the positions of any underscore characters
# (post-translational modification)
# so 'S_AS_' will give the result (1, 3)

sub pos_underscore {
  my $seq = shift;
  my $pos = -1;
  my @result;
  while (($pos = index($seq, '_', $pos)) > -1) {
    # adjust the position to allow for the number of other underscores
    # before this one
    push @result, $pos - @result; 
    $pos++;
  }
  return @result;
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




# Add perl documentation in POD format
# This should expand on your brief description above and 
# add details of any options that can be used with the program.  
# Such documentation can be viewed using the perldoc command.


__END__

=pod

=head2 NAME - make_mass_spec_ace.pl

=head1 USAGE

=over 4

=item make_mass_spec_ace.pl  [-options]

=back

This script writes the ace file for a set of mass spec data.

It reads in a file with columns consisting of the following data

1) matching CDS sequence name or gene WGN (CGC) name (or for those
   peptides which are found by searching aganist translated ORFs and
   don't match an existing gene, 'clone:'clone-name)
2) mass-spec peptide sequence
3) protein probability (if known, otherwise leave as blank or '.')
4) peptide probability (if known, otherwise leave as blank or '.')
5) experiment ID (generally the experimentor's initials and a number to make this unique)

The output file may have to be edited to add further details of the experiment
after the line "Mass_spec_experiment : "

script_template.pl MANDATORY arguments:

=over 4

=item -input file to read the mass-spec data from

=back

=over 4

=item -output file to write ace results to.

=back



=head1 script_template.pl  OPTIONAL arguments:

=over 4

=item -database database used to map mass spec peptides to, default is currentDB

=back

=over 4

=item -natural 

=back

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


=head1 REQUIREMENTS

=over 4

=item None at present.

=back

=head1 AUTHOR

=over 4

=item Gary Williams

=back

=cut
