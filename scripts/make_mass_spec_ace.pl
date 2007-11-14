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
# Last updated on: $Date: 2007-11-14 16:12:13 $      

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

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($input, $output, $database, $experiment_id, $natural);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "input:s"    => \$input,
	    "output:s"   => \$output,
	    "database:s" => \$database,
	    "experiment_id:s" => \$experiment_id,
	    "natural"    => \$natural,
	    );

# always in test mode
$test = 1;

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

# check parameters
if (! defined $database) {$database = $wormbase->database('current')}
if (! defined $input) {die "-input file not specified\n";}
if (! defined $output) {die "-output file not specified\n";}
if (! defined $experiment_id) {die "-experiment_id not specified\n";}

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


# open an ACE connection to parse details for mapping to genome
my $tace            = $wormbase->tace;        # TACE PATH
print "Connecting to Ace\n";
my $db = Ace->connect (-path => $database,
                       -program => $tace) || die "cannot connect to database at $database\n";



# this gives a hash with key=experiment ID value=list of details
my ($experiment_hashref, $peptide_hashref) = parse_data();

my %experiment = %{$experiment_hashref};
my %peptide = %{$peptide_hashref};
my $peptide_id_count=1;
my %peptides_used;		# quick record of the peptides that have been used, for debugging to produce small test outputs
my %unique_peptides;		# non-redundant set of peptides used in all the experiments
my %unique_clones;  		# non-redundant set of clones mapped to in all the experiments
			

# open the output ACE file
open (OUT, "> $output") || die "Can't open output file $output\n";
print OUT "\n";			# blank line at start just in case we concatenate this to another ace file

# get the set of CDS, history versions and resulting protein IDs
my $protein_history_aref = &get_protein_history;

# get the peptides used in each protein
my %proteins;
foreach my $experiment_id (keys %experiment) {

  my %peptide_count = ();		# hash to hold count of times a peptide maps to proteins in this experiment
  print "experiment_id = $experiment_id\n";

  my @peptides = @{ $experiment{$experiment_id}->{'PEPTIDES'} };
  # get the unique non-redundant set of peptides in all experiments
  foreach my $pep (@peptides) {
    $unique_peptides{$pep} = 1;
  }

  # get the proteins and their peptides
  foreach my $peptide_sequence (@peptides) {
    my @protein_names = @{ $peptide{$peptide_sequence}->{'PROTEINS'} };
    foreach my $protein (@protein_names) {
      push @{ $proteins{$protein} }, $peptide_sequence;
    }
  }

  # go through the proteins mapping peptides to them
  foreach my $protein_name (keys %proteins) {

    print "$protein_name " . @{ $proteins{$protein_name} } . "\n" if ($verbose);
    my @peptides_in_protein = @{ $proteins{$protein_name} };

    # go through each of the historical wormpep proteins for this gene
    # looking for the most recent one that all the peptides map to
    # and get the version of that protein so we can then get the history CDS if it is not the latest version

    my $mapped_ok = 0;		# flag for mapping peptides to this protein sequence OK
    my ($wormpep_ids_aref, $versions_aref) = &get_previous_wormpep_ids($protein_name, $protein_history_aref);
    my @wormpep_ids = reverse @{$wormpep_ids_aref};
    my @versions = reverse @{$versions_aref};

    my $latest_cds = 0;		# count of the previous revisions of genes that we have to go back through before finding a match
    foreach my $wormpep_id (@wormpep_ids) {

      my $wormpep_seq = &get_wormpep_by_wp_id($wormpep_id);

      # see if the peptides all map to this version of the wormpep protein sequence
      my ($ok, %protein_positions) = &map_peptides_to_protein($wormpep_seq, @peptides_in_protein);

      if ($ok) {
	print "We have found all the peptides in CDS: $protein_name wormpep_id: $wormpep_id version=$latest_cds\n" if ($verbose);

	# flag is true if we can map this to the genome
	my $maps_to_genome_ok = 1;

	# see if we are using the latest version of this protein or not
	# if not, then construct the history ID name of the CDS to get
	if (defined $versions[$latest_cds] && $latest_cds == 0) {
	  $log->write_to("*** $protein_name has been changed to isoforms - it won't have the history CDS for exons, so can't map its peptides to the genome\n");
	  $maps_to_genome_ok = 0;

	} elsif (defined $versions[$latest_cds] || $latest_cds != 0) {
	  $protein_name = $protein_name . ":wp" . ($versions[$latest_cds] - 1);
	  print "$protein_name is using history version $versions[$latest_cds]\n" if ($verbose);
	}

	my ($clone, $cds_start, $cds_end, $exons_start_ref, $exons_end_ref);
	if ($maps_to_genome_ok) {
	  # get the details for mapping to the genome
	  # we can't do this outside of this $wormpep_id loop because we don't know until now whether we must use a history CDS
	  ($clone, $cds_start, $cds_end, $exons_start_ref, $exons_end_ref) = &get_cds_details($protein_name);
	  # if we couldn't get the details of the CDS, then we just have to skip writing the genome mapping data for this protein
	  if (! defined $clone) {
	    $log->write_to("*** Couldn't get the mapping details of the CDS $protein_name\n");
	    $maps_to_genome_ok = 0;
	  } 
	}

        # only write this genome mapping stuff if we could find a CDS for the protein the peptides matched
	if ($maps_to_genome_ok) {
	  # note which clones were used for genome mapping
	  $unique_clones{$clone} = 1;
	}

	# output the details for each of the peptides for this protein
	foreach my $ms_peptide (@peptides_in_protein) {
	  
	  # only write this genome mapping stuff if we could find a CDS for the protein this peptide matched
	  if ($maps_to_genome_ok) {
	    # get the details for mapping to the genome
	    my (@homol_lines) = &get_genome_mapping($ms_peptide, $protein_name, $clone, $cds_start, $cds_end, $exons_start_ref, $exons_end_ref, \%protein_positions); # get part of the MSPeptide_homol line like "clone_start clone_end 1 peptide_length AlignPepDNA"
	    # write the Homol_data for the genome mapping
	    print OUT "\n";
	    print OUT "Homol_data : \"$clone:Mass-spec\"\n";
	    foreach my $homol_line (@homol_lines) {
	      print OUT "MSPeptide_homol \"MSP:$ms_peptide\" mass_spec_genome 100.0 $homol_line\n";
	    }
	  }

	  # count the number of times this peptide maps to a protein in this experiment
	  $peptide_count{$ms_peptide}++;

	  # output the MS_peptide_results hash in the Mass_spec_peptide object
	  print OUT "\n";
	  print OUT "Mass_spec_peptide : \"MSP:$ms_peptide\"\n";
	  print OUT "Mass_spec_experiments \"$experiment_id\" Protein_probability ", $experiment{$experiment_id}{'PROTEIN_PROBABILITY'}{$ms_peptide} ,"\n" if (defined $experiment{$experiment_id}{'PROTEIN_PROBABILITY'}{$ms_peptide});
          print OUT "Mass_spec_experiments \"$experiment_id\" Peptide_probability ", $experiment{$experiment_id}{'PEPTIDE_PROBABILITY'}{$ms_peptide} ,"\n" if (defined $experiment{$experiment_id}{'PEPTIDE_PROBABILITY'}{$ms_peptide});
	  print OUT "Mass_spec_experiments \"$experiment_id\" Protein \"WP:$wormpep_id\"\n"; # the protein that this peptide in this experiment matches

	  # output the protein object for this peptide
	  print OUT "\n";
	  print OUT "Protein : \"MSP:$ms_peptide\"\n";
	  print OUT "Peptide \"MSP:$ms_peptide\"\n";

	  # output the homology of the wormpep protein to the peptide
	  print OUT "\n";
	  print OUT "Protein : \"WP:$wormpep_id\"\n";                                          
	  my $pos = $protein_positions{$ms_peptide}; # position in protein
	  my $len = length($ms_peptide); # length of peptide sequence
	  my $end = $pos+$len-1;
          print OUT "Pep_homol \"MSP:$ms_peptide\" mass-spec 1 $pos $end 1 $len\n"; # note where the peptide maps to a wormpep protein

	  # note that this peptide has been used so that we need to output it when only doing a small debugging output
	  $peptides_used{$ms_peptide} = 1;
	}

	$mapped_ok = 1;
	last;
      }

      # count the version back we have to go
      $latest_cds++;

    }
    if ($latest_cds != 1) {
      print "Protein : $protein_name required revision $latest_cds\n";
    }
    if (!$mapped_ok) {
      $log->write_to("*** Not mapped all the peptides for CDS $protein_name\n");
    }
  }

  # output the experiment data object
  print OUT "\n";
  print OUT "Mass_spec_experiment : \"$experiment_id\"\n";
  print OUT "Species \"Caenorhabditis elegans\"\n";
  print OUT "Life_stage \"mixed stage\"\n";
  print OUT "Sub_cellular_localization \"whole organism\"\n";
  print OUT "Digestion \"Trypsin\"\n";
  print OUT "Instrumentation \"QTOF\"\n";
  print OUT "Multiple_ambiguous_IDs_allowed\n";
  print OUT "//Remark \"(some description of the experimental technique etc.)\"\n";
  print OUT "//Person \"(WBPerson ID)\" \n";
  print OUT "//Author \"(Author if there is no WBPerson ID you can use)\" \n";
  print OUT "//Reference \"(WBPaper ID)\" \n";
  print OUT "//Database \"(A brief description of the sequence database used to search against)\"\n";
  print OUT "//Program \"(Mascot or SEQUEST etc.)\"\n";
  print OUT "//Laboratory \" \n";
  print OUT "//Strain \" \n";
  print OUT "//Genotype \" \n";
  print OUT "//Anatomy_term \"\n";
  print OUT "//Cell_type \" \n";
  print OUT "//Ionisation_source \" \n";
  print OUT "//Minimum_ion_proportion \n";
  print OUT "//False_discovery_rate \n";


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
      print OUT "Peptide_is_natural\n" if ($natural);
      print OUT "Mass_spec_experiments \"$experiment_id\" Matches_database_uniquely\n" if (exists $peptide_count{$ms_peptide} && $peptide_count{$ms_peptide} == 1);
      
      # and output the Peptide sequence object 
      print OUT "\n";	
      print OUT "Peptide : \"MSP:$ms_peptide\"\n";
      print OUT "$ms_peptide\n";
    }
  }
}

# write out the Homol_data lines for the clones and superlinks used
my %clone_sizes = $wormbase->FetchData("clonesize", undef, $database . "/COMMON_DATA");
my $coords      = Coords_converter->invoke($database, 0, $wormbase);

foreach my $clone (keys %unique_clones) {
  my $size;
  if ($clone =~ /^SUPERLINK/) {
    $size = $coords->Superlink_length($clone);
  } else {
    $size = $clone_sizes{$clone};
  }
  print OUT "\n";	
  print OUT "Sequence : $clone\n";
  print OUT "Homol_data $clone:Mass-spec 1 $size\n";
}



close (OUT);

# close the ACE connection
$db->close;

# Close log files and exit
$log->write_to("\n\nStatistics\n");
$log->write_to("----------\n\n");
$log->write_to("Put some statistics here.\n");

$log->mail();
print "Finished.\n" if ($verbose);
exit(0);






##############################################################
#
# Subroutines
#
##############################################################

# example peptide records
#FMRFamide-related peptides or FaRPs
#-LRFamide family
# flp-1  SADPNFLRF . .
# flp-1  AAADPNFLRF . .
# flp-18 EIPGVLRF . .
# flp-18 SEVPGVLRF . .
# flp-18 SYFDEKKSVPGVLRF . .
# flp-18 SVPGVLRF . .
# flp-18 DFDGAMPGVLRF . .
# flp-18 GAMPGVLRF . .
# 
#                 
##-IRFamide family
# flp-5  GAKFIRF . .
# flp-13 APEASPFIRF . .
# flp-13 AMDSPLIRF  . .

# flp-13 ASPSAPLIRF . .
# flp-13 SPSAVPLIRF . .
# flp-13 SAAAPLIRF . .
# flp-13 AADGAPLIRF . .



##########################################
# this gives a hash with key=experiment ID value=list of details
# %experiment = parse_data();

sub parse_data {

  my %experiment;
  my %peptide;
  my $peptide_sequence;
  my $peptide_probability;
  my $protein_probability;
  my $strain_name;
  my $stage_name;
  my $sub_cellular_fraction;
  my $purification_type;
  my $quantification_type;
  my $protein_name;
  my @protein_names;
  my @experiments; # the list of experiments that found the current peptide


  my %cgc_names = $wormbase->FetchData("cds2cgc", undef, $database . "/COMMON_DATA");
  my %cds_cgc_names;
  foreach my $cgc (keys %cgc_names) { # invert the cgc_names hash
    my $cds = $cgc_names{$cgc};
    push @{$cds_cgc_names{$cds}}, $cgc;
  }


  push @experiments, $experiment_id;


  open (DATA, "< $input") || die "Can't open $input\n";
  while (my $line = <DATA>) {
    chomp $line;
    if ($line =~ /^\#/) {next;} 
    if ($line =~ /^\s*$/) {next;}

    my $cgc_name;
    if ($line =~ /\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)/) {
      $cgc_name = $1;
      $peptide_sequence = $2;
      $protein_probability = $3;
      if ($protein_probability eq '.') {$protein_probability = undef}
      $peptide_probability = $4;
      if ($peptide_probability eq '.') {$peptide_probability = undef}
      @protein_names = ();
    }

    if ($cgc_name =~ /\S+\-\d+/) { # standard CGC-style name e.g. daf-1
      if (!exists $cds_cgc_names{$cgc_name}) {
	print "line $line has a non-existent cds_cgc name\n";
	next;
      }
      # start getting the protein data
      push @protein_names, @{$cds_cgc_names{$cgc_name}}; # get the list of CDS names for this CGC name
  
    } elsif ($cgc_name =~ /\S+\.\d+/) {	# standard Sequence style CDS name e.g. C25A7.2
      push @protein_names, $cgc_name; # use CDS name as-is

    } else {
      die "gene name '$cgc_name' not recognised\n";
    }

    # construct the peptide data
    $peptide{$peptide_sequence}->{'PROTEINS'} = [@protein_names];
    #print "$peptide_sequence protein = @protein_names\n";

    # save the peptide data in all its experiments
    foreach my $experiment_id (@experiments) {
      #print "experiment $experiment_id = $peptide_sequence\n";
      push @{ $experiment{$experiment_id}->{'PEPTIDES'} }, $peptide_sequence;
      $experiment{$experiment_id}{'PROTEIN_PROBABILITY'}{$peptide_sequence} = $protein_probability;
      $experiment{$experiment_id}{'PEPTIDE_PROBABILITY'}{$peptide_sequence} = $peptide_probability;
    }
  }
  close (DATA);

  return (\%experiment, \%peptide);
}


##########################################
# take a protein and a set of peptides and map the peptides onto the protein
# returns a hash key=peptide sequence, value=position in protein

sub map_peptides_to_protein {
  my ($protein, @peptides) = @_;
  my %position;
  my $ok = 1;			# true if all the peptides map to this protein

  foreach my $peptide (@peptides) {
    my $pos = index($protein, $peptide);
    $position{$peptide} = $pos+1; # convert from offset=0 to offset=1;
    if ($pos == -1) {$ok = 0;}	# set the flag if a peptide is not found in the protein
  }

  return ($ok, %position);
}

##########################################
# get the set of CDS, history versions and resulting protein IDs
# my %protein_history = &get_protein_history;

sub get_protein_history {
  my @protein_history;

  # get database version
  my $version = `grep "NAME WS" $database/wspec/database.wrm`;
  chomp($version);
  $version =~ s/.*WS//;
  my $data_file = "/nfs/disk100/wormpub/BUILD/WORMPEP/wormpep${version}/wormpep.history$version";

  open (HIST, "< $data_file") || die "Can't open $data_file\n";
  while (my $line = <HIST>) {
    chomp $line;
    my @f = split /\s+/, $line;
    # $cds_id, $wormpep_id, $version1, $version2 <- ($version2 is undef if current version)
    push @protein_history, [@f];
  }
  close(HIST);

  return \@protein_history;
}


##########################################
# get all previous wormpep protein IDs for a given CDS
# returned in the order in which they occur in the data file
# which seems to be in historical order (oldest first)
# @worpep_ids = get_previous_wormpep_ids($gene_name);

sub get_previous_wormpep_ids {

  my ($gene_name, $protein_history_aref) = @_;

  my @protein_history = @{$protein_history_aref};

  my @wormpep_ids;
  my @versions;

  foreach my $line (@protein_history) {
    my ($cds_id, $wormpep_id, $version1, $version2) = @{$line};
    if ($cds_id eq $gene_name) {
      print "WP ID $wormpep_id\n";
      push @wormpep_ids, $wormpep_id;
      if (defined $version2) {	# if this is an old version, save the Build release number of the change

# I think this should be: 
#	push @versions, $version2;
# not
#        push @versions, $version1;
# because version2 is the CDS version number at which the CDS was
# changed and the protein is defined by CDS:wp(version2-1)

	push @versions, $version2;
      } else {			# otherwise indicate that it is the current version by an 'undef'
	push @versions, undef;
      }
    }
  }

  return \@wormpep_ids, \@versions;
}

##########################################
# get the wormpep protein from the wormpep fasta file by wormpep ID
# $seq = get_wormpep_by_wp_id($wormpep_id);

sub get_wormpep_by_wp_id {

  my ($id) = @_;
  my $data_file = "/nfs/disk100/wormpub/WORMPEP/wp.fasta_current";
  my $title="";
  my $seq="";

  open (SEQ, "<$data_file") || die "Can't open $data_file\n";
  while (my $line = <SEQ>) {
    chomp $line;
    if ($line =~ />(\S+)/) {
      if ($title eq $id) { # NB this is the previous title, not the one on this line
	return $seq;
      }
      $seq = "";
      $title = $1;
    } else {
      $seq .= $line;
    }
  }
  close (SEQ);

  if ($title eq $id) { # NB this is the last title in the file
    return $seq;
  }
}

##########################################
# my ($clone, $cds_start, $exons_start_ref, $exons_end_ref) = &get_cds_details($protein_name);
# get the clone, clone start position and exons list from a CDS
#

sub get_cds_details {
  my ($protein_name) = @_;

  print "Ace fetch->" if ($verbose);
  print "$protein_name\n" if ($verbose);
  my $cds_obj = $db->fetch("CDS" => $protein_name);
  # debug
  if (! defined $cds_obj) {
    print "Can't fetch CDS object for $protein_name\n";
    return (undef, undef, undef, undef, undef);
  }

  my $clone = $cds_obj->Sequence;
  if (! defined $clone) {
    print "Can't fetch clone for $protein_name\n";
    return (undef, undef, undef, undef, undef);
  }

  # get the named object
  # then get the data following a tag in that object
  # and the NEXT data following the tag
  my @exons_start = $cds_obj->Source_exons;
  my @exons_end = $cds_obj->Source_exons(2);
  $cds_obj->DESTROY();
  if (! @exons_start || ! @exons_end) {
    print "Can't fetch exons for $protein_name\n";
    return (undef, undef, undef, undef, undef);
  }

  # get the start position of the CDS in the clone's SMap
  print "Ace fetch->clone = $clone\n" if ($verbose);
  my $clone_obj = $db->fetch("Sequence" => $clone);
  if (! defined $clone_obj) {
    print "Can't fetch clone object for $clone\n";
    return (undef, undef, undef, undef, undef);
  }

  # the original routine is no longer used because it didn't
  # deal with CDS names that had missing start/end positions which
  # messed up the synchrony between the names and positions from that
  # point onwards.
  #
  # You have to creep up on the data you want by using the $obj->at()
  # structure, in which you can specify exactly the path of tags and
  # data through the object that you wish to traverse and then you can
  # get the list of resulting positions under the sub-path that you
  # want. This means that you must break up the query into a loop
  # looking at every score value of every wormbase homology and then
  # you can get the positions under that score.

  my $foundit = 0;
  my ($cds, $cds_start, $cds_end);
  foreach $clone_obj ($clone_obj->at('SMap.S_child.CDS_child')) {
    ($cds, $cds_start, $cds_end) = $clone_obj->row;
    if ($cds eq $protein_name) {
      $foundit = 1;
      last;
    }
  }

  if (! $foundit || ! defined $cds_start || ! defined $cds_end) {
    print "Can't fetch start/end positions for $protein_name in $clone\n";
    return (undef, undef, undef, undef, undef);
  }

  print "**** Found $protein_name in $clone:$cds_start..$cds_end\n";

  # if the CDS is on the reverse sense, then $cds_end > $cds_start
  return ($clone, $cds_start, $cds_end, \@exons_start, \@exons_end);

}


##########################################
# my @homol_lines = &get_genome_mapping($ms_peptide, $protein_name, $clone, $cds_start, $exons_start_ref, $exons_end_ref, \%protein_positions); 
# get part of the MSPeptide_homol line like "clone_start clone_end 1 peptide_length AlignPepDNA"
#
# take the information describing the position of the CDS's exons in
# the clone and the peptides position in the CDS and construct
# required part of the MSPeptide_homol data line for output.
#
# Check that the peptide maps to the genome correctly by translating the mapped part

sub get_genome_mapping {
  my ($ms_peptide, $protein_name, $clone, $cds_start, $cds_end, $exons_start_ref, $exons_end_ref, $protein_positions_href) = @_;
  my %protein_positions = %{$protein_positions_href};

  # the positions of the exons on the clone relative to the start of the CDS
  my @exons_start = @{$exons_start_ref};
  my @exons_end = @{$exons_end_ref};

  my @output;			# returned list of output lines

  # the position in the CDS of the start and end of exons
  my @exons_internal_start = ();
  my @exons_internal_end = ();

  # get the peptide length
  my $peptide_length = length($ms_peptide);

  print "ms_peptide = $ms_peptide protein_name = $protein_name\n" if ($verbose);

  # get the start and end position of the peptide in the protein
  my $peptide_cds_start = $protein_positions{$ms_peptide}; 
  print "position of peptide in protein = $peptide_cds_start\n" if ($verbose);
  my $peptide_cds_end = $peptide_cds_start + $peptide_length - 1;
  # convert protein positions to CDS coding positions
  $peptide_cds_start = ($peptide_cds_start * 3) - 2; # start at the beginning of the codon, not the end :-)
  $peptide_cds_end *= 3;
  print "peptide_cds_start = $peptide_cds_start\n" if ($verbose);

  # find the positions in the cDNA of the start and end of the exons (i.e. the splice sites on the cDNA)
  my $exon_count = 0;
  my $prev_end = 0;
  print "number of exons = $#exons_start and $#exons_end\n" if ($verbose);
  foreach my $exon_start (@exons_start) {
    print "exon_start $exon_start exon_end $exons_end[$exon_count]\n" if ($verbose);
    my $exon_length = $exons_end[$exon_count] - $exon_start + 1;
    print "exon_length $exon_length\n" if ($verbose);
    my $start = $prev_end + 1;
    my $end = $start + $exon_length - 1;
    print "start $start end $end\n" if ($verbose);
    push @exons_internal_start, $start;
    push @exons_internal_end, $end;
    $prev_end = $end;
    $exon_count++;
  }
  print "number of internal coords exons = $#exons_internal_start+1 and $#exons_internal_end+1\n" if ($verbose);


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

  foreach $exon_count (0..scalar @exons_start) {

    print "exon_count = $exon_count\n" if ($verbose);
    print "peptide_cds_start = $peptide_cds_start\n" if ($verbose);
    print "peptide_cds_end = $peptide_cds_end\n" if ($verbose);
    print "exon CDS boundaries = $exons_internal_start[$exon_count] $exons_internal_end[$exon_count]\n" if ($verbose);

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
#      $peptide_start = int (($exons_internal_start[$exon_count] - $peptide_cds_start) / 3) + 1;
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
#      $peptide_start = int (($exons_internal_start[$exon_count] - $peptide_cds_start) / 3) + 1;
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

1) matching CDS sequence name or gene WGN (CGC) name
2) mass-spec peptide sequnce
3) protein probability (if known, otherwise leave as blank or '.')
4) peptide probability (if known, otherwise leave as blank or '.')

The output file may have to be edited to add further details of the experiment
after the line "Mass_spec_experiment : "

script_template.pl MANDATORY arguments:

=over 4

=item -input file to read the mass-spec data from

=back

=over 4

=item -output file to write ace results to.

=back

=over 4

=item -experiment_id name of the experiment, usually formed from the initials of the author, plus a number, e.g. 'SH_1'

=back

script_template.pl  OPTIONAL arguments:

=over 4

=item -database database used to map mass spec peptides to, default is currentDB

=back

=over 4

=item -natural specifies that the peptides in this data set are naturally occuring.

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
