#!/usr/bin/env perl
#
# RNASeq_make_analysis_objects.pl
# Small script to take the pain out of making Analysis objects
# Give it the species and the name of the Study and it will help make the analysis objects

# by Gary Williams                        
#
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2015-03-20 14:37:25 $      

# Run as:
# perl ~gw3/bin/RNASeq_make_analysis_objects.pl -species elegans -study SRP021083 -paper WBPaper00044426 -out SRP021083.ace
# Or:
# perl ~gw3/bin/RNASeq_make_analysis_objects.pl -species elegans -suggest

use strict;
use lib $ENV{'CVS_DIR'};
use Carp;
use Modules::RNASeq;

use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

my ($help, $debug, $test, $verbose, $store, $wormbase, $species, $study, $outfile, $paper, $output, $redo, $suggest);
GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
            "test"       => \$test,
            "verbose"    => \$verbose,
            "store:s"    => \$store,
            "species:s"  => \$species, # the default is elegans
	    "study:s"    => \$study,   # the study to make experiment analysis objects for
	    "paper:s"    => \$paper,   # the WBPaper ID to use
	    "output:s"   => \$output,  # the output ace file ready for you to read into geneace
	    "redo"       => \$redo,    # do even though an analysis already exists
	    "suggest"    => \$suggest, # suggest potential Studies that could now be curated (are the correct types and have a paper, etc.)
);

$debug = "gw3";

$species = 'elegans' unless $species;

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
                             -organism => $species,
                             );
}

# establish log file.
my $log = Log_files->make_build_log($wormbase);

my $database = $wormbase->wormpub . "/DATABASES/current_DB";

my $RNASeq = RNASeq->new($wormbase, $log, 0, 0);

if (defined $suggest) {
  suggest(); 
  $log->mail();
  exit(0);
}

if (!defined $output) {$log->log_and_die("-output not defined\n")}
if (!defined $study) {$log->log_and_die("-study not defined\n")}
if (!defined $paper) {$log->log_and_die("-paper not defined\n")}
if (!defined $redo) {$redo=0}

######################################
# variables and command-line options # 
######################################

my %strain_session;
my %sex_session;
my %tissue_session;
my %tissue_ontology;
my %life_stage_session;
my %life_stage_ontology;

my %data;
my $data;
my $name;
my $life_stage;
my $life_stage_name;
my $strain;
my $sex;
my $tissue;
my $tissue_name;
my $treatment;
my $genotype;
my $temperature;

my $study_title;
my $study_alias;
my $experiment_title;
my $library_name;
my $run_alias;

my $written_Study_analysis = 0;

open(OUT, ">$output") || $log->log_and_die("Can't open $output\n");

my %anatomy = get_anatomy();

# find all experiments in the Study
my %hash;
my $hashref = \%hash;
my $studies_ini = $RNASeq->read_all_studies_config();
my @studies = $studies_ini->Sections;
foreach my $study_accession (@studies) {
  if ($study eq $study_accession) {
    my $experiments_ini = $RNASeq->read_experiments_from_study_config($study_accession);
    $RNASeq->convert_ini_to_hashref($experiments_ini, $hashref);

    # Study short description displayed on website box-plots
    my $study_display = get_study_display();

    foreach my $expt (keys %{$hashref}) {
      if ($hashref->{$expt}{library_source} ne 'TRANSCRIPTOMIC') {print "\n$expt is not TRANSCRIPTOMIC - skipping\n"; next}
      if ($hashref->{$expt}{library_strategy} ne 'RNA-Seq' && 
	  $hashref->{$expt}{library_strategy} ne 'FL-cDNA' && 
	  $hashref->{$expt}{library_strategy} ne 'EST') {print "\n$expt is not RNA-Seq or FL-cDNA or EST - skipping\n"; next}
      if ($hashref->{$expt}{library_selection} eq'size fractionation') {print "\n$expt is size selected - skipping\n"; next}
      if (exists $hashref->{$expt}{ignore}) {print "\n$expt is marked as 'ignore' - skipping\n"; next}
      if (exists $hashref->{$expt}{analysis} && !$redo) {print "\n$expt has an Analysis noted already in '${study_accession}.ini' (",$hashref->{$expt}{analysis},") - skipping\n"; next}

      # get the study_title, study_alias, experiment_title, library_name, run_alias and do things with them
      $study_title      = $hashref->{$expt}{study_title};
      $study_alias      = $hashref->{$expt}{study_alias};
      $experiment_title = $hashref->{$expt}{experiment_title};
      $library_name     = $hashref->{$expt}{library_name};
      $run_alias        = $hashref->{$expt}{run_alias};
      my ($dummy, $ENA_dev_stage, $ENA_sex, $ENA_strain, $ENA_temperature, $ENA_tissue, $ENA_description) = query_ENA($expt);

      print "\nExperiment: $expt\n";
      print "Study title: $study_title\n" if ($study_title ne '');
      print "Study alias: $study_alias\n" if ($study_alias ne '');
      print "Experiment title: $experiment_title\n" if ($experiment_title ne '');
      print "Library name: $library_name\n" if ($library_name ne '');
      print "Run alias: $run_alias\n" if ($run_alias ne '');
      print "ENA description: $ENA_description\n" if ($ENA_description);
      print "\n";


      # look for temperature
      # default is nothing
      $temperature = get_temperature($ENA_temperature);

      # look for genotype
      # default is nothing
      $genotype = get_genotype();

      # look for treatment
      # default is nothing
      $treatment = get_treatment();

      # look for strain
      # default is initially N2 in elegans
      if (!defined $genotype || $genotype eq "") {
	$strain = get_strain($ENA_strain);
      }

      # look for sex
      # one probable synonym is 'him-8:Male'
      # synonym / lookup will be different in different species
      $sex = get_sex($ENA_sex);

      # look for tissue
      # use synonym / lookup hash with recently found synonyms being checked first in the title they were last found in
      # default is initially whole organism
      ($tissue, $tissue_name) = get_tissue($ENA_tissue, %anatomy);
      
      # look for life stage
      # use synonym / lookup hash with recently found synonyms being checked first in the title they were last found in
      # synonym / lookup will be different in different species
      ($life_stage, $life_stage_name) = get_life_stage($ENA_dev_stage);

      # backslash double quotes in $study_title, "$experiment_title
      $study_title =~ s/"/\"/;
      $experiment_title =~ s/"/\"/;

      # write the Study Analysis
      my $full_species = $wormbase->full_name;
      my $study_name = "RNASeq_Study.$study";

      if (!$written_Study_analysis) {
	print OUT "\n";
	print OUT "Analysis : \"RNASeq_Study.$study\"\n";
	print OUT "Database SRA Study \"$study\"\n";
	print OUT "Title \"$study_display\"\n";
	print OUT "Description \"$study_title\"\n" if ($study_title ne '');
	print OUT "Reference \"$paper\"\n";
	print OUT "Species_in_analysis \"$full_species\"\n";
	print OUT "Project \"RNASeq.$species\"\n";
	$written_Study_analysis = 1;
      }

      my $name = "RNASeq.$species.$strain.${life_stage}.$sex.$tissue.$study.$expt";

      print OUT "\n";
      print OUT "Analysis : \"$name\"\n";
      print OUT "Database SRA SRA \"$expt\"\n";
      print OUT "Database SRA Study \"$study\"\n";
      print OUT "Title \"$study_title\"\n";
      print OUT "Description \"$experiment_title\"\n";
      print OUT "Reference \"$paper\"\n";
      print OUT "Sample \"$name\"\n";
      print OUT "Project \"$study_name\"\n";
      print OUT "Species_in_analysis \"$full_species\"\n";
      print OUT "\n";
      print OUT "Condition : \"$name\"\n";
      print OUT "Life_stage \"$life_stage\"\n";
      print OUT "Strain \"$strain\"\n" if (!defined $genotype || $genotype eq "");
      print OUT "Sex \"$sex\"\n";
      print OUT "Reference \"$paper\"\n";
      print OUT "Tissue \"$tissue\"\n";
      print OUT "Genotype \"$genotype\"\n" if (defined $genotype && $genotype ne "");
      print OUT "Treatment \"$treatment\"\n" if (defined $treatment && $treatment ne "");
      print OUT "Temperature \"$temperature\"\n" if (defined $temperature && $temperature ne "");
      print OUT "Species \"$full_species\"\n";

      # store the analysis name etc. in the INI file
      $experiments_ini->newval($expt, 'analysis', $name);
      $experiments_ini->RewriteConfig;
      if (!defined $genotype || $genotype eq "") {
	$experiments_ini->newval($expt, 'WBstrain', $strain);
	$experiments_ini->RewriteConfig;
      }
      $experiments_ini->newval($expt, 'WBls_code', $life_stage);
      $experiments_ini->RewriteConfig;
      $experiments_ini->newval($expt, 'WBls_name', $life_stage_name);
      $experiments_ini->RewriteConfig;
      $experiments_ini->newval($expt, 'WBsex', $sex);
      $experiments_ini->RewriteConfig;
      $experiments_ini->newval($expt, 'WBanatomy_code', $tissue);
      $experiments_ini->RewriteConfig;
      $experiments_ini->newval($expt, 'WBanatomy_term', $tissue_name);
      $experiments_ini->RewriteConfig;
      if (defined $genotype && $genotype ne "") {
	$experiments_ini->newval($expt, 'WBgenotype', $genotype);
	$experiments_ini->RewriteConfig;
      }
      if (defined $treatment && $treatment ne "") {
	$experiments_ini->newval($expt, 'WBtreatment', $treatment);
	$experiments_ini->RewriteConfig;
      }

    }
  }
}

close(OUT);

$log->mail();
print "Finished.\n";
exit(0);

#############################################################################
# for temperature, we just accept what is given with default always being blank
sub get_temperature {
  my ($ENA_temperature) = @_;


  my $candidate_temp = '';
  my $why = '';

  # assume we are using the same sex as in the last experiment
  if (defined $ENA_temperature) {
    $candidate_temp = $ENA_temperature;
    $why = '(from ENA)';
  }

  print "Temp [$candidate_temp] $why > ";
  my $error;
  my $input;
  do {
    $input =  <STDIN>;
    chomp ($input);
    if ($input eq '') {
      $input = $candidate_temp;
    }
    if ($input ne '' && $input !~ /^[\d\.]+$/) {
      print "Temperature must be numeric\n";
      $error = 1;
    }

  } while ($error);
  return $input;
}


#############################################################################
# for treatment, we just accept what is given with default always being blank
sub get_treatment {

  print "Treatment > ";
  my $input =  <STDIN>;
  chomp ($input);

  return $input;
}


#############################################################################
# for genotype, we just accept what is given with default always being blank
sub get_genotype {

  print "Genotype > ";
  my $input =  <STDIN>;
  chomp ($input);

  return $input;
}


#############################################################################
# We don't use a stored hash of synonyms and codes etc for strains
# we just look for likely words and remember what was used before in this session
# but we have a default of N2 for elegans

sub get_strain {
  my ($ENA_strain) = @_;

  my $strain; 
  my $candidate_strain = undef;
  my $line;
  my $why = '';

  my %default = (
		 'elegans'   => 'N2',
		 'brenneri'  => 'PB2801',
		 'briggsae'  => 'AF16',
		 'brugia'    => 'FR3',
		 'sratti'    => 'ED321',
		 'japonica'  => 'DF5081',
		 'ovolvulus' => 'O_volvulus_Cameroon_isolate',
		 'remanei'   => 'SB146',
		 'pristionchus'   => 'PS312',
		 'tmuris'    => 'Edinburgh',
		);
  my $default = $default{$species};

  # see if we can find a strain already seen in this study
 NAMES:
  foreach my $test_name (@{$strain_session{names_found}}) {
    foreach my $line ($run_alias, $experiment_title, $study_alias, $study_title) {
      if ($line =~ /$test_name/) {
	$candidate_strain = $test_name;
	$why = 'in description';
	last NAMES;
      }
    }    
  }

  # find an uppercase word preceeding or following 'strain'
  if (!defined $candidate_strain) {
    foreach my $line ($run_alias, $experiment_title, $study_alias, $study_title) {
      my ($before) = ($line =~ /(\w+)\s[Ss]train/);
      my $tmp = $before;
      if ($tmp eq uc $before) {
	$candidate_strain = $before;
	$why = "before 'strain'";
	last;
      }
      my ($after) = ($line =~ /[Ss]train\s(\w+)/);
      $tmp = $after;
      if ($tmp eq uc $after) {
	$candidate_strain = $after;
	$why = "after 'strain'";
	last;
      }
    }
  }

  if (!defined  $candidate_strain) {
    if (defined $ENA_strain) {
      $candidate_strain = $ENA_strain;
      $why = 'from ENA';
    }
  }

  # assume we are using the same strain as in the last experiment
  if (!defined $candidate_strain && defined $strain_session{previous}) {
    $candidate_strain = $strain_session{previous};
    $why = 'repeat of last expt';
  }

  # else if this is the first time in, use the default
  if (!defined $candidate_strain) {
    $candidate_strain = $default;
    $why = 'default';
  }

  # confirm and update status
  print "Strain [$candidate_strain] ($why) > ";
  my $input =  <STDIN>;
  chomp ($input);
  if ($input eq '') {
    $strain = $candidate_strain;
  } else {
    $strain = $input;
    push @{$strain_session{names_found}}, $strain;
  }
  $strain_session{previous} = $strain;

  return $strain;
}

#############################################################################
# We don't use a stored hash of synonyms and codes etc for sex
# we just look for likely words and remember what was used before in this session
# Females are Hermaphrodites in elegans

sub get_sex {
  my ($ENA_sex) = @_;

  my $sex; 
  my $candidate_sex = undef;
  my $line;
  my $why = '';

  my %default = (
		 'elegans'   => 'Hermaphrodite',
		 'brenneri'  => 'Unknown',
		 'briggsae'  => 'Hermaphrodite',
		 'sratti'    => 'Unknown',
		 'brugia'    => 'Unknown',
		 'japonica'  => 'Unknown',
		 'ovolvulus' => 'Unknown',
		 'remanei'   => 'Unknown',
		 'tmuris'    => 'Unknown',
		);
  my $default = $default{$species};

  # see if we can find a sex already seen in this study
 NAMES:
  foreach my $test_name (@{$sex_session{names_found}}) {
    foreach my $line ($run_alias, $experiment_title, $study_alias, $study_title) {
      if ($line =~ /$test_name/) {
	$candidate_sex = $test_name;
	$why = 'in description';
	last NAMES;
      }
    }    
  }

  if (!defined  $candidate_sex) {
    if (defined $ENA_sex) {
      $candidate_sex = $ENA_sex;
      $why = 'from ENA';
    }
  }

  # look for sex in the usual places
  if (!defined $candidate_sex) {
    foreach my $line ($run_alias, $experiment_title, $study_alias, $study_title) {
      if ($line =~ /female/i) {
	$candidate_sex = 'Female';
	if ($species eq 'elegans') {$candidate_sex = 'Hermaphrodite';} # probably true in most cases
	$why = 'in description';
      } elsif ($line =~ /hermaphrodite/i) {
	$candidate_sex = 'Hermaphrodite';
	$why = 'in description';
      } elsif ($line =~ /male/i) {
	$candidate_sex = 'Male';
	$why = 'in description';

      }
    }
  }


  # in elegans him-8 indicates a Male
  if (!defined $candidate_sex) {
    if ($species eq 'elegans') {
      foreach my $line ($run_alias, $experiment_title, $study_alias, $study_title) {
	if ($line =~ /him-8/) {
	  $candidate_sex = 'Male';
	  $why = 'him-8';
	}
      }
    }
  }

  # assume we are using the same sex as in the last experiment
  if (!defined $candidate_sex && defined $sex_session{previous}) {
    $candidate_sex = $sex_session{previous};
    $why = 'repeat of last expt';
  }

  # else if this is the first time in, use the default
  if (!defined $candidate_sex) {
    $candidate_sex = $default;
    $why = 'default';
  }

  if ($candidate_sex eq 'mixed sex') {
    $candidate_sex = 'Unknown';
  }

  # confirm and update status
  print "Sex [$candidate_sex] ($why) > ";
  my $input =  <STDIN>;
  chomp ($input);
  if ($input eq '') {
    $sex = $candidate_sex;
  } else {
    $sex = $input;
    push @{$sex_session{names_found}}, $sex;
  }
  $sex_session{previous} = $sex;

  return $sex;
}

#############################################################################
# We use a stored hash of synonyms and codes etc for tissue

sub get_tissue {
  my ($ENA_tissue, %tissue_ontology) = @_;

  my $tissue = undef; 
  my $tissue_name = undef; 
  my $candidate_tissue = undef;
  my $line;
  my $why = '';

  my $default = 'organism'; # WBbt:0007833'
 

  # see if we can find a tissue already seen in this study
 NAMES:
  foreach my $test_name (@{$tissue_session{names_found}}) {
    foreach my $line ($run_alias, $experiment_title, $study_alias, $study_title) {
      if ($line =~ /\b$test_name\b/) {
	$candidate_tissue = $test_name;
	$why = 'in description';
	last NAMES;
      }
    }
  }
  
  # search for tissue names from the ontology in the descriptions
  # search for the longest names first to get the most exact term matching first
 ONTOLOGIES:
  foreach my $tissue_term (sort {length $b <=> length $a} (keys %tissue_ontology)) {
    foreach my $line ($run_alias, $experiment_title, $study_alias, $study_title) {
      if ($line =~ /\b$tissue_term\b/) {
	$candidate_tissue = $tissue_term;
	$why = 'in description';
	last ONTOLOGIES;
      }
    }    
  }

  if (!defined  $candidate_tissue) {
    if (defined $ENA_tissue) {
      $candidate_tissue = $ENA_tissue;
      $why = 'from ENA';
    }
  }

  # look for sex in the usual places  # assume we are using the same tissue as in the last experiment
  if (!defined $candidate_tissue && defined $tissue_session{previous}) {
    $candidate_tissue = $tissue_session{previous};
    $why = 'repeat of last expt';
  }

  # else if this is the first time in, use the default
  if (!defined $candidate_tissue || $candidate_tissue eq '') {
    $candidate_tissue = $default;
    $why = 'default';
  }

  # confirm and update status
  do {
    print "Tissue [$candidate_tissue] ($why) or '?substring' > ";
    my $input =  <STDIN>;
    chomp ($input);
    if ($input eq '') {
      $tissue = $tissue_ontology{$candidate_tissue};
      $tissue_name = $candidate_tissue;
    } elsif ($input =~ /\?/) {
      print "Matches found:\n";
      $input =~ s/\?//g;
      $input = lc $input;
      if ($input ne '') {
	foreach my $term (keys %tissue_ontology) {
	  if (index(lc $term, $input) != -1) {
	    print "\t$term\n";
	  }
	}
      }
    } else {
      $tissue = $tissue_ontology{$input};
      if (!defined $tissue) {
	print "WARNING: $input is not known in the ontology.\n"
      } else {
	push @{$tissue_session{names_found}}, $input;
      $tissue_name = $input;
      }
    }
  } until (defined $tissue);
  $tissue_session{previous} = $tissue_name;

  return ($tissue, $tissue_name);
}

#############################################################################
# We use a stored hash of synonyms and codes etc for life_stage

sub get_life_stage {
  my ($ENA_dev_stage) = @_;

  my $life_stage; 
  my $life_stage_name; 
  my $candidate_life_stage = undef;
  my $line;
  my $why = '';
  
  my $default = 'all stages'; # WBls:0000002
  
  my %life_stage_ontology;
  

  if ($species ne 'brugia' && $species ne 'sratti' && $species ne 'tmuris') {
    
    %life_stage_ontology = (
			    'all stages' => 'WBls:0000002',
			    'embryo' => 'WBls:0000003',
			    'proliferating embryo' => 'WBls:0000004',
			    'blastula embryo' => 'WBls:0000005',

			    'oocyte' => 'WBls:0000669',
			    'unfertilized egg' => 'WBls:0000669',

			    '1-cell embryo' => 'WBls:0000006',
			    'zygote' => 'WBls:0000006',
			    'fertilized egg' => 'WBls:0000006',

			    '2-cell embryo' => 'WBls:0000007',
			  '4-cell embryo' => 'WBls:0000008',
			  '28-cell embryo' => 'WBls:0000009',
			  'gastrulating embryo' => 'WBls:0000010',
			  '51-cell embryo' => 'WBls:0000011',
			  '88-cell embryo' => 'WBls:0000012',
			  'enclosing embryo' => 'WBls:0000013',
			  'late cleavage stage embryo' => 'WBls:0000014',
			  'elongating embryo' => 'WBls:0000015',
			  'bean embryo' => 'WBls:0000016',
			  'comma embryo' => 'WBls:0000017',
			  '1.5-fold embryo' => 'WBls:0000018',
			  '2-fold embryo' => 'WBls:0000019',
			  '3-fold embryo' => 'WBls:0000020',
			  'fully-elongated embryo' => 'WBls:0000021',
			  'postembryonic' => 'WBls:0000022',
			  'larva' => 'WBls:0000023',
			  'L1 larva' => 'WBls:0000024',
			  'L1-L2 lethargus' => 'WBls:0000025',
			  'L1-L2 molt' => 'WBls:0000026',
			  'L2 larva' => 'WBls:0000027',
			  'L2-L3 lethargus' => 'WBls:0000028',
			  'L2-L3 molt' => 'WBls:0000029',
			  'L2d-dauer lethargus' => 'WBls:0000030',
			  'L2d-dauer molt' => 'WBls:0000031',
			  'dauer larva' => 'WBls:0000032',
			  'postdauer-L4 lethargus' => 'WBls:0000033',
			  'postdauer-L4 molt' => 'WBls:0000034',
			  'L3 larva' => 'WBls:0000035',
			  'L3-L4 lethargus' => 'WBls:0000036',
			  'L3-L4 molt' => 'WBls:0000037',
			  'L4 larva' => 'WBls:0000038',
			  'L4-adult lethargus' => 'WBls:0000039',
			  'L4-adult molt' => 'WBls:0000040',
			  'adult' => 'WBls:0000041',
			  'L1-L2 ecdysis' => 'WBls:0000042',
			  'L1-L2d molt' => 'WBls:0000043',
			  'L1-L2d lethargus' => 'WBls:0000044',
			  'L1-L2d ecdysis' => 'WBls:0000045',
			  'L2d larva' => 'WBls:0000046',
			  'L2-L3 ecdysis' => 'WBls:0000047',
			  'L2d-dauer ecdysis' => 'WBls:0000048',
			  'L3-L4 ecdysis' => 'WBls:0000049',
			  'L4-adult ecdysis' => 'WBls:0000050',
			  'postdauer-L4 ecdysis' => 'WBls:0000051',
			  'post dauer stage' => 'WBls:0000052',
			  'L2d-L3 molt' => 'WBls:0000053',
			  'L2d-L3 lethargus' => 'WBls:0000054',
			  'L2d-L3 ecdysis' => 'WBls:0000055',
			  'adult male' => 'WBls:0000056',
			  'adult hermaphrodite' => 'WBls:0000057',
			  'pre-reproductive stage adult hermaphrodite' => 'WBls:0000058',
			  'reproductive stage adult hermaphrodite' => 'WBls:0000060',
			  'oocyte-laying stage adult hermaphrodite' => 'WBls:0000061',
			  'post-reproductive stage adult hermaphrodite' => 'WBls:0000062',
			  'newly molted young adult hermaphrodite' => 'WBls:0000063',
			  '1-day post-L4 adult hermaphrodite' => 'WBls:0000064',
			  '2-days post-L4 adult hermaphrodite' => 'WBls:0000065',
			  '3-days post-L4 adult hermaphrodite' => 'WBls:0000066',
			  '4-days post-L4 adult hermaphrodite' => 'WBls:0000067',
			  '5-days post-L4 adult hermaphrodite' => 'WBls:0000068',
			  '6-days post-L4 adult hermaphrodite' => 'WBls:0000670',
			  '7-days post-L4 adult hermaphrodite' => 'WBls:0000671',
			  '8-days post-L4 adult hermaphrodite' => 'WBls:0000672',
			  '9-days post-L4 adult hermaphrodite' => 'WBls:0000673',
			  '10-days post-L4 adult hermaphrodite' => 'WBls:0000674',
			  '15-days post-L4 adult hermaphrodite' => 'WBls:0000675',
			  '20-days post-L4 adult hermaphrodite' => 'WBls:0000676',

			  '4-7 days post-L4 adult hermaphrodite' => 'WBls:0000069',
			  '7-10 days post-L4 adult hermaphrodite' => 'WBls:0000070',

			  '8-cell embryo' => 'WBls:0000071',
			  '12-cell embryo' => 'WBls:0000072',
			  'L4 larva male' => 'WBls:0000073',
			  '11-15 days post-L4 adult hermaphrodite' => 'WBls:0000074',
			  '14-cell embryo' => 'WBls:0000084',
			  '24-cell embryo' => 'WBls:0000085',
			  '44-cell embryo' => 'WBls:0000086',
			  '68-cell embryo' => 'WBls:0000087',
			  '86-cell embryo' => 'WBls:0000088',
			  '190-cells embryo' => 'WBls:0000089',
			  '96-cell embryo' => 'WBls:0000090',
			  '16-18 days post-L4 adult hermaphrodite' => 'WBls:0000111',

			  '1 min post first-cleavage' => 'WBls:0000112',
			  '2 min post first-cleavage' => 'WBls:0000113',
			  '3 min post first-cleavage' => 'WBls:0000114',
			  '4 min post first-cleavage' => 'WBls:0000115',
			  '5 min post first-cleavage' => 'WBls:0000116',
			  '6 min post first-cleavage' => 'WBls:0000117',
			  '7 min post first-cleavage' => 'WBls:0000118',
			  '8 min post first-cleavage' => 'WBls:0000119',
			  '9 min post first-cleavage' => 'WBls:0000120',
			  '10 min post first-cleavage' => 'WBls:0000121',
			  '11 min post first-cleavage' => 'WBls:0000122',
			  '12 min post first-cleavage' => 'WBls:0000123',
			  '13 min post first-cleavage' => 'WBls:0000124',
			  '14 min post first-cleavage' => 'WBls:0000125',
			  '15 min post first-cleavage' => 'WBls:0000126',
			  '16 min post first-cleavage' => 'WBls:0000127',
			  '17 min post first-cleavage' => 'WBls:0000128',
			  '18 min post first-cleavage' => 'WBls:0000129',
			  '19 min post first-cleavage' => 'WBls:0000130',
			  '20 min post first-cleavage' => 'WBls:0000131',
			  '21 min post first-cleavage' => 'WBls:0000132',
			  '22 min post first-cleavage' => 'WBls:0000133',
			  '23 min post first-cleavage' => 'WBls:0000134',
			  '24 min post first-cleavage' => 'WBls:0000135',
			  '25 min post first-cleavage' => 'WBls:0000136',
			  '26 min post first-cleavage' => 'WBls:0000137',
			  '27 min post first-cleavage' => 'WBls:0000138',
			  '28 min post first-cleavage' => 'WBls:0000139',
			  '29 min post first-cleavage' => 'WBls:0000140',
			  '30 min post first-cleavage' => 'WBls:0000141',
			  '31 min post first-cleavage' => 'WBls:0000142',
			  '32 min post first-cleavage' => 'WBls:0000143',
			  '33 min post first-cleavage' => 'WBls:0000144',
			  '34 min post first-cleavage' => 'WBls:0000145',
			  '35 min post first-cleavage' => 'WBls:0000146',
			  '36 min post first-cleavage' => 'WBls:0000147',
			  '37 min post first-cleavage' => 'WBls:0000148',
			  '38 min post first-cleavage' => 'WBls:0000149',
			  '39 min post first-cleavage' => 'WBls:0000150',
			  '40 min post first-cleavage' => 'WBls:0000151',
			  '41 min post first-cleavage' => 'WBls:0000152',
			  '42 min post first-cleavage' => 'WBls:0000153',
			  '43 min post first-cleavage' => 'WBls:0000154',
			  '44 min post first-cleavage' => 'WBls:0000155',
			  '45 min post first-cleavage' => 'WBls:0000156',
			  '46 min post first-cleavage' => 'WBls:0000157',
			  '47 min post first-cleavage' => 'WBls:0000158',
			  '48 min post first-cleavage' => 'WBls:0000159',
			  '49 min post first-cleavage' => 'WBls:0000160',
			  '50 min post first-cleavage' => 'WBls:0000161',
			  '60 min post first-cleavage' => 'WBls:0000171',
			  '70 min post first-cleavage' => 'WBls:0000181',
			  '80 min post first-cleavage' => 'WBls:0000191',
			  '90 min post first-cleavage' => 'WBls:0000201',
			  '100 min post first-cleavage' => 'WBls:0000211',
			  '110 min post first-cleavage' => 'WBls:0000221',
			  '120 min post first-cleavage' => 'WBls:0000231',
			  '130 min post first-cleavage' => 'WBls:0000241',
			  '140 min post first-cleavage' => 'WBls:0000251',
			  '150 min post first-cleavage' => 'WBls:0000261',
			  '160 min post first-cleavage' => 'WBls:0000271',
			  '170 min post first-cleavage' => 'WBls:0000281',
			  '180 min post first-cleavage' => 'WBls:0000291',
			  '190 min post first-cleavage' => 'WBls:0000301',
			  '200 min post first-cleavage' => 'WBls:0000311',
			  '210 min post first-cleavage' => 'WBls:0000321',
			  '220 min post first-cleavage' => 'WBls:0000331',
                          '230 min post first-cleavage' => 'WBls:0000341',
			  '240 min post first-cleavage' => 'WBls:0000351',
                          '250 min post first-cleavage' => 'WBls:0000361',
			  '270 min post first-cleavage' => 'WBls:0000381',
                          '280 min post first-cleavage' => 'WBls:0000391',
                          '290 min post first-cleavage' => 'WBls:0000401',
			  '300 min post first-cleavage' => 'WBls:0000411',
                          '310 min post first-cleavage' => 'WBls:0000421',
                          '320 min post first-cleavage' => 'WBls:0000431',
			  '330 min post first-cleavage' => 'WBls:0000441',
                          '340 min post first-cleavage' => 'WBls:0000451',
                          '350 min post first-cleavage' => 'WBls:0000461',
			  '360 min post first-cleavage' => 'WBls:0000471',
                          '370 min post first-cleavage' => 'WBls:0000481',
                          '380 min post first-cleavage' => 'WBls:0000491',
			  '390 min post first-cleavage' => 'WBls:0000501',
                          '400 min post first-cleavage' => 'WBls:0000511',
                          '410 min post first-cleavage' => 'WBls:0000521',
			  '420 min post first-cleavage' => 'WBls:0000531',
                          '430 min post first-cleavage' => 'WBls:0000541',
                          '440 min post first-cleavage' => 'WBls:0000551',
			  '450 min post first-cleavage' => 'WBls:0000561',
                          '470 min post first-cleavage' => 'WBls:0000581',
			  '480 min post first-cleavage' => 'WBls:0000591',
			  '510 min post first-cleavage' => 'WBls:0000621',
                          '530 min post first-cleavage' => 'WBls:0000641',
			  '540 min post first-cleavage' => 'WBls:0000651',
			  '560 min post first-cleavage' => 'WBls:0000693',
			  '570 min post first-cleavage' => 'WBls:0000694',
			  '640 min post first-cleavage' => 'WBls:0000695',
			  '650 min post first-cleavage' => 'WBls:0000696',
			  '660 min post first-cleavage' => 'WBls:0000697',
			  '720 min post first-cleavage' => 'WBls:0000698',
			  '770 min post first-cleavage' => 'WBls:0000699',
			  '780 min post first-cleavage' => 'WBls:0000700',
			  '820 min post first-cleavage' => 'WBls:0000701',
			  '830 min post first-cleavage' => 'WBls:0000702',
			  '850 min post first-cleavage' => 'WBls:0000703',


			    # my extra terms
			    'early embryo' => 'WBls:0000004', # 'proliferating embryo' => 'WBls:0000004',
			    'late embryo' => 'WBls:0000021', # 'fully-elongated embryo' => 'WBls:0000021',
			    'young adult' => 'WBls:0000063', # 'newly molted young adult hermaphrodite' => 'WBls:0000063',
			    'Young Adult' => 'WBls:0000063', # 'newly molted young adult hermaphrodite' => 'WBls:0000063',
			    'mixed stage' => 'WBls:0000002', # all stages

			    'L1 larvae' => 'WBls:0000024',
			    'L2 larvae' => 'WBls:0000027',
			    'L3 larvae' => 'WBls:0000035',
			    'L4 larvae' => 'WBls:0000038',

			    'L1' => 'WBls:0000024',
			    'L2' => 'WBls:0000027',
			    'L3' => 'WBls:0000035',
			    'L4' => 'WBls:0000038',


			 );



  }

  if ($species eq 'brugia') {
    %life_stage_ontology = (
			    'all stages' => 'WBls:0000101',
			    'adult' => 'WBls:0000104',
			    'dauer larva' => 'WBls:0000032', # using the Ce code, not the nematode one

			    'sheathed microfilaria' => 'WBls:0000077',
			    'unsheathed microfilariae' => 'WBls:0000078',
			    'L1 larva' => 'WBls:0000079',
			    'L2 larva' => 'WBls:0000080',
			    'L3 larva' => 'WBls:0000081',
			    'L4 larva' => 'WBls:0000082',
			    'adult' => 'WBls:0000083',
			    'all stages' => 'WBls:0000091',
			    'embryo' => 'WBls:0000092',
			    'postembryonic' => 'WBls:0000093',
			    'early embryo' => 'WBls:0000094',
			    'middle embryo' => 'WBls:0000095',
			    'late embryo' => 'WBls:0000096',
			    'larva' => 'WBls:0000097',
			    'vector-derived L3' => 'WBls:0000098',
			    'post-infection L3' => 'WBls:0000099',
			    'young adult' => 'WBls:0000100',

			    'L1 larvae' => 'WBls:0000079',
			    'L2 larvae' => 'WBls:0000080',
			    'L3 larvae' => 'WBls:0000081',
			    'L4 larvae' => 'WBls:0000082',

			    'L1' => 'WBls:0000079',
			    'L2' => 'WBls:0000080',
			    'L3' => 'WBls:0000081',
			    'L4' => 'WBls:0000081',

			   );
  } elsif ($species eq 'sratti') {
    %life_stage_ontology = (
			    'all stages' => 'WBls:0000101',
			    'embryo' => 'WBls:0000102',
			    'postembryonic' => 'WBls:0000103',
			    'adult' => 'WBls:0000104',
			    'larva' => 'WBls:0000105',
			    'L1 larva' => 'WBls:0000106',
			    'L2 larva' => 'WBls:0000107',
			    'L3 larva' => 'WBls:0000108',
			    'L4 larva' => 'WBls:0000109',
			    'sheathed microfilaria' => 'WBls:0000110',
			    'dauer larva' => 'WBls:0000032', # using the Ce code, not the nematode one

			    'parasitic females'	 =>	'WBls:0000678',
			    'free-living females'	 =>	'WBls:0000677',
			    'free living adults'	 =>	'WBls:0000682',
			    'infective larvae (iL3)'	=>	'WBls:0000680',

			    'L1 larvae' => 'WBls:0000106',
			    'L2 larvae' => 'WBls:0000107',
			    'L3 larvae' => 'WBls:0000108',
			    'L4 larvae' => 'WBls:0000109',

			    'L1' => 'WBls:0000106',
			    'L2' => 'WBls:0000107',
			    'L3' => 'WBls:0000108',
			    'L4' => 'WBls:0000109',

			   );
    
  } elsif ($species eq 'tmuris') {
    %life_stage_ontology = (
			    'all stages' => 'WBls:0000101',
			    'embryo' => 'WBls:0000102',
			    'postembryonic' => 'WBls:0000103',
			    'adult' => 'WBls:0000104',
			    'larva' => 'WBls:0000105',
			    'L1 larva' => 'WBls:0000106',
			    'L2 larva' => 'WBls:0000107',
			    'L3 larva' => 'WBls:0000108',
			    'L4 larva' => 'WBls:0000109',
			    'sheathed microfilaria' => 'WBls:0000110',
			    'dauer larva' => 'WBls:0000032', # using the Ce code, not the nematode one

			    'intestinal'                => 'WBls:0000721',
			    'adult intestinal'          => 'WBls:0000721',
			    'lumen'                     => 'WBls:0000722',
			    'lumenal'                   => 'WBls:0000722',
			    'adult lumenal'             => 'WBls:0000722',

			    'L1 larvae' => 'WBls:0000106',
			    'L2 larvae' => 'WBls:0000107',
			    'L3 larvae' => 'WBls:0000108',
			    'L4 larvae' => 'WBls:0000109',

			    'L1' => 'WBls:0000106',
			    'L2' => 'WBls:0000107',
			    'L3' => 'WBls:0000108',
			    'L4' => 'WBls:0000109',

			    );
  }

# see if we can find a life_stage already seen in this study
 NAMES:
  foreach my $test_name (@{$life_stage_session{names_found}}) {
    foreach my $line ($run_alias, $experiment_title, $study_alias, $study_title) {
      if ($line =~ /\b$test_name\b/) {
	$candidate_life_stage = $test_name;
	$why = 'in description';
	last NAMES;
      }
    }    
  }
  
  # search for life_stage names from the ontology in the descriptions
  # search for the longest names first to get the most exact term matching first
 ONTOLOGIES:
  foreach my $life_stage_term (sort {length $b <=> length $a} (keys %life_stage_ontology)) {
    foreach my $line ($run_alias, $experiment_title, $study_alias, $study_title) {
      if ($line =~ /\b$life_stage_term\b/) {
	$candidate_life_stage = $life_stage_term;
	$why = 'in description';
	last ONTOLOGIES;
      }
    }    
  }

  if (!defined  $candidate_life_stage || $candidate_life_stage eq '') {
    if (defined $ENA_dev_stage) {
      $candidate_life_stage = $ENA_dev_stage;
      $why = 'from ENA';
    }
  }

  # assume we are using the same life_stage as in the last experiment
  if (!defined $candidate_life_stage && defined $life_stage_session{previous}) {
    $candidate_life_stage = $life_stage_session{previous};
    $why = 'repeat of last expt';
  }

  # else if this is the first time in, use the default
  if (!defined $candidate_life_stage || $candidate_life_stage eq '') {
    $candidate_life_stage = $default;
    $why = 'default';
  }

  # confirm and update status
  do {
    print "Life_Stage [$candidate_life_stage] ($why) > ";
    my $input =  <STDIN>;
    chomp ($input);
    if ($input eq '') {
      $life_stage = $life_stage_ontology{$candidate_life_stage};
      $life_stage_name = $candidate_life_stage;
    } else {
      $life_stage = $life_stage_ontology{$input};
      if (!defined $life_stage) {
	print "WARNING: $input is not known in the ontology.\n"
      } else {
	push @{$life_stage_session{names_found}}, $input;
      $life_stage_name = $input;
      }
    }
  } until (defined $life_stage);
  $life_stage_session{previous} = $life_stage_name;

  return ($life_stage, $life_stage_name);
}


###################################################################

sub get_study_display {

  my $display='';
  # get a brief title and citation to display in the website
  while ($display eq '') {
    print "Study display title like 'Between species comparison (Morazavi 2013)' > ";
    $display =  <STDIN>;
    chomp ($display);
  }
  return $display;

}

###################################################################

# look through the .ini config files and find studies with a pubmed ID
# that have experiments that have not been curated and which are
# TRANSCRIPTOMIC, etc. Display details of each Study.

sub suggest {

  my %expt_hashref;
  my %study_hashref;
  my $count; # number of experiments in this study that could be curated

  my $studies_ini = $RNASeq->read_all_studies_config();
  if (defined $studies_ini) {
    $RNASeq->convert_ini_to_hashref($studies_ini, \%study_hashref);
    foreach my $study (keys %study_hashref) {
      %expt_hashref = ();
      $count = 0;
      if (exists $study_hashref{ $study}{ignore} ) {next}
      if (!exists $study_hashref{$study}{pubmed}) {next}
      my $pubmed = $study_hashref{$study}{pubmed};
      my $wbpaper = (exists $study_hashref{$study}{wbpaper}) ? $study_hashref{$study}{wbpaper} : '';
      my $experiments_ini = $RNASeq->read_experiments_from_study_config($study);
      $RNASeq->convert_ini_to_hashref($experiments_ini, \%expt_hashref);
      foreach my $expt (keys %expt_hashref) {
	if (!exists $expt_hashref{$expt}{analysis} ) {
	  $count++;
	}
      }
      if ($count) {
	print "$study: $count experiments pubmed: $pubmed wbpaper: $wbpaper\n";
      }
    }
  }
  print "Inform Kimberley of any pubmed papers that do not have a WBPaper ID\n";
}

###################################################################
sub get_anatomy {
  my $def = &write_anatomy_def();
  my $query = $wormbase->table_maker_query($database, $def);

  my %anatomy;

  while(<$query>) {
    chomp;
    s/\"//g;
    next if (/acedb/ or /\/\//);
    next if /^\s*$/;

    my ($annot_name, $WBbt) = split(/\t/, $_);

    $anatomy{$annot_name} = $WBbt;
  }

  return %anatomy;
}
###################################################################
sub write_anatomy_def {
  my $txt = <<"ENDE";
Sortcolumn 1

Colonne 1 
Width 64 
Optional 
Visible 
Class 
Class Anatomy_name 
From 1 
 
Colonne 2 
Width 64 
Mandatory 
Visible 
Class 
Class Anatomy_term 
From 1 
Tag Name_for_anatomy_term 
ENDE

  return &write_tm_def("current_go_annots", $txt);
}
###################################################################
sub write_tm_def {
  my ($fname, $string) = @_;
  
  my $file = "/tmp/$fname.def";

  open(my $fh, ">$file") or $log->log_and_die("Could not open $fname for writing\n");
  print $fh $string;
  close($fh) or $log->log_and_die("Could not close $fname after writing\n");

  return $file;
}

###################################################################
sub query_ENA {
  my ($experiment) = @_;

my $query = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${experiment}&result=read_experiment&fields=dev_stage,sex,strain,temperature,tissue_type,description";


  open (DATA, "wget -q -O - '$query' |") || die("RNASeq: Can't get information on SRA entry $experiment in read_accession()\n");

  my $line_count=0;
  while (my $line = <DATA>) {
    if (++$line_count == 1) {next;} # skip the title line
    chomp $line;
    
    my @f = split /\t/, $line;

    return @f; # dev_stage, sex,strain, temperature, tissue_type, description
  }


}
