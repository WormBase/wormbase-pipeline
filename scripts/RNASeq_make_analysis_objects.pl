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

my %default_strain = (
		      'elegans'   => 'N2',
		      'brenneri'  => 'PB2801',
		      'briggsae'  => 'AF16',
		      'brugia'    => 'FR3',
		      'sratti'    => 'ED321',
		      'japonica'  => 'DF5081',
		      'ovolvulus' => 'O_volvulus_Cameroon_isolate',
		      'remanei'   => 'SB146',
		      'pristionchus'   => 'PS312',
		     );


open(OUT, ">$output") || $log->log_and_die("Can't open $output\n");


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
      if (exists $hashref->{$expt}{analysis} && !$redo) {print "\n$expt has an Analysis already (",$hashref->{$expt}{analysis},") - skipping\n"; next}

      # get the study_title, study_alias, experiment_title, library_name, run_alias and do things with them
      $study_title      = $hashref->{$expt}{study_title};
      $study_alias      = $hashref->{$expt}{study_alias};
      $experiment_title = $hashref->{$expt}{experiment_title};
      $library_name     = $hashref->{$expt}{library_name};
      $run_alias        = $hashref->{$expt}{run_alias};

      print "\nExperiment: $expt\n";
      print "Study title: $study_title\n" if ($study_title ne '');
      print "Study alias: $study_alias\n" if ($study_alias ne '');
      print "Experiment title: $experiment_title\n" if ($experiment_title ne '');
      print "Library name: $library_name\n" if ($library_name ne '');
      print "Run alias: $run_alias\n" if ($run_alias ne '');
      print "\n";


      # look for temperature
      # default is nothing
      $temperature = get_temperature();

      # look for genotype
      # default is nothing
      $genotype = get_genotype();

      # look for treatment
      # default is nothing
      $treatment = get_treatment();

      # look for strain
      # default is initially N2 in elegans
      if (!defined $genotype || $genotype eq "") {
	$strain = get_strain();
      } else {
	$strain = $default_strain{$species};
      }

      # look for sex
      # one probable synonym is 'him-8:Male'
      # synonym / lookup will be different in different species
      $sex = get_sex();

      # look for tissue
      # use synonym / lookup hash with recently found synonyms being checked first in the title they were last found in
      # default is initially whole organism
      ($tissue, $tissue_name) = get_tissue();
      
      # look for life stage
      # use synonym / lookup hash with recently found synonyms being checked first in the title they were last found in
      # synonym / lookup will be different in different species
      ($life_stage, $life_stage_name) = get_life_stage();

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

  print "Temperature > ";
  my $error;
  my $input;
  do {
    $input =  <STDIN>;
    chomp ($input);
    $error = 0;
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

  my $tissue = undef; 
  my $tissue_name = undef; 
  my $candidate_tissue = undef;
  my $line;
  my $why = '';

  my $default = 'organism'; # WBbt:0007833'

  my %tissue_ontology = (
			 'organism' => 'WBbt:0007833',

			 'M2 neuron' => 'WBbt:0003634',
			 'MC neuron' => 'WBbt:0003638',
			 'marginal cell neuron' => 'WBbt:0003638',
			 'I2 neuron'	 => 'WBbt:0003645',
			 'I1 neuron'	 => 'WBbt:0003649',
			 'MI neuron' => 'WBbt:0003664',
			 'pharyngeal motor-interneuron' => 'WBbt:0003664',
			 'NSM' => 'WBbt:0003666',
			 'neurosecretory-motor neuron' => 'WBbt:0003666',
			 'pharyngeal interneuron'	 => 'WBbt:0003668',
			 'gland cell'	 => 'WBbt:0003670',
			 'epithelial cell'	 => 'WBbt:0003672',
			 'epithelium' => 'WBbt:0003672',
			 'marginal cell'	 => 'WBbt:0003673',
			 'muscle cell'	 => 'WBbt:0003675',
			 'pharyngeal motor neuron'	 => 'WBbt:0003677',
			 'neuron'	 => 'WBbt:0003679',
			 'neurone' => 'WBbt:0003679',
			 'pharynx'	 => 'WBbt:0003681',
			 'esophagus' => 'WBbt:0003681',
			 'metacorpus'	 => 'WBbt:0003711',
			 'procorpus'	 => 'WBbt:0003713',
			 'terminal bulb'	 => 'WBbt:0003732',
			 'corpus'	 => 'WBbt:0003733',
			 'isthmus'	 => 'WBbt:0003734',
			 'M3 neuron'	 => 'WBbt:0003754',
			 'M3' => 'WBbt:0003754',
			 'organ'	 => 'WBbt:0003760',
			 'caudal longitudinal muscle right'	 => 'WBbt:0003783',
			 'cdlR' => 'WBbt:0003783',
			 'caudal longitudinal muscle right'	 => 'WBbt:0003783',
			 'posterior inner longitudinal muscle right'	 => 'WBbt:0003784',
			 'posterior inner longitudinal muscle right'	 => 'WBbt:0003784',
			 'anterior inner longitudinal muscle right'	 => 'WBbt:0003785',
			 'anterior inner longitudinal muscle right'	 => 'WBbt:0003785',
			 'posterior outer longitudinal muscle right'	 => 'WBbt:0003786',
			 'posterior outer longitudinal muscle right' => 'WBbt:0003786',
			 'anterior outer longitudinal muscle right' => 'WBbt:0003787',
			 'anterior outer longitudinal muscle right' => 'WBbt:0003787',
			 'caudal longitudinal muscle left' => 'WBbt:0003788',
			 'caudal longitudinal muscle left' => 'WBbt:0003788',
			 'posterior inner longitudinal muscle left' => 'WBbt:0003789',
			 'posterior inner longitudinal muscle left' => 'WBbt:0003789',
			 'anterior inner longitudinal muscle left' => 'WBbt:0003790',
			 'anterior inner longitudinal muscle left' => 'WBbt:0003790',
			 'posterior outer longitudinal muscle left' => 'WBbt:0003791',
			 'posterior outer longitudinal muscle left' => 'WBbt:0003791',
			 'anterior outer longitudinal muscle left' => 'WBbt:0003792',
			 'anterior outer longitudinal muscle left' => 'WBbt:0003792',
			 'gubernacular retractor muscle right' => 'WBbt:0003793',
			 'gubernacular retractor muscle right' => 'WBbt:0003793',
			 'gubernacular retractor muscle left' => 'WBbt:0003794',
			 'gubernacular retractor muscle left' => 'WBbt:0003794',
			 'gubernacular erector right' => 'WBbt:0003795',
			 'gubernacular erector right' => 'WBbt:0003795',
			 'gubernacular erector left' => 'WBbt:0003796',
			 'gubernacular erector left' => 'WBbt:0003796',
			 'AVF'	 => 'WBbt:0003851',
			 'AVFL/R' => 'WBbt:0003851',
			 'AB'	 => 'WBbt:0004015',
			 'AB blastomere' => 'WBbt:0004015',
			 'anal depressor muscle'	 => 'WBbt:0004292',
			 'dilator muscle' => 'WBbt:0004292',
			 'anal depressor muscle'	 => 'WBbt:0004292',
			 'PDA'	 => 'WBbt:0004386',
			 'M5 neuron'	 => 'WBbt:0004465',
			 'M4 neuron'	 => 'WBbt:0004467',
			 'M1 neuron'	 => 'WBbt:0004488',
			 'M cell'	 => 'WBbt:0004489',
			 'K cell'	 => 'WBbt:0004497',
			 'Anchor cell'	 => 'WBbt:0004522',
			 'AC' => 'WBbt:0004522',
			 'Anchor cell'	 => 'WBbt:0004522',
			 'excretory socket cell' => 'WBbt:0004534',
			 'excretory socket cell' => 'WBbt:0004534',
			 'excretory duct cell' => 'WBbt:0004540',
			 'Y cell'	 => 'WBbt:0004578',
			 'postembryonic C cell' => 'WBbt:0004578',
			 'I6 neuron'	 => 'WBbt:0004739',
			 'I5 neuron'	 => 'WBbt:0004740',
			 'I4 neuron'	 => 'WBbt:0004741',
			 'I3 neuron'	 => 'WBbt:0004742',
			 'dorsal spicule protractor right' => 'WBbt:0004885',
			 'dorsal spicule protractor right' => 'WBbt:0004885',
			 'dorsal spicule protractor left' => 'WBbt:0004886',
			 'dorsal spicule protractor left' => 'WBbt:0004886',
			 'dorsal spicule retractor muscle right' => 'WBbt:0004907',
			 'dorsal spicule retractor muscle right' => 'WBbt:0004907',
			 'ventral spicule retractor muscle right' => 'WBbt:0004908',
			 'ventral spicule retractor muscle right' => 'WBbt:0004908',
			 'ventral spicule protrctor muscle right' => 'WBbt:0004909',
			 'ventral spicule protrctor muscle right' => 'WBbt:0004909',
			 'dorsal spicule retractor muscle left' => 'WBbt:0004910',
			 'dorsal spicule retractor muscle left' => 'WBbt:0004910',
			 'ventral spicule retractor muscle left' => 'WBbt:0004911',
			 'ventral spicule retractor muscle left' => 'WBbt:0004911',
			 'ventral spicule protractor muscle left' => 'WBbt:0004912',
			 'ventral spicule protractor muscle left' => 'WBbt:0004912',
			 'posterior oblique muscle right' => 'WBbt:0004913',
			 'posterior oblique muscle right' => 'WBbt:0004913',
			 'anterior oblique muscle right' => 'WBbt:0004914',
			 'anterior oblique muscle right' => 'WBbt:0004914',
			 'posterior oblique muscle left' => 'WBbt:0004915',
			 'posterior oblique muscle left' => 'WBbt:0004915',
			 'anterior oblique muscle left' => 'WBbt:0004916',
			 'anterior oblique muscle left' => 'WBbt:0004916',
			 'U cell'	 => 'WBbt:0004942',
			 'postembryonic E cell' => 'WBbt:0004942',
			 'RIPR' => 'WBbt:0005049',
			 'linker cell' => 'WBbt:0005062',
			 'lumbar left ganglion neuron'	 => 'WBbt:0005097',
			 'lumbar ganglion left'	 => 'WBbt:0005098',
			 'lateral ganglion right neuron'	 => 'WBbt:0005100',
			 'lateral ganglion right'	 => 'WBbt:0005101',
			 'lateral ganglion left neuron'	 => 'WBbt:0005102',
			 'lateral ganglion left'	 => 'WBbt:0005103',
			 'lateral ganglion'	 => 'WBbt:0005105',
			 'lateral ganglia' => 'WBbt:0005105',
			 'labial sensillum'	 => 'WBbt:0005107',
			 'interneuron'	 => 'WBbt:0005113',
			 'inner labial sensillum'	 => 'WBbt:0005116',
			 'inner labial neuron'	 => 'WBbt:0005117',
			 'IL neuron' => 'WBbt:0005117',
			 'IL2 neuron'	 => 'WBbt:0005118',
			 'IL1 neuron'	 => 'WBbt:0005119',
			 'IL sensillum ventral right'	 => 'WBbt:0005120',
			 'IL sensillum ventral left'	 => 'WBbt:0005121',
			 'IL sensillum right'	 => 'WBbt:0005122',
			 'IL sensillum left'	 => 'WBbt:0005123',
			 'IL sensillum dorsal right'	 => 'WBbt:0005124',
			 'IL sensillum dorsal left'	 => 'WBbt:0005126',
			 'head ganglion'	 => 'WBbt:0005135',
			 'gonad'	 => 'WBbt:0005175',
			 'GLR'	 => 'WBbt:0005176',
			 'GLR scaffold cell' => 'WBbt:0005176',
			 'HO neuron'	 => 'WBbt:0005177',
			 'hermaphrodite gonad'	 => 'WBbt:0005178',
			 'ovotestis' => 'WBbt:0005178',
			 'ganglion'	 => 'WBbt:0005189',
			 'GABAergic neuron'	 => 'WBbt:0005190',
			 'dorso-rectal ganglion'	 => 'WBbt:0005212',
			 'dorsal ganglion'	 => 'WBbt:0005214',
			 'URA'	 => 'WBbt:0005225',
			 'touch receptor neuron' => 'WBbt:0005237',
			 'touch receptor neuron'	 => 'WBbt:0005237',
			 'microtubule cell' => 'WBbt:0005237',
			 'touch receptor neuron'	 => 'WBbt:0005237',
			 'touch cell' => 'WBbt:0005237',
			 'touch receptor neuron'	 => 'WBbt:0005237',
			 'touch receptor' => 'WBbt:0005237',
			 'CEP'	 => 'WBbt:0005244',
			 'CEM'	 => 'WBbt:0005246',
			 'buccal cavity'	 => 'WBbt:0005255',
			 'ventral ganglion (post)'	 => 'WBbt:0005264',
			 'ventral ganglion (ant)'	 => 'WBbt:0005266',
			 'DD neuron'	 => 'WBbt:0005270',
			 'DB neuron'	 => 'WBbt:0005274',
			 'DA neuron'	 => 'WBbt:0005278',
			 'ventral ganglion'	 => 'WBbt:0005298',
			 'ventral ganglia' => 'WBbt:0005298',
			 'ventral cord neuron'	 => 'WBbt:0005300',
			 'ventral cord motoneuron' => 'WBbt:0005300',
			 'VD neuron'	 => 'WBbt:0005303',
			 'VC neuron'	 => 'WBbt:0005304',
			 'spicule socket cell'	 => 'WBbt:0005309',
			 'CP neuron'	 => 'WBbt:0005310',
			 'spicule sheath cell'	 => 'WBbt:0005311',
			 'copulatory spicule'	 => 'WBbt:0005312',
			 'copulatory spicule right'	 => 'WBbt:0005314',
			 'copulatory spicule left'	 => 'WBbt:0005316',
			 'spermatheca'	 => 'WBbt:0005319',
			 'body ganglion'	 => 'WBbt:0005332',
			 'AS neuron'	 => 'WBbt:0005334',
			 'VB neuron'	 => 'WBbt:0005336',
			 'vas deferens' => 'WBbt:0005337',
			 'vas deferens'	 => 'WBbt:0005337',
			 'sperm duct' => 'WBbt:0005337',
			 'VA neuron'	 => 'WBbt:0005339',
			 'uterine muscle'	 => 'WBbt:0005342',
			 'URY'	 => 'WBbt:0005346',
			 'cloacal ganglion'	 => 'WBbt:0005348',
			 'SMD'	 => 'WBbt:0005353',
			 'SMB'	 => 'WBbt:0005355',
			 'SIB'	 => 'WBbt:0005359',
			 'SIA'	 => 'WBbt:0005361',
			 'anus'	 => 'WBbt:0005364',
			 'anus interface' => 'WBbt:0005364',
			 'anterior spermatheca' => 'WBbt:0005368',
			 'anterior ganglion (post)'	 => 'WBbt:0005371',
			 'anterior ganglion (ant)'	 => 'WBbt:0005372',
			 'anterior hypodermis'	 => 'WBbt:0005373',
			 'head hypodermis' => 'WBbt:0005373',
			 'anterior gonad arm'	 => 'WBbt:0005374',
			 'anterior ganglion'	 => 'WBbt:0005375',
			 'amphid sensillum'	 => 'WBbt:0005391',
			 'amphid' => 'WBbt:0005391',
			 'amphid right sensillum'	 => 'WBbt:0005392',
			 'amphid process'	 => 'WBbt:0005393',
			 'amphid neuron'	 => 'WBbt:0005394',
			 'amphid sensory neuron' => 'WBbt:0005394',
			 'amphid left sensillum'	 => 'WBbt:0005395',
			 'SAB'	 => 'WBbt:0005396',
			 'SAA'	 => 'WBbt:0005397',
			 'RME'	 => 'WBbt:0005399',
			 'RMD'	 => 'WBbt:0005400',
			 'retrovesicular ganglion neuron'	 => 'WBbt:0005403',
			 'retrovesicular ganglion' => 'WBbt:0005403',
			 'retrovesicular ganglion neuron'	 => 'WBbt:0005403',
			 'RVG neuron' => 'WBbt:0005403',
			 'ALM'	 => 'WBbt:0005406',
			 'motor neuron'	 => 'WBbt:0005409',
			 'motoneuron' => 'WBbt:0005409',
			 'AIY'	 => 'WBbt:0005413',
			 'ADE'	 => 'WBbt:0005415',
			 'phasmid sensillum'	 => 'WBbt:0005425',
			 'phasmid right sensillum'	 => 'WBbt:0005427',
			 'phasmid left sensillum'	 => 'WBbt:0005431',
			 'pharyngeal nerve process'	 => 'WBbt:0005437',
			 'AB cells' => 'WBbt:0005437',
			 'pharyngeal neuron'	 => 'WBbt:0005439',
			 'pharyngeal nervous system'	 => 'WBbt:0005440',
			 'pharyngeal nerve'	 => 'WBbt:0005441',
			 'preanal ganglion neuron'	 => 'WBbt:0005447',
			 'preanal ganglion neurons' => 'WBbt:0005447',
			 'preanal ganglion'	 => 'WBbt:0005448',
			 'pharyngeal muscle cell'	 => 'WBbt:0005451',
			 'pharyngeal epithelial cell'	 => 'WBbt:0005459',
			 'pharyngeal epithelium' => 'WBbt:0005459',
			 'pharyngeal cell'	 => 'WBbt:0005460',
			 'Posterior spermatheca' => 'WBbt:0005462',
			 'posterior lateral right ganglion'	 => 'WBbt:0005463',
			 'posterior lateral left ganglion'	 => 'WBbt:0005464',
			 'posterior lateral ganglion'	 => 'WBbt:0005465',
			 'posterior lateral ganglia' => 'WBbt:0005465',
			 'posterior gonad arm'	 => 'WBbt:0005466',
			 'postdeirid sensillum'	 => 'WBbt:0005471',
			 'PC neuron'	 => 'WBbt:0005473',
			 'postcloacal sensillum right'	 => 'WBbt:0005485',
			 'postcloacal sensillum left'	 => 'WBbt:0005486',
			 'postcloacal sensillum'	 => 'WBbt:0005487',
			 'PLM'	 => 'WBbt:0005490',
			 'outer labial sensillum'	 => 'WBbt:0005501',
			 'outer labial quadrant sensillum'	 => 'WBbt:0005502',
			 'outer labial lateral sensillum'	 => 'WBbt:0005503',
			 'OL quadrant ventral right sensillum'	 => 'WBbt:0005505',
			 'OL quadrant ventral left sensillum'	 => 'WBbt:0005506',
			 'OL quadrant dorsal right sensillum'	 => 'WBbt:0005507',
			 'OL quadrant dorsal left sensillum'	 => 'WBbt:0005508',
			 'OL lateral right sensillum'	 => 'WBbt:0005509',
			 'OL lateral left sensillum'	 => 'WBbt:0005510',
			 'lumbar right ganglion neuron'	 => 'WBbt:0005600',
			 'lumbar ganglion right'	 => 'WBbt:0005601',
			 'retrovesicular ganglion'	 => 'WBbt:0005656',
			 'AVFL'	 => 'WBbt:0005657',
			 'AVFR'	 => 'WBbt:0005658',
			 'PHBR'	 => 'WBbt:0005659',
			 'PHRB' => 'WBbt:0005659',
			 'ADF'	 => 'WBbt:0005660',
			 'ADL'	 => 'WBbt:0005661',
			 'AFD'	 => 'WBbt:0005662',
			 'ASE'	 => 'WBbt:0005663',
			 'ASG'	 => 'WBbt:0005664',
			 'ASH'	 => 'WBbt:0005665',
			 'ASI'	 => 'WBbt:0005666',
			 'ASJ'	 => 'WBbt:0005667',
			 'ASK'	 => 'WBbt:0005668',
			 'AWA'	 => 'WBbt:0005670',
			 'AWB'	 => 'WBbt:0005671',
			 'AWC'	 => 'WBbt:0005672',
			 'epithelial system'	 => 'WBbt:0005730',
			 'rectal lumen'	 => 'WBbt:0005731',
			 'extracellular component'	 => 'WBbt:0005732',
			 'hypodermis' => 'WBbt:0005733',
			 'epidermis' => 'WBbt:0005733',
			 'hyp7 syncytium' => 'WBbt:0005734',
			 'nervous system'	 => 'WBbt:0005735',
			 'excretory system'	 => 'WBbt:0005736',
			 'muscular system'	 => 'WBbt:0005737',
			 'body region'	 => 'WBbt:0005738',
			 'head'	 => 'WBbt:0005739',
			 'midbody'	 => 'WBbt:0005740',
			 'tail'	 => 'WBbt:0005741',
			 'bodywall'	 => 'WBbt:0005742',
			 'digestive tract'	 => 'WBbt:0005743',
			 'reproductive tract'	 => 'WBbt:0005744',
			 'pseudocoelom'	 => 'WBbt:0005745',
			 'body cavity' => 'WBbt:0005745',
			 'Organ system'	 => 'WBbt:0005746',
			 'reproductive system'	 => 'WBbt:0005747',
			 'alimentary system'	 => 'WBbt:0005748',
			 'coelomic system'	 => 'WBbt:0005749',
			 'socket cell'	 => 'WBbt:0005750',
			 'coelomocyte'	 => 'WBbt:0005751',
			 'seam cell'	 => 'WBbt:0005753',
			 'lateral hypodermis' => 'WBbt:0005753',
			 'interfacial epithelial cell'	 => 'WBbt:0005754',
			 'interfacial cell' => 'WBbt:0005754',
			 'interfacial epithelial cell'	 => 'WBbt:0005754',
			 'interfacial epithelium' => 'WBbt:0005754',
			 'interfacial epithelial cell'	 => 'WBbt:0005754',
			 'transitional epithelium' => 'WBbt:0005754',
			 'cuticle'	 => 'WBbt:0005755',
			 'exoskeleton' => 'WBbt:0005755',
			 'basal lamina'	 => 'WBbt:0005756',
			 'basement membrane' => 'WBbt:0005756',
			 'basal lamina'	 => 'WBbt:0005756',
			 'ECM' => 'WBbt:0005756',
			 'basal lamina'	 => 'WBbt:0005756',
			 'extracellular matrix' => 'WBbt:0005756',
			 'basal lamina'	 => 'WBbt:0005756',
			 'lamina' => 'WBbt:0005756',
			 'male-specific'	 => 'WBbt:0005757',
			 'hermaphrodite-specific'	 => 'WBbt:0005758',
			 'sensory neuron'	 => 'WBbt:0005759',
			 'somatic nervous system'	 => 'WBbt:0005760',
			 'accessory cell'	 => 'WBbt:0005762',
			 'support cell' => 'WBbt:0005762',
			 'pharyngeal-intestinal valve'	 => 'WBbt:0005767',
			 'cardia' => 'WBbt:0005767',
			 'pharyngeal-intestinal valve'	 => 'WBbt:0005767',
			 'esophago-intestinal valve' => 'WBbt:0005767',
			 'pharyngeal-intestinal valve'	 => 'WBbt:0005767',
			 'pharyngeal valve' => 'WBbt:0005767',
			 'pharyngeal-intestinal valve'	 => 'WBbt:0005767',
			 'pharyngo-intestinal valve' => 'WBbt:0005767',
			 'axis'	 => 'WBbt:0005768',
			 'anterior-posterior'	 => 'WBbt:0005769',
			 'dorsal-ventral'	 => 'WBbt:0005770',
			 'left-right'	 => 'WBbt:0005771',
			 'intestine'	 => 'WBbt:0005772',
			 'gut' => 'WBbt:0005772',
			 'intestine'	 => 'WBbt:0005772',
			 'mesenteron' => 'WBbt:0005772',
			 'rectum'	 => 'WBbt:0005773',
			 'cloaca'	 => 'WBbt:0005774',
			 'excretory canal'	 => 'WBbt:0005775',
			 'canal' => 'WBbt:0005775',
			 'excretory gland cell' => 'WBbt:0005776',
			 'Excretory duct'	 => 'WBbt:0005777',
			 'excretory pore'	 => 'WBbt:0005778',
			 'pore' => 'WBbt:0005778',
			 'striated muscle'	 => 'WBbt:0005779',
			 'non-striated muscle'	 => 'WBbt:0005780',
			 'smooth muscle'	 => 'WBbt:0005781',
			 'germ line'	 => 'WBbt:0005784',
			 'germline' => 'WBbt:0005784',
			 'somatic gonad'	 => 'WBbt:0005785',
			 'somatic germline' => 'WBbt:0005785',
			 'muscle of the reproductive system'	 => 'WBbt:0005786',
			 'Pharyngeal gland cell'	 => 'WBbt:0005788',
			 'pharyngeal segment'	 => 'WBbt:0005789',
			 'pharyngeal lumen'	 => 'WBbt:0005790',
			 'intestinal lumen'	 => 'WBbt:0005791',
			 'intestinal cell'	 => 'WBbt:0005792',
			 'arcade cell'	 => 'WBbt:0005793',
			 'anterior arcade cell' => 'WBbt:0005794',
			 'posterior arcade cell' => 'WBbt:0005795',
			 'intestinal muscle'	 => 'WBbt:0005796',
			 'stomato-intestinal muscle' => 'WBbt:0005796',
			 'rectal valve cell'	 => 'WBbt:0005797',
			 'intestinal-rectal valve cell' => 'WBbt:0005797',
			 'anal sphincter muscle' => 'WBbt:0005798',
			 'rectal gland cell'	 => 'WBbt:0005799',
			 'rectal epithelium'	 => 'WBbt:0005800',
			 'dorsal-rectal ganglion neuron'	 => 'WBbt:0005801',
			 'dorso-rectal ganglion neuron' => 'WBbt:0005801',
			 'rectal muscle' => 'WBbt:0005803',
			 'anal muscle' => 'WBbt:0005803',
			 'cloacal rectal valve cell'	 => 'WBbt:0005804',
			 'cloacal sphincter muscle'	 => 'WBbt:0005805',
			 'cloacal rectal gland cell'	 => 'WBbt:0005806',
			 'cloacal neuron'	 => 'WBbt:0005807',
			 'cloacal epithelium'	 => 'WBbt:0005808',
			 'cloacal muscle'	 => 'WBbt:0005809',
			 'cloacal dilator muscle'	 => 'WBbt:0005810',
			 'neuronal sheath cell'	 => 'WBbt:0005811',
			 'pocket cell' => 'WBbt:0005811',
			 'excretory cell'	 => 'WBbt:0005812',
			 'excretory canal cell' => 'WBbt:0005812',
			 'excretory cell'	 => 'WBbt:0005812',
			 'body wall musculature'	 => 'WBbt:0005813',
			 'body muscle' => 'WBbt:0005813',
			 'diagonal muscle' => 'WBbt:0005815',
			 'diagonal muscle'	 => 'WBbt:0005815',
			 'male tail diagonal muscle' => 'WBbt:0005815',
			 'diagonal muscle' => 'WBbt:0005815',
			 'dorsal left quadrant body wall muscle'	 => 'WBbt:0005816',
			 'dorsal right quadrant body wall muscle'	 => 'WBbt:0005817',
			 'ventral left quadrant body wall muscle'	 => 'WBbt:0005818',
			 'ventral right quadrant body wall muscle'	 => 'WBbt:0005819',
			 'alimentary muscle'	 => 'WBbt:0005820',
			 'vulval muscle'	 => 'WBbt:0005821',
			 'anterior oblique muscle'	 => 'WBbt:0005822',
			 'posterior oblique muscle'	 => 'WBbt:0005823',
			 'gubernacular erector muscle'	 => 'WBbt:0005824',
			 'gubernacular retractor muscle'	 => 'WBbt:0005825',
			 'spicule protractor muscle'	 => 'WBbt:0005826',
			 'spicule retractor muscle'	 => 'WBbt:0005827',
			 'gonadal sheath cell'	 => 'WBbt:0005828',
			 'ventral nerve cord'	 => 'WBbt:0005829',
			 'ventral cord' => 'WBbt:0005829',
			 'lumbar ganglion'	 => 'WBbt:0005830',
			 'lumbar lateral ganglia' => 'WBbt:0005830',
			 'lumbar ganglion'	 => 'WBbt:0005830',
			 'lumbar lateral ganglion' => 'WBbt:0005830',
			 'spicule neuron'	 => 'WBbt:0005831',
			 'SP neuron' => 'WBbt:0005831',
			 'AWC-ON'	 => 'WBbt:0005832',
			 'AWC-OFF'	 => 'WBbt:0005833',
			 'RIB'	 => 'WBbt:0005834',
			 'RIM'	 => 'WBbt:0005835',
			 'AIZ'	 => 'WBbt:0005836',
			 'chemosensory neuron'	 => 'WBbt:0005837',
			 'thermosensory neuron'	 => 'WBbt:0005838',
			 'odorsensory neuron'	 => 'WBbt:0005839',
			 'PVC'	 => 'WBbt:0005840',
			 'AVB'	 => 'WBbt:0005841',
			 'AVA'	 => 'WBbt:0005842',
			 'eggshell'	 => 'WBbt:0005843',
			 'dopaminergic neuron'	 => 'WBbt:0006746',
			 'PDE'	 => 'WBbt:0006747',
			 'postdeirid neuron' => 'WBbt:0006747',
			 'vulva'	 => 'WBbt:0006748',
			 'nerve ring'	 => 'WBbt:0006749',
			 'circumpharyngeal nerve ring' => 'WBbt:0006749',
			 'dorsal nerve cord'	 => 'WBbt:0006750',
			 'dorsal cord' => 'WBbt:0006750',
			 'head neuron'	 => 'WBbt:0006751',
			 'somatic neuron'	 => 'WBbt:0006752',
			 'body neuron' => 'WBbt:0006752',
			 'phasmid neuron'	 => 'WBbt:0006753',
			 'amphid sheath cell'	 => 'WBbt:0006754',
			 'phasmid sheath cell'	 => 'WBbt:0006755',
			 'spermathecal-uterine junction' => 'WBbt:0006756',
			 'spermathecal-uterine junction'	 => 'WBbt:0006756',
			 'spermatheca-uterus valve' => 'WBbt:0006756',
			 'spermathecal-uterine junction'	 => 'WBbt:0006756',
			 'spermathecal valve' => 'WBbt:0006756',
			 'anterior spermathecal-uterine junction' => 'WBbt:0006757',
			 'posterior spermathecal-uterine junction' => 'WBbt:0006758',
			 'tail neuron'	 => 'WBbt:0006759',
			 'uterus'	 => 'WBbt:0006760',
			 'head muscle'	 => 'WBbt:0006761',
			 'head body wall muscle' => 'WBbt:0006761',
			 'vulA'	 => 'WBbt:0006762',
			 'vulB1'	 => 'WBbt:0006763',
			 'vulB2'	 => 'WBbt:0006764',
			 'vulC'	 => 'WBbt:0006765',
			 'vulD'	 => 'WBbt:0006766',
			 'vulE'	 => 'WBbt:0006767',
			 'vulF'	 => 'WBbt:0006768',
			 'lateral nerve cord'	 => 'WBbt:0006769',
			 'P1'	 => 'WBbt:0006770',
			 'P1 post-embryonic blast cell' => 'WBbt:0006770',
			 'P2'	 => 'WBbt:0006771',
			 'P3'	 => 'WBbt:0006772',
			 'P4'	 => 'WBbt:0006773',
			 'P5'	 => 'WBbt:0006774',
			 'P6'	 => 'WBbt:0006775',
			 'P7'	 => 'WBbt:0006776',
			 'P8'	 => 'WBbt:0006777',
			 'P9'	 => 'WBbt:0006778',
			 'P10'	 => 'WBbt:0006779',
			 'ventral uterine precursor'	 => 'WBbt:0006780',
			 'VU' => 'WBbt:0006780',
			 'dorsal uterine precursor'	 => 'WBbt:0006781',
			 'DU' => 'WBbt:0006781',
			 'dorsal uterine cell' => 'WBbt:0006782',
			 'du cell' => 'WBbt:0006782',
			 'blast cell'	 => 'WBbt:0006783',
			 'uterine toroidal epithelial cell'	 => 'WBbt:0006784',
			 'ut cell' => 'WBbt:0006784',
			 'uterine toroidal epithelial cell'	 => 'WBbt:0006784',
			 'uterine toroid' => 'WBbt:0006784',
			 'uterine seam cell'	 => 'WBbt:0006789',
			 'utse' => 'WBbt:0006789',
			 'uterine-vulval cell'	 => 'WBbt:0006790',
			 'uv cell' => 'WBbt:0006790',
			 'male gonad'	 => 'WBbt:0006794',
			 'proctodeum'	 => 'WBbt:0006795',
			 'germ cell'	 => 'WBbt:0006796',
			 'oocyte'	 => 'WBbt:0006797',
			 'ovum' => 'WBbt:0006797',
			 'sperm'	 => 'WBbt:0006798',
			 'spermatozoa' => 'WBbt:0006798',
			 'sperm'	 => 'WBbt:0006798',
			 'spermatozoon' => 'WBbt:0006798',
			 'spermatocyte'	 => 'WBbt:0006799',
			 'spermatid'	 => 'WBbt:0006800',
			 'outer labial neuron'	 => 'WBbt:0006801',
			 'OL neuron' => 'WBbt:0006801',
			 'lumbar neuron'	 => 'WBbt:0006802',
			 'lumbar ganglion neuron' => 'WBbt:0006802',
			 'body wall muscle cell'	 => 'WBbt:0006804',
			 'body muscle cell' => 'WBbt:0006804',
			 'body wall muscle cell'	 => 'WBbt:0006804',
			 'body wall muscle' => 'WBbt:0006804',
			 'body wall muscle cell'	 => 'WBbt:0006804',
			 'bodywall muscle' => 'WBbt:0006804',
			 'body wall muscle cell'	 => 'WBbt:0006804',
			 'adult seam cell left'	 => 'WBbt:0006809',
			 'adult seam cell right'	 => 'WBbt:0006810',
			 'ADA'	 => 'WBbt:0006811',
			 'AIA'	 => 'WBbt:0006812',
			 'AIB'	 => 'WBbt:0006813',
			 'AIM'	 => 'WBbt:0006814',
			 'AIN'	 => 'WBbt:0006815',
			 'ciliated neuron'	 => 'WBbt:0006816',
			 'AUA'	 => 'WBbt:0006817',
			 'AVD'	 => 'WBbt:0006818',
			 'AVE'	 => 'WBbt:0006819',
			 'AVH'	 => 'WBbt:0006821',
			 'AVJ'	 => 'WBbt:0006822',
			 'AVK'	 => 'WBbt:0006823',
			 'BAG'	 => 'WBbt:0006825',
			 'interneuron motoneuron' => 'WBbt:0006825',
			 'BDU'	 => 'WBbt:0006826',
			 'CAN'	 => 'WBbt:0006827',
			 'FLP'	 => 'WBbt:0006828',
			 'glutamatergic neuron'	 => 'WBbt:0006829',
			 'HSN'	 => 'WBbt:0006830',
			 'hermaphrodite specific neuron' => 'WBbt:0006830',
			 'PVD'	 => 'WBbt:0006831',
			 'PVP'	 => 'WBbt:0006832',
			 'RIA'	 => 'WBbt:0006833',
			 'RIC'	 => 'WBbt:0006834',
			 'RIF'	 => 'WBbt:0006835',
			 'RIG'	 => 'WBbt:0006836',
			 'serotonergic neuron'	 => 'WBbt:0006837',
			 '5-HT neuron' => 'WBbt:0006837',
			 'RIP'	 => 'WBbt:0006838',
			 'RIV'	 => 'WBbt:0006839',
			 'cholinergic neuron'	 => 'WBbt:0006840',
			 'ACh neuron' => 'WBbt:0006840',
			 'RMF'	 => 'WBbt:0006841',
			 'RMG'	 => 'WBbt:0006842',
			 'RMH'	 => 'WBbt:0006843',
			 'SDQ'	 => 'WBbt:0006844',
			 'URB'	 => 'WBbt:0006845',
			 'URX'	 => 'WBbt:0006846',
			 'OLL'	 => 'WBbt:0006847',
			 'OLQ'	 => 'WBbt:0006848',
			 'germline precursor cell'	 => 'WBbt:0006849',
			 'PGC' => 'WBbt:0006849',
			 'germline precursor cell'	 => 'WBbt:0006849',
			 'primordial germ cell' => 'WBbt:0006849',
			 'excretory secretory system'	 => 'WBbt:0006850',
			 'ALN'	 => 'WBbt:0006851',
			 'CA neuron'	 => 'WBbt:0006852',
			 'PLN'	 => 'WBbt:0006855',
			 'R5A'	 => 'WBbt:0006856',
			 'R7A'	 => 'WBbt:0006857',
			 'R9A'	 => 'WBbt:0006858',
			 'PCB'	 => 'WBbt:0006859',
			 'PCC'	 => 'WBbt:0006860',
			 'SPC'	 => 'WBbt:0006861',
			 'SPV'	 => 'WBbt:0006862',
			 'hermaphrodite distal tip cell'	 => 'WBbt:0006863',
			 'male distal tip cell' => 'WBbt:0006864',
			 'DTC' => 'WBbt:0006865',
			 'distal tip cell' => 'WBbt:0006865',
			 'anterior spermathecal-uterine valve cell'	 => 'WBbt:0006866',
			 'posterior spermathecal-uterine valve cell'	 => 'WBbt:0006867',
			 'anterior spermathecal-uterine core cell' => 'WBbt:0006868',
			 'posterior spermathecal-uterine core cell'	 => 'WBbt:0006869',
			 'seminal vesicle' => 'WBbt:0006870',
			 'anterior gonadal sheath cell'	 => 'WBbt:0006871',
			 'posterior gonadal sheath cell'	 => 'WBbt:0006872',
			 'EMS'	 => 'WBbt:0006876',
			 'developing hermaphrodite gonad'	 => 'WBbt:0006903',
			 'hook hypodermis' => 'WBbt:0006906',
			 'gubernacular muscle' => 'WBbt:0006908',
			 'longitudinal male muscle' => 'WBbt:0006909',
			 'ventral longitudinal muscle' => 'WBbt:0006909',
			 'oblique male muscle' => 'WBbt:0006910',
			 'spicule muscle' => 'WBbt:0006911',
			 'male tail seam cell'	 => 'WBbt:0006912',
			 'male posterior seam' => 'WBbt:0006912',
			 'male tail seam cell' => 'WBbt:0006912',
			 'hermaphrodite alae seam cell' => 'WBbt:0006913',
			 'male alae seam cell' => 'WBbt:0006914',
			 'uterine muscle, type 1' => 'WBbt:0006915',
			 'uterine muscle, type 2' => 'WBbt:0006916',
			 'vulval muscle, type 1' => 'WBbt:0006917',
			 'vulval muscle, type 2' => 'WBbt:0006918',
			 'anal region'	 => 'WBbt:0006919',
			 'cephalic sensillum'	 => 'WBbt:0006920',
			 'cephalic sensillum dorsal left'	 => 'WBbt:0006921',
			 'cephalic sensillum dorsal right'	 => 'WBbt:0006922',
			 'cephalic sensillum ventral right'	 => 'WBbt:0006923',
			 'cephalic sensillum ventral left'	 => 'WBbt:0006924',
			 'apoptotic cell'	 => 'WBbt:0006925',
			 'deirid sensillum'	 => 'WBbt:0006926',
			 'anterior deirid sensillum'	 => 'WBbt:0006927',
			 'deirid neuron'	 => 'WBbt:0006928',
			 'sensillum'	 => 'WBbt:0006929',
			 'hook sensillum'	 => 'WBbt:0006930',
			 'ray precursor cell'	 => 'WBbt:0006931',
			 'R1'	 => 'WBbt:0006932',
			 'R2'	 => 'WBbt:0006933',
			 'R3'	 => 'WBbt:0006934',
			 'R4'	 => 'WBbt:0006935',
			 'R5'	 => 'WBbt:0006936',
			 'R6'	 => 'WBbt:0006937',
			 'R7'	 => 'WBbt:0006938',
			 'R8'	 => 'WBbt:0006939',
			 'R9'	 => 'WBbt:0006940',
			 'ray'	 => 'WBbt:0006941',
			 'male sensory ray' => 'WBbt:0006941',
			 'ray 1'	 => 'WBbt:0006942',
			 'ray 1 right'	 => 'WBbt:0006943',
			 'ray 1 left'	 => 'WBbt:0006944',
			 'ray 2'	 => 'WBbt:0006945',
			 'ray 2 left'	 => 'WBbt:0006946',
			 'ray 2 right'	 => 'WBbt:0006947',
			 'ray 3'	 => 'WBbt:0006948',
			 'ray 4'	 => 'WBbt:0006949',
			 'ray 5'	 => 'WBbt:0006950',
			 'ray 6'	 => 'WBbt:0006951',
			 'ray 7'	 => 'WBbt:0006952',
			 'ray 8'	 => 'WBbt:0006953',
			 'ray 9'	 => 'WBbt:0006954',
			 'ray 3 left'	 => 'WBbt:0006955',
			 'ray 3 right'	 => 'WBbt:0006956',
			 'ray 4 left'	 => 'WBbt:0006957',
			 'ray 4 right'	 => 'WBbt:0006958',
			 'ray 5 left'	 => 'WBbt:0006959',
			 'ray 5 right'	 => 'WBbt:0006960',
			 'ray 6 left'	 => 'WBbt:0006961',
			 'ray 6 right'	 => 'WBbt:0006962',
			 'ray 7 left'	 => 'WBbt:0006963',
			 'ray 7 right'	 => 'WBbt:0006964',
			 'ray 8 left'	 => 'WBbt:0006965',
			 'ray 8 right'	 => 'WBbt:0006966',
			 'ray 9 left'	 => 'WBbt:0006967',
			 'ray 9 right'	 => 'WBbt:0006968',
			 'postdeirid sensillum left'	 => 'WBbt:0006969',
			 'postdeirid sensillum right'	 => 'WBbt:0006970',
			 'ray neuron'	 => 'WBbt:0006971',
			 'ray neuron type A'	 => 'WBbt:0006972',
			 'RnA' => 'WBbt:0006972',
			 'ray neuron type B'	 => 'WBbt:0006973',
			 'RnB' => 'WBbt:0006973',
			 'nerve ring neuron'	 => 'WBbt:0006974',
			 'PVQ'	 => 'WBbt:0006976',
			 'tail ganglion'	 => 'WBbt:0006977',
			 'tail hypodermis'	 => 'WBbt:0006978',
			 'tail spike'	 => 'WBbt:0006979',
			 'embryonic cell'	 => 'WBbt:0007028',
			 'post-embryonic cell'	 => 'WBbt:0007030',
			 'engulfing cell'	 => 'WBbt:0007803',
			 'command interneuron'	 => 'WBbt:0007804',
			 'command neuron' => 'WBbt:0007804',
			 'SPD'	 => 'WBbt:0007805',
			 'PCA'	 => 'WBbt:0007806',
			 'PHA'	 => 'WBbt:0007807',
			 'PHB'	 => 'WBbt:0007808',
			 'vulval precursor cell'	 => 'WBbt:0007809',
			 'VPC' => 'WBbt:0007809',
			 'vulval precursor cell'	 => 'WBbt:0007809',
			 'vulval precursor' => 'WBbt:0007809',
			 'body muscle cell'	 => 'WBbt:0007810',
			 'uterine pi cell'	 => 'WBbt:0007813',
			 'uterine rho cell'	 => 'WBbt:0007814',
			 'hermaphrodite somatic gonadal cell'	 => 'WBbt:0007815',
			 'spermathecal cell' => 'WBbt:0007816',
			 'SM3L'	 => 'WBbt:0007817',
			 'left SM3' => 'WBbt:0007817',
			 'SM3R'	 => 'WBbt:0007818',
			 'right SM3' => 'WBbt:0007818',
			 'B_gamma'	 => 'WBbt:0007830',
			 'vulval cell'	 => 'WBbt:0007831',
			 'cuticular ala'	 => 'WBbt:0007832',
			 'alae' => 'WBbt:0007832',
			 'organism'	 => 'WBbt:0007833',
			 'DX neuron'	 => 'WBbt:0007834',
			 'EF neuron'	 => 'WBbt:0007835',
			 'PVU'	 => 'WBbt:0007836',
			 'PVPL' => 'WBbt:0007836',
			 'PVS'	 => 'WBbt:0007837',
			 'PVPR' => 'WBbt:0007837',
			 'EF1'	 => 'WBbt:0007838',
			 'EF2'	 => 'WBbt:0007839',
			 'EF3'	 => 'WBbt:0007840',
			 'EF4'	 => 'WBbt:0007841',
			 'DX1'	 => 'WBbt:0007842',
			 'DX2'	 => 'WBbt:0007843',
			 'DX3'	 => 'WBbt:0007844',
			 'DX4'	 => 'WBbt:0007845',
			 'hypodermal cell'	 => 'WBbt:0007846',
			 'hook hypodermal cell'	 => 'WBbt:0007847',
			 'pharyngeal-intestinal valve cell'	 => 'WBbt:0007848',
			 'SPsoL'	 => 'WBbt:0007851',
			 'SPsoR'	 => 'WBbt:0007852',
			 'somatic gonad precursor'	 => 'WBbt:0007854',
			 'SGP' => 'WBbt:0007854',
			 'somatic gonad precursor'	 => 'WBbt:0007854',
			 'somatic gonadal precursor cell' => 'WBbt:0007854',
			 'XXX cell'	 => 'WBbt:0007855',
			 'LUA'	 => 'WBbt:0008031',
			 'R1L_hyp' => 'WBbt:0008043',
			 'R1R_hyp' => 'WBbt:0008044',
			 'R2L_hyp' => 'WBbt:0008045',
			 'R2R_hyp' => 'WBbt:0008046',
			 'R3L_hyp' => 'WBbt:0008047',
			 'R3R_hyp' => 'WBbt:0008048',
			 'R4L_hyp' => 'WBbt:0008049',
			 'R4R_hyp' => 'WBbt:0008050',
			 'R5L_hyp' => 'WBbt:0008051',
			 'R5R_hyp' => 'WBbt:0008052',
			 'R6L_hyp' => 'WBbt:0008053',
			 'R6R_hyp' => 'WBbt:0008054',
			 'R7L_hyp' => 'WBbt:0008055',
			 'R7R_hyp' => 'WBbt:0008056',
			 'R8L_hyp' => 'WBbt:0008057',
			 'R8R_hyp' => 'WBbt:0008058',
			 'R9L_hyp' => 'WBbt:0008059',
			 'R9R_hyp' => 'WBbt:0008060',
			 'B cell male'	 => 'WBbt:0008064',
			 'hyp7 precursor cell'	 => 'WBbt:0008069',
			 'hyp7 syncytium male'	 => 'WBbt:0008070',
			 'hyp7 syncytium hermaphrodite'	 => 'WBbt:0008071',
			 'spike' => 'WBbt:0008072',
			 'tail spike cell' => 'WBbt:0008072',
			 'syncytium precursor'	 => 'WBbt:0008073',
			 'syncytium'	 => 'WBbt:0008074',
			 'syncitium' => 'WBbt:0008074',
			 'syncytium'	 => 'WBbt:0008074',
			 'syncycia' => 'WBbt:0008074',
			 'syncytium'	 => 'WBbt:0008074',
			 'syncytial cell' => 'WBbt:0008074',
			 'ventral cord blast cell'	 => 'WBbt:0008115',
			 'P cell' => 'WBbt:0008115',
			 'adult seam cell left male'	 => 'WBbt:0008136',
			 'adult seam cell left hermaphrodite'	 => 'WBbt:0008137',
			 'adult seam cell right male'	 => 'WBbt:0008138',
			 'adult seam cell right hermaphrodite'	 => 'WBbt:0008139',
			 'se seam cell' => 'WBbt:0008156',
			 'rectal epithelial cell'	 => 'WBbt:0008182',
			 'rectal epithelium' => 'WBbt:0008182',
			 'male gonadal cell'	 => 'WBbt:0008209',
			 'male somatic gonadal cell' => 'WBbt:0008209',
			 'anterior spermathecal-uterine valve cell 5L'	 => 'WBbt:0008211',
			 'anterior spermathecal-uterine valve cell 5R'	 => 'WBbt:0008212',
			 'posterior spermathecal-uterine valve cell 5L'	 => 'WBbt:0008214',
			 'posterior spermathecal-uterine valve cell 5R'	 => 'WBbt:0008215',
			 'spermathecal-uterine valve cell' => 'WBbt:0008217',
			 'spermathecal-uterine core cell' => 'WBbt:0008218',
			 'uterine seam cell 5L' => 'WBbt:0008219',
			 'uterine seam cell 5R' => 'WBbt:0008220',
			 'linker cell Z1'	 => 'WBbt:0008239',
			 'linker cell Z4'	 => 'WBbt:0008240',
			 'F cell hermaphrodite'	 => 'WBbt:0008243',
			 'F cell male'	 => 'WBbt:0008244',
			 'U cell hermaphrodite'	 => 'WBbt:0008245',
			 'U cell male'	 => 'WBbt:0008246',
			 'Y cell hermaphrodite'	 => 'WBbt:0008247',
			 'Y cell male'	 => 'WBbt:0008248',
			 'VLL muscle'	 => 'WBbt:0008255',
			 'VML muscle'	 => 'WBbt:0008256',
			 'VL muscle'	 => 'WBbt:0008257',
			 'VR muscle'	 => 'WBbt:0008259',
			 'DL muscle'	 => 'WBbt:0008260',
			 'DR muscle'	 => 'WBbt:0008261',
			 'body wall muscle quadrant'	 => 'WBbt:0008262',
			 'DLR muscle'	 => 'WBbt:0008263',
			 'DMR muscle'	 => 'WBbt:0008264',
			 'VLR'	 => 'WBbt:0008265',
			 'VMR'	 => 'WBbt:0008266',
			 'DLL muscle'	 => 'WBbt:0008267',
			 'DML muscle'	 => 'WBbt:0008268',
			 'pm1 precursor'	 => 'WBbt:0008269',
			 'pm2 precursor'	 => 'WBbt:0008270',
			 'pm3 precursor'	 => 'WBbt:0008271',
			 'pm4 precursor'	 => 'WBbt:0008272',
			 'pm5 precursor'	 => 'WBbt:0008273',
			 'proctodeal cell' => 'WBbt:0008274',
			 'PVN'	 => 'WBbt:0008275',
			 'ray structural cell' => 'WBbt:0008276',
			 'vulval syncytia precursor'	 => 'WBbt:0008289',
			 'uterine pi cell (5L)'	 => 'WBbt:0008308',
			 'uterine pi cell (5R)'	 => 'WBbt:0008309',
			 'uterine rho cell (5L)'	 => 'WBbt:0008310',
			 'uterine rho cell (5R)'	 => 'WBbt:0008311',
			 'mc3 precursor'	 => 'WBbt:0008349',
			 'g1ARP'	 => 'WBbt:0008350',
			 'g1AR' => 'WBbt:0008350',
			 'g1ARP'	 => 'WBbt:0008350',
			 'g1P' => 'WBbt:0008350',
			 'set left'	 => 'WBbt:0008351',
			 'set right'	 => 'WBbt:0008352',
			 'set precursor'	 => 'WBbt:0008353',
			 'male tail tip cell'	 => 'WBbt:0008354',
			 'vas deferens precursor'	 => 'WBbt:0008357',
			 'vas deferens precursor 1L'	 => 'WBbt:0008358',
			 'vas deferens precursor 4L'	 => 'WBbt:0008359',
			 'gonadal primordium'	 => 'WBbt:0008366',
			 'gonadal primordium male'	 => 'WBbt:0008367',
			 'gonadal primordium hermaphrodite'	 => 'WBbt:0008368',
			 'sex myoblast'	 => 'WBbt:0008373',
			 'stoma'	 => 'WBbt:0008374',
			 'arcade ring'	 => 'WBbt:0008375',
			 'anterior arcade ring'	 => 'WBbt:0008376',
			 'posterior arcade ring'	 => 'WBbt:0008377',
			 'somatic cell' => 'WBbt:0008378',
			 'amphid socket cell'	 => 'WBbt:0008379',
			 'AMso' => 'WBbt:0008379',
			 'diagonal muscle right'	 => 'WBbt:0008380',
			 'dglR' => 'WBbt:0008380',
			 'diagonal muscle left'	 => 'WBbt:0008381',
			 'dglL' => 'WBbt:0008381',
			 'outer logitudinal muscle'	 => 'WBbt:0008382',
			 'inner longitudinal muscle'	 => 'WBbt:0008383',
			 'caudal longitudinal muscle'	 => 'WBbt:0008384',
			 'cdl' => 'WBbt:0008384',
			 'anal depressor muscle male'	 => 'WBbt:0008385',
			 'anal depressor muscle hermaphrodite'	 => 'WBbt:0008386',
			 'anal sphincter muscle male'	 => 'WBbt:0008387',
			 'anal sphincter muscle hermaphrodite'	 => 'WBbt:0008388',
			 'dorsal body wall muscle'	 => 'WBbt:0008389',
			 'ventral body wall muscle'	 => 'WBbt:0008390',
			 'R1A'	 => 'WBbt:0008391',
			 'R2A'	 => 'WBbt:0008392',
			 'R3A'	 => 'WBbt:0008393',
			 'R4A'	 => 'WBbt:0008394',
			 'R6A'	 => 'WBbt:0008395',
			 'R8A'	 => 'WBbt:0008396',
			 'R1B'	 => 'WBbt:0008397',
			 'R2B'	 => 'WBbt:0008398',
			 'R3B'	 => 'WBbt:0008399',
			 'R4B'	 => 'WBbt:0008400',
			 'R5B'	 => 'WBbt:0008401',
			 'R6B'	 => 'WBbt:0008402',
			 'R7B'	 => 'WBbt:0008403',
			 'R8B'	 => 'WBbt:0008404',
			 'R9B'	 => 'WBbt:0008405',
			 'cephalic sheath cell'	 => 'WBbt:0008406',
			 'grinder'	 => 'WBbt:0008408',
			 'tail precursor cell'	 => 'WBbt:0008409',
			 'T cell' => 'WBbt:0008409',
			 'phasmid socket cell'	 => 'WBbt:0008410',
			 'ADE sheath cell'	 => 'WBbt:0008411',
			 'IL sheath cell'	 => 'WBbt:0008412',
			 'OL sheath cell'	 => 'WBbt:0008413',
			 'PC sheath cell'	 => 'WBbt:0008414',
			 'PDE sheath cell'	 => 'WBbt:0008415',
			 'ADE socket cell'	 => 'WBbt:0008416',
			 'CEP socket cell'	 => 'WBbt:0008417',
			 'IL socket cell'	 => 'WBbt:0008418',
			 'OL socket cell'	 => 'WBbt:0008419',
			 'PC socket cell'	 => 'WBbt:0008420',
			 'PDE socket cell'	 => 'WBbt:0008421',
			 'sex organ'	 => 'WBbt:0008422',
			 'genitalia' => 'WBbt:0008422',
			 'male genital'	 => 'WBbt:0008423',
			 'gubernaculum'	 => 'WBbt:0008424',
			 'hook'	 => 'WBbt:0008425',
			 'nose'	 => 'WBbt:0008426',
			 'blastocoel'	 => 'WBbt:0008427',
			 'blastopore'	 => 'WBbt:0008428',
			 'gastrulation pore' => 'WBbt:0008428',
			 'polar body'	 => 'WBbt:0008429',
			 'carbon dioxide sensory neuron'	 => 'WBbt:0008430',
			 'CO2 sensory neuron' => 'WBbt:0008430',
			 'mechanosensory neuron'	 => 'WBbt:0008431',
			 'oxygene sensory neuron'	 => 'WBbt:0008432',
			 'osmosensory neuron'	 => 'WBbt:0008433',
			 'nociceptor neuron'	 => 'WBbt:0008434',
			 'uterine epithelial cell'	 => 'WBbt:0008435',
			 'gon herm vut (5L)'	 => 'WBbt:0008436',
			 'gon herm vut (5R)'	 => 'WBbt:0008437',
			 'PVW'	 => 'WBbt:0008438',
			 'vas deferens valve region'	 => 'WBbt:0008440',
			 'vas deferens valve cell'	 => 'WBbt:0008441',
			 'vas deferens cuboidal cell region'	 => 'WBbt:0008442',
			 'vas deferens elongated cell region'	 => 'WBbt:0008443',
			 'vas deferens cuboidal cell'	 => 'WBbt:0008444',
			 'vas deferens elongated cell'	 => 'WBbt:0008445',
			 'PHC'	 => 'WBbt:0008446',
			 'somatic blastomere'	 => 'WBbt:0008447',
			 'germline blastomere'	 => 'WBbt:0008448',
			 'B_delta'	 => 'WBbt:0008451',
			 'B_alpha'	 => 'WBbt:0008506',
			 'B_beta'	 => 'WBbt:0008507',
			 'egg-laying apparatus'	 => 'WBbt:0008587',
			 'primary spermatocyte'	 => 'WBbt:0008591',
			 'secondary spermatocyte'	 => 'WBbt:0008592',
			 'neuroblast'	 => 'WBbt:0008594',
			 'H cell'	 => 'WBbt:0008596',
			 'V cell'	 => 'WBbt:0008597',
			 'Q cell'	 => 'WBbt:0008598',
			 'G cell'	 => 'WBbt:0008599',
			 'enteric muscle'	 => 'WBbt:0008600',
			 'male sex myoblast'	 => 'WBbt:0008604',
			 'hermaphrodite sex myoblast'	 => 'WBbt:0008605',
			 'C cell'        => 'WBbt:0003810',
			 'E cell'        => 'WBbt:0004804',
			 'MS cell'       => 'WBbt:0004458',
			);

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

  # assume we are using the same tissue as in the last experiment
  if (!defined $candidate_tissue && defined $tissue_session{previous}) {
    $candidate_tissue = $tissue_session{previous};
    $why = 'repeat of last expt';
  }

  # else if this is the first time in, use the default
  if (!defined $candidate_tissue) {
    $candidate_tissue = $default;
    $why = 'default';
  }

  # confirm and update status
  do {
    print "Tissue [$candidate_tissue] ($why) > ";
    my $input =  <STDIN>;
    chomp ($input);
    if ($input eq '') {
      $tissue = $tissue_ontology{$candidate_tissue};
      $tissue_name = $candidate_tissue;
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
  my $life_stage; 
  my $life_stage_name; 
  my $candidate_life_stage = undef;
  my $line;
  my $why = '';
  
  my $default = 'all stages'; # WBls:0000002
  
  my %life_stage_ontology;
  
  if ($species eq 'brugia') {
    %life_stage_ontology = (
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
			   );
  } elsif ($species eq 'sratti') {
    %life_stage_ontology = (
			    'parasitic females'	 =>	'WBls:0000678',
			    'free-living females'	 =>	'WBls:0000677',
			    'free living adults'	 =>	'WBls:0000682',
			    'infective larvae (iL3)'	=>	'WBls:0000680',
			   );
    
  } elsif ($species eq 'elegans' || $species eq 'brenneri' || $species eq 'briggsae' || $species eq 'remanei' || $species eq 'japonica') {
    
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

  } else {
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

  # assume we are using the same life_stage as in the last experiment
  if (!defined $candidate_life_stage && defined $life_stage_session{previous}) {
    $candidate_life_stage = $life_stage_session{previous};
    $why = 'repeat of last expt';
  }

  # else if this is the first time in, use the default
  if (!defined $candidate_life_stage) {
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
