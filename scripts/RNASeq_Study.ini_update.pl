#!/usr/bin/env perl
#
# RNASeq_Study.ini_update.pl
#
# This searches the ENA SRA for Studies from the
# nematode/platyhelminth species.
#
# For each species it updates the .ini files with details of new
# Studies and then makes an INI file of the details of Experiments in
# that Study.

use strict;
use Carp;
use Getopt::Long;
use Config::IniFiles;

my ($help, $dir, $stats_parasite, $stats_wormbase);
GetOptions ("help"       => \$help,
            "dir:s"      => \$dir, # the ouput directory
	    "stats_parasite" => \$stats_parasite,
	    "stats_wormbase" => \$stats_wormbase,
);


my @wormbase_species = qw(
		Caenorhabditis_elegans	
		Caenorhabditis_brenneri
		Caenorhabditis_briggsae
		Caenorhabditis_japonica
		Caenorhabditis_remanei
		Brugia_malayi
		Onchocerca_volvulus
		Pristionchus_pacificus
);

# this is an estimate - these species should be updated as and when required
my @parasite_species = qw (
			   Brugia_malayi
			   Echinococcus_granulosus
			   Echinococcus_multilocularis
			   Hymenolepis_microstoma
			   Onchocerca_volvulus
			   Strongyloides_ratti
			   Strongyloides_stercoralis
);

my $report_stats = 0;
my @check_species;
if ($stats_parasite) {$report_stats = 1; @check_species = @parasite_species}
if ($stats_wormbase) {$report_stats = 1; @check_species = @wormbase_species}

if (!defined $dir) {die "-dir output directory is not defined\n";}
mkdir $dir, 0755;
mkdir "$dir/nematoda", 0755;
mkdir "$dir/platyhelminthes", 0755;

my %updated_files;

# Retrieve all studies from ENA where the library strategy is RNASeq
# and the taxon is Nematoda or Platyhelminth. 
# %20 is space
# %22 is "
# %28 is (
# %29 is )

# To get all nematode and platyhelminth Studies:
# http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22tax_tree(6231)%20AND%20library_strategy=%22RNA-Seq%22%22&domain=read
# and
# http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22tax_tree(6157)%20AND%20library_strategy=%22RNA-Seq%22%22&domain=read

# tabulated data is returned by:
# http://www.ebi.ac.uk/ena/data/warehouse/search?query=<query string>&result=<result>&fields=<fields>&display=report[&sortfields=<sortfields>][&download=txt][Pagination options]

# So for all nematode species:
# http://www.ebi.ac.uk/ena/data/warehouse/search?query=tax_tree(6231)%20AND%20library_strategy=%22RNA-Seq%22&result=read_run&fields=study_accession,secondary_study_accession,study_title,experiment_accession,library_source,library_selection,center_name,scientific_name,tax_id&display=report

my %species_files_seen; # keep track on which species have been pulled out of the ENA
my %studies_not_seen; # note which studies have not been pulled out of the ENA
my %experiments_not_seen; # note which experiments have not been pulled out of the ENA
my %novel_experiments; # note new experiments in existing Studies

my %studies_stats;
my %experiments_stats;
my %selection_stats;
my %layout_stats;
my %size_stats;
my %reads_stats;


my %phyla=('6231' => 'nematoda', '6157' => 'platyhelminthes');
foreach my $phylum (keys %phyla) { # nematode and platyhelminth
  my %data = &get_phylup_data($phylum);
  my $phylum_name = $phyla{$phylum};

  # Create an ini file for all species with at least one study
  foreach my $species (keys %data) {
    my $species_name = $species;
    $species_name =~ s/ /_/g;
    my $species_file = $species_name . ".ini";
    my $path = "$dir/${phylum_name}/$species_file";
    $species_files_seen{"${phylum_name}/$species_file"} = 1;
    if (!-e $path) {system("touch $path")}
    my $ini = Config::IniFiles->new( -file => $path, -allowempty => 1);
    my @existing_studies = $ini->Sections;

    foreach my $secondary_study_accession (keys $data{$species}) {

      # Each ini file should contain the following data from ENA, in per-study stanzas:
      #
      #	- study title
      #	- study accession
      #	- secondary study accession
      #	- library source (if multiple for study list all)
      #	- library selection (if multiple for study list all)
      #	- scientific name
      #	- tax id
      #	- centre name
      #
      # In addition, there should be the following tags in the ini file:
      #
      #	-remark
      #	-use
      #	-status
      
      my @experiments = @{$data{$species}{$secondary_study_accession}{experiment}};

      # If the *only* library source is GENOMIC, don’t create a stanza.
      my $not_genomic = 0;
      foreach my $expt (@experiments) {
	my ($experiment_accession, $library_source, $library_selection, $secondary_sample_accession) = @{$expt};
	if ($library_source ne 'GENOMIC') {$not_genomic = 1}
      }
      
      if ($not_genomic) {
	if (grep /$species_name/, @check_species) {$studies_stats{$species_name}++}
	# only write these if this is a new study, otherwise they would be overwritten
	if (!grep /^$secondary_study_accession$/, @existing_studies) {
	  print "New Study: $species\t$secondary_study_accession\n";
	  $updated_files{"${phylum_name}/$species_file"} = 1;
	  $ini->AddSection($secondary_study_accession);
	  $ini->newval($secondary_study_accession, 'use', '');
	  $ini->newval($secondary_study_accession, 'remark', '');
	  $ini->newval($secondary_study_accession, 'pubmed', '');
	  $ini->newval($secondary_study_accession, 'status', '');
	  $ini->newval($secondary_study_accession, 'ArrayExpress_ID', '');
	  
	  $ini->newval($secondary_study_accession, 'scientific_name', $species);
	  $ini->newval($secondary_study_accession, 'secondary_study_accession', $secondary_study_accession);
	  foreach my $param (qw(study_accession study_title center_name tax_id)) {
	    my $value = $data{$species}{$secondary_study_accession}{$param};
	    $ini->newval($secondary_study_accession, $param, $value);
	  }
	  $ini->newval($secondary_study_accession, "Colour", '');

	  # we can't get this field from the 'read_run' results section of the ENA warehouse, it has to be a separate query
	  my $study_description = &get_study_description($secondary_study_accession);
	  $ini->newval($secondary_study_accession, 'study_description', $study_description);
	}
	
	my %samples = ();
	my %new_samples = ();
	foreach my $expt (@experiments) {
	  my ($experiment_accession, $library_source, $library_selection, $secondary_sample_accession, $library_layout) = @{$expt};
	  
	  if (grep /$species_name/, @check_species) {$experiments_stats{$species_name}++}
	  if (grep /$species_name/, @check_species) {$selection_stats{$species_name}{$library_selection}++}
	  if (grep (/$species_name/, @check_species) && ($library_selection eq 'cDNA')) {$layout_stats{$species_name}{$library_layout}++}

	  # see if this is a new Experiment in an existing study
	  if (! $ini->exists($secondary_study_accession, "library_selection_$experiment_accession") && grep /^$secondary_study_accession$/, @existing_studies) {
	    push @{$novel_experiments{"${phylum_name}/$species_file"}}, $experiment_accession;
	  }
	  
	  if (! $ini->exists($secondary_study_accession, "library_selection_$experiment_accession")) {
	    $new_samples{$secondary_sample_accession} = 1; # note this experiment is being written and so we will need the sample tags as well
	  }
	  
	  $ini->newval($secondary_study_accession, "library_selection_$experiment_accession", "$library_selection");      
	  $ini->newval($secondary_study_accession, "library_sample_$experiment_accession", "$secondary_sample_accession");      	  
	  $samples{$secondary_sample_accession} = 1;
	}
	
	foreach my $secondary_sample_accession (keys %samples) {
	  if (exists $new_samples{$secondary_sample_accession}) {
	    $ini->newval($secondary_study_accession, "sample_longLabel_$secondary_sample_accession", '');
	    $ini->newval($secondary_study_accession, "sample_shortLabel_$secondary_sample_accession", '');
	    $ini->newval($secondary_study_accession, "sample_ChEBI_ID_$secondary_sample_accession", '');
	    $ini->newval($secondary_study_accession, "sample_WormBaseLifeStage_$secondary_sample_accession", '');
	  }
	}
	
      }
    }
    $ini->RewriteConfig;

    # check to see if all the Studies pulled out of the ENA match the Studies in the .ini file for this species
    foreach my $study (@existing_studies) {
      if (! exists $data{$species}{$study}) {
	push @{$studies_not_seen{"${phylum_name}/$species_file"}}, $study; 
      }
    }
  }
}

print "Finished updating\n";
print "You may now do a git commit and git push on the following changed files:\n\n";
foreach my $file (sort {$a cmp $b} keys %updated_files) {
  print "\t$file\n";
}
print "\n";

# check that all the species files are still in use
foreach my $phylum (keys %phyla) { # nematode and platyhelminth
  my $phylum_name = $phyla{$phylum};
  opendir(my $dh, "$dir/$phylum_name") || die "can't opendir $dir/$phylum_name: $!";
  while(my $file = readdir $dh) {
    if ($file eq '.' || $file eq '..') {next}
    if ($file =~ /~$/) {next} #  ignore emacs backup files
    if (! exists $species_files_seen{"$phylum_name/$file"}) {
      print "WARNING: File '$phylum_name/$file' no longer holds data found in the ENA (has the species name changed?)\n";
    }
  }
  closedir $dh;
}

# check that all studies are still in use
foreach my $sp (keys %studies_not_seen) {
  foreach my $st (@{$studies_not_seen{$sp}}) {
    print "WARNING: Study '$st' in '$sp' is no longer found in the ENA\n";
  }
}

# report novel experiments in existing Studies
foreach my $sp (keys %novel_experiments) {
  foreach my $ex (@{$novel_experiments{$sp}}) {
    print "WARNING: Experiment '$ex' in '$sp' has been recently added to an existing Study\n";
  }
}

# report statistics
if ($report_stats) {
  print "Species\t\tStudies\tExperiments\n";
  foreach my $sp (@check_species) {
    print "$sp\t$studies_stats{$sp}\t$experiments_stats{$sp}";
    foreach my $selection (keys %{$selection_stats{$sp}}) {
      print "\t$selection: $selection_stats{$sp}{$selection}";
    }
     print "\tLayout_PAIRED: $layout_stats{$sp}{PAIRED} Layout_SINGLE: $layout_stats{$sp}{SINGLE}"; # this is only reporting values where library_selection is 'cDNA'
    print "\n";
  }


}




###########################################################

sub get_phylup_data {
  my ($phylum) = @_;
  
  my $query = 'http://www.ebi.ac.uk/ena/data/warehouse/search?query=tax_tree(PHYLUM)%20AND%20(library_strategy=%22RNA-Seq%22%20OR%20library_strategy=%22FL-cDNA%22%20OR%20library_strategy=%22EST%22)&result=read_run&fields=study_accession,secondary_study_accession,study_title,secondary_sample_accession,experiment_accession,library_source,library_selection,library_strategy,center_name,scientific_name,tax_id,library_layout&display=report';
  $query =~ s/PHYLUM/$phylum/;
  
  my %data;
  my $first=1;
  open (DATA, "wget -q -O - '$query' |") || die("Can't get information on SRA entries in get_phylum_data()\n");
#  print "$query\n";
  while (my $line = <DATA>) {
    if ($first) {$first=0; next;} # skip title line
    chomp $line;
    my (
	$run_accession,               # don't want to store this, we just get it because we queried against 'read_run'
	$study_accession,
	$secondary_study_accession,
	$study_title,
	$secondary_sample_accession,
	$experiment_accession,
	$library_source,
	$library_selection,
	$library_stategy,
	$center_name,
	$scientific_name,
	$tax_id,
	$library_layout,
       ) = split /\t/, $line;
#    if ($library_stategy ne 'RNA-Seq' && $scientific_name !~ /^Caenorhabditis/) {next}; # The Parasite project is only interested in $library_stategy = 'RNA-Seq', but WormBase likes the others as well
    if ($library_stategy ne 'RNA-Seq') {next}; # The Parasite project is only interested in $library_stategy = 'RNA-Seq'

    $data{$scientific_name}{$secondary_study_accession}{study_accession} = $study_accession;
    $data{$scientific_name}{$secondary_study_accession}{study_title} = $study_title;
    $data{$scientific_name}{$secondary_study_accession}{center_name} = $center_name;
    $data{$scientific_name}{$secondary_study_accession}{tax_id} = $tax_id;
    push @{$data{$scientific_name}{$secondary_study_accession}{experiment}}, [$experiment_accession, $library_source, $library_selection, $secondary_sample_accession, $library_layout];
  }
  close(DATA);

  if (! scalar keys %data) {die("Can't get information on SRA entries in get_phylum_data()\n")}
  return %data;
}

###########################################################
# get the study_description from the 'study' results
# http://www.ebi.ac.uk/ena/data/warehouse/search?query=secondary_study_accession=SRP004901&result=study&fields=study_description&display=report

sub get_study_description {
  my ($secondary_study_accession) = @_;

  my $query = 'http://www.ebi.ac.uk/ena/data/warehouse/search?query=secondary_study_accession=SECONDARY_STUDY_ACCESSION&result=study&fields=study_description&display=report';
  $query =~ s/SECONDARY_STUDY_ACCESSION/$secondary_study_accession/;
  
  my $result;
  my $first=1;
  open (DATA, "wget -q -O - '$query' |") || die("Can't get information on SRA entries in get_study_description()\n");
#  print "$query\n";
  while (my $line = <DATA>) {
    if ($first) {$first=0; next;} # skip title line
    chomp $line;
    my (
	$study_accession,
	$study_description,
       ) = split /\t/, $line;
    $result  = $study_description;
  }
  close(DATA);

  return $result;
}
