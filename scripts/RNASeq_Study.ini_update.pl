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

my ($help, $dir);
GetOptions ("help"       => \$help,
            "dir:s"      => \$dir, # the ouput directory
);

if (!defined $dir) {$dir = glob("~/Parasite_config");}
mkdir $dir, 0755;

my @updated_files;

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


foreach my $phylum ('6231', '6157') { # nematode and platyhelminth
  my %data = &get_phylup_data($phylum);


###  # Exclude all C. elegans studies.
###  delete $data{'Caenorhabditis elegans'};

  # Create an ini file for all species with at least one study
  foreach my $species (keys %data) {
    my $species_file = $species;
    $species_file =~ s/ /_/g;
    $species_file .= ".ini";
    my $path = "$dir/$species_file";
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
	my ($experiment_accession, $library_source, $library_selection) = @{$expt};
	if ($library_source ne 'GENOMIC') {$not_genomic = 1}
      }

      if ($not_genomic) {
	if (!grep /^$secondary_study_accession$/, @existing_studies) {
	  print "New Study: $species\t$secondary_study_accession\n";
	  $ini->AddSection($secondary_study_accession);
	  $ini->newval($secondary_study_accession, 'remark', '');
	  $ini->newval($secondary_study_accession, 'use', '');
	  $ini->newval($secondary_study_accession, 'status', '');
	}

	foreach my $param (qw(study_accession study_title center_name tax_id)) {
	  my $value = $data{$species}{$secondary_study_accession}{$param};
	  $ini->newval($secondary_study_accession, $param, $value);
	}

	# we can't get this field from the 'read_run' results section of the ENA warehouse, it has to be a separate query
	my $study_description = &get_study_description($secondary_study_accession);
	$ini->newval($secondary_study_accession, 'study_description', $study_description);

	foreach my $expt (@experiments) {
	  my ($experiment_accession, $library_source, $library_selection) = @{$expt};
	  #	$ini->newval($secondary_study_accession, 'experiment_accession', "$experiment_accession");
	  $ini->newval($secondary_study_accession, "library_source_$experiment_accession", "$library_source");
	  $ini->newval($secondary_study_accession, "library_selection_$experiment_accession", "$library_selection");      
	}

	push @updated_files, $species_file;
      }
    }
    $ini->RewriteConfig;
  }
}

print "Finished updating\n";
print "You may now do a git commit and git push on the following changed file:\n\n";
foreach my $file (@updated_files) {
  print "\t$file\n";
}
print "\n";

###########################################################

sub get_phylup_data {
  my ($phylum) = @_;
  
  my $query = 'http://www.ebi.ac.uk/ena/data/warehouse/search?query=tax_tree(PHYLUM)%20AND%20library_strategy=%22RNA-Seq%22&result=read_run&fields=study_accession,secondary_study_accession,study_title,experiment_accession,library_source,library_selection,center_name,scientific_name,tax_id&display=report';
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
	$experiment_accession,
	$library_source,
	$library_selection,
	$center_name,
	$scientific_name,
	$tax_id,
       ) = split /\t/, $line;
    $data{$scientific_name}{$secondary_study_accession}{study_accession} = $study_accession;
    $data{$scientific_name}{$secondary_study_accession}{study_title} = $study_title;
    $data{$scientific_name}{$secondary_study_accession}{center_name} = $center_name;
    $data{$scientific_name}{$secondary_study_accession}{tax_id} = $tax_id;
    push @{$data{$scientific_name}{$secondary_study_accession}{experiment}}, [$experiment_accession, $library_source, $library_selection];
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
