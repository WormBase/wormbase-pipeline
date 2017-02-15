#!/software/bin/perl -w
#
# RNASeq.pm                       
# 
# by Gary Williams
#
# Methods for running the RNAseq pipeline and other useful things like searching the ENA warehouse
#
# Last updated by: $Author: gw3 $     
# Last updated on: $Date: 2015-06-09 14:32:44 $      

=pod

=head1 NAME

 RNASeq

=head1 SYNOPSIS

  use RNASeq;

  # initialise
  my $RNASeq = RNASeq->new($wormbase, $log, $new_genome, $check);

  # query the ENA warehouse data
  $data = $RNASeq->read_accession($accession);
  $data = $self->find_studies();

  # manipulate the hash ref of ENA data
  $data = $self->get_experiments_from_run_data($run_data_hash);

  # get data stored in the INI config files
  $ini = $self->read_all_studies_config();
  $ini = $self->read_experiments_from_study_config($study_accession);
  $data = $self->get_all_experiments();

  # write new data to the INI config files
  @new_study_accessions = $self->find_new_studies($study_ini, $studies);
  $ini = $self->add_new_studies($study_ini, $studies, @new_study_accessions);
  $self->add_new_experiments_in_study($study_accession, $study_accession)
  @studies_changed = $self->update_study_config_data();
  $self->add_param_to_experiment_config( $experiment_accession, $param, $value)

  # manipulate the INI config data
  $data = $self->convert_ini_to_hashref($ini_object);


=head1 DESCRIPTION

Methods for running the RNAseq pipeline and other useful things like
searching the ENA warehouse.

=head1 CONTACT

Gary gw3@ebi.ac.uk

=head1 APPENDIX

 The rest of the documentation details each of the object methods.
 Internal methods are preceded with a _

=cut




package RNASeq;

use strict;             
#use XML::LibXML;
use LWP;
use lib $ENV{'CVS_DIR'};
use Carp;
use Data::Dumper;
use Config::IniFiles;
use Coords_converter;
use File::Basename;

=head2 

    Title   :   new
    Usage   :   my $RNASeq = RNASeq->new($wormbase, $log);
    Function:   initialises the object
    Returns :   RNASeq object;
    Args    :   

=cut

sub new {
  my $class = shift;
  my $self = {};
  bless $self, $class;

  $self->{'wormbase'} = shift;    
  $self->{'log'} = shift;
  $self->{'new_genome'} = shift; # true if the libraries should be mapped to the genome in this session (e.g. if the assembly changed)
  $self->{'check'} = shift;      # true if existing GTF and cufflinks data should be left untouched (for checkpointing and restarting)

  # set up useful paths etc.
  $self->{'RNASeqBase'}      = "/nfs/nobackup/ensemblgenomes/wormbase/BUILD/RNASeq/" . $self->{wormbase}->{species};
  $self->{'RNASeqSRADir'}    = $self->RNASeqBase . "/SRA";
  $self->{'RNASeqGenomeDir'} = $self->RNASeqBase . "/Genome";
  $self->{'Software'}        = "/nfs/panda/ensemblgenomes/wormbase/software/packages";
#  $self->{'alignmentDir'}    = "tophat_out"; # use this for tophat
  $self->{'alignmentDir'}    = "star_out"; # use this for STAR
  return $self;
}

# getter routines for paths etc.
sub RNASeqBase        { my $self = shift; return $self->{'RNASeqBase'}; }
sub RNASeqSRADir      { my $self = shift; return $self->{'RNASeqSRADir'}; }
sub RNASeqGenomeDir   { my $self = shift; return $self->{'RNASeqGenomeDir'}; }
sub new_genome        { my $self = shift; return $self->{'new_genome'}; }
sub check             { my $self = shift; return $self->{'check'}; }


#####################################################################################################
# Config file methods

=head2 

    Title   :   read_accession
    Usage   :   $data = $self->read_accession($accession);
    Function:   Returns data about the runs in a run, experiment or study by their accession by querying the ENA
    Returns :   ref to hash of hash: primary keys are the run accessions, secondary are types of information about each run
    Args    :   Accession ID

=cut


sub read_accession {
  my ($self, $accession) = @_;
  if (! defined $accession) {die("RNASeq: accession not given in read_run()\n")};

  my %data;

  # For searchable fields, see the last table on http://www.ebi.ac.uk/ena/data/warehouse/usage
  my @fields = (
		"run_accession",
#		"study_accession", # don't read the study_accession treat the secondary_study_accession as the study_accession
		"secondary_study_accession",
		"sample_accession",
		"experiment_accession",
		"submission_accession",
		"center_name",
		"experiment_alias",
		"experiment_title",
		"fastq_ftp",           # may contain two paired end filenames ;-delimited
		"instrument_model",
		"instrument_platform",
		"library_layout",
		"library_name",
		"library_selection",
		"library_source",
		"library_strategy",
		"run_alias",
		"scientific_name",
		"study_alias",
		"study_title",
		"tax_id",
	       );

  my $fields = join ',', @fields;

  my $query = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$accession&result=read_run&fields=$fields";
    
  open (DATA, "wget -q -O - '$query' |") || die("RNASeq: Can't get information on SRA entry $accession in read_accession()\n");

  my $line_count=0;
  while (my $line = <DATA>) {
    if (++$line_count == 1) {next;} # skip the title line
    chomp $line;
    
    my @f = split /\t/, $line;

    # check the data if something goes wrong with the ENA warehouse
    if ($#f != $#fields) {
      my $log = $self->{log};
      $log->write_to("RNASeq: problem found in read_accession() when querying accession: $accession\n");
      if ($#f < $#fields) {
	$log->write_to("RNASeq: Fewer fields written than expected:\n$line\n");
      } else {
	$log->write_to("RNASeq: More fields written than expected:\n$line\n");
      }
      # do the search again using the run accession and see which fields are not returned
      $log->write_to("RNASeq: Checking each field ...\n");
      my $run_accession = $f[0];
      my @working_fields;
      foreach my $test_field (@fields) {
	my $test_query = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$run_accession&result=read_run&fields=$test_field";    
	open (TEST_DATA, "wget -q -O - '$test_query' |") || die("RNASeq: Can't get test_query information on SRA entry $accession in read_accession()\n");	
	my $test_line_count=0;
	my $field_ok = 1;
	while(my $test_line = <TEST_DATA>) {
	  if (++$test_line_count == 1) {next;} # skip the title line    
	  chomp $test_line;
	  my @test_f = split /\t/, $test_line;
	  if ($#test_f != 0) {
	    $log->write_to("RNASeq: Found: '$test_line' returned for field: $test_field\n");
	    $field_ok = 0; # problem
	  }
	}
	if ($field_ok) {
	  push @working_fields, $test_field;
	}
      }
      # now do the query again with the fields that worked
      $fields = join ',', @working_fields;
      $query = "http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$run_accession&result=read_run&fields=$fields";
      open (WORKING_DATA, "wget -q -O - '$query' |") || die("RNASeq: Can't get information on SRA entry $accession in read_accession()\n");
      my $working_line_count=0;
      while(my $working_line = <WORKING_DATA>) {
	if (++$working_line_count == 1) {next;} # skip the title line
	chomp $working_line;    
	my @working_f = split /\t/, $working_line;
	my $working_field_count = 0;
	foreach my $working_field_name (@working_fields) {
	  $data{$f[0]}{$working_field_name} = $f[$working_field_count];
	  if ($working_field_name eq 'secondary_study_accession') {$data{$f[0]}{'study_accession'} = $f[$working_field_count];} # change secondary_study_accession to study_accession
	  $working_field_count++;
	}
      }
      $log->write_to("RNASeq: Problem with missing values now fixed.\n");

    } else { # no problem with the number of fields returned

      # the primary key used is the Run_accession
      my $field_count = 0;
      foreach my $field_name (@fields) {
	$data{$f[0]}{$field_name} = $f[$field_count];
	if ($field_name eq 'secondary_study_accession') {$data{$f[0]}{'study_accession'} = $f[$field_count];} # change secondary_study_accession to study_accession
	# debug Study IDs
	if ($field_name eq 'secondary_study_accession' && $f[$field_count] !~ /^(D|E|S)RP/) {print "Have a secondary_study_accession field in read_accession('$accession') = $f[$field_count]\n"}
	$field_count++;
      }
    }
    if (!exists $data{$f[0]}{'study_accession'} ) {die "No study accession found for $f[0]\n"} # sanity check that we used secondary_study_accession correctly to replace study_accession
  }

  return \%data;
}


=head2 

    Title   :   get_pubmed
    Usage   :   %pubmed = $self->get_pubmed();
    Function:   get the xref of study ID and PubMed ID from the ENA
    Returns :   hash of pubmed IDs keyed by primary study ID (study_alias=PRJNA259320)
    Args    :   

  Get Pubmed for each Study (primary 'PRJN' accession)
 http://www.ebi.ac.uk/ena/data/xref/search?source=PubMed&target=study
  Get EuropePMC for each study gives no results, so ignore this
 http://www.ebi.ac.uk/ena/data/xref/search?source=EuropePMC&target=study

=cut


sub get_pubmed {

  my $query = "http://www.ebi.ac.uk/ena/data/xref/search?source=PubMed&target=study";
    
  open (DATA, "wget -q -O - '$query' |") || die("RNASeq: Can't get information on PubMed xrefs in get_pubmed()\n");

  my %pubmed;
  my $line_count=0;
  while (my $line = <DATA>) {
    if (++$line_count == 1) {next;} # skip the title line

#Source	Source primary accession	Source secondary accession	Target	Target primary accession	Target secondary accession
#PubMed	10023641		study	PRJNA75501	
    my ($source, $pubmed, $secondary, $target, $study) = split /\s/, $line;
    $pubmed{$study}=$pubmed;

  }
  close(DATA);

  return %pubmed;
}

=head2 

    Title   :   get_WBPaper
    Usage   :   %wbpaper = $self->get_WBPaper();
    Function:   get the xref of WBPaper ID and PubMed ID from Caltech
    Returns :   hash of WBPaper IDs keyed by PubMed ID
    Args    :   

=cut

sub get_WBPaper {

  my $query = "http://tazendra.caltech.edu/~azurebrd/cgi-bin/forms/generic.cgi?action=WpaXref";
    
  open (DATA, "wget -q -O - '$query' |") || die("RNASeq: Can't get information on PubMed xrefs in get_WBPaper()\n");

  my %wbpaper;
  while (my $line = <DATA>) {

#WBPaper00000044	cgc44
#WBPaper00000044	doi10.1163/187529267X00058
#WBPaper00000045	cgc45
#WBPaper00000045	doi10.1007/BF01896577
#WBPaper00000045	pmid5800143
#WBPaper00000046	cgc46
#WBPaper00000046	pmid19322369
    my ($wbpaper, $pubmed) = split /\s+/, $line;
    if ($pubmed =~ /pmid(\d+)/) { # ignore anything that is not a PubMed ID
      $wbpaper{$1} = $wbpaper;
    }
  }
  close(DATA);

  return %wbpaper;
}

=head2 

    Title   :   get_experiments_from_run_data
    Usage   :   $data = $self->get_experiments_from_run_data($run_data_hash)
    Function:   converts the run_data_hash from read_accession() into an experiment_hash (each experiment hash is simply pointing to one of its constituent run hashes)
    Returns :   ref to hash of experiment data keyed by the experiment accession
    Args    :   ref to hash of the run data from read_accession()

=cut


sub get_experiments_from_run_data {
  my ($self, $run_hashref) = @_;

  my %experiment_data;
  
  foreach my $run (keys %{$run_hashref}) {
    my $run_count = 0;
    my $expt = $run_hashref->{$run}{'experiment_accession'};
    if (exists $experiment_data{$expt}) { # get the count of the number of times this experiment has been seen before == the number of runs for this experiment seen already
      $run_count = $experiment_data{$expt}{'run_count'}
    }
    $run_count++;
    $experiment_data{$expt} = $run_hashref->{$run};
    $experiment_data{$expt}{'run_count'} = $run_count;
  }

  return \%experiment_data;
}

=head2 

    Title   :   find_studies
    Usage   :   $data = $self->find_studies
    Function:   Finds all studies for the current species by querying the ENA
    Returns :   ref to hash of hash: primary keys are the study accessions, secondary are types of information about each study
    Args    :   

=cut


sub find_studies {
  my ($self) = @_;
  my $taxon = $self->{wormbase}->ncbi_tax_id;
  if (! defined $taxon) {die("RNASeq: taxon ID not given in find_studies()\n")};

  my %data;
  my $data = \%data; # hashref
  $self->find_studies_search($taxon, $data);

  # C. briggsae has some libraries with the AF16 strain taxon_id so search for this as well
  if ($taxon == 6238) {$self->find_studies_search(473542, $data);} 

  return $data;
}

=head2 

    Title   :   find_studies_search
    Usage   :   $data = $self->find_studies_search($taxon_id, $hashref)
    Function:   Finds all studies for the given species taxon_id by querying the ENA
    Returns :   ref to hash of hash: primary keys are the study accessions, secondary are types of information about each study
    Args    :   taxon_id to search for, hashref to populate

=cut


sub find_studies_search {
  my ($self, $taxon_id, $data) = @_;

  my $query = 'http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22tax_eq%28TAXON%29%22&result=read_study&display=XML';
  $query =~ s/TAXON/$taxon_id/;

  my $ua = LWP::UserAgent->new();
  my $response = $ua->get($query);
  my $xmlstring = $response->content;

# the following XML parser stuff doesn't work as the xmlstring is not parsed correctly - it doesn't like having two lines like:
# <?xml version="1.0" encoding="UTF-8"?>

#  my $parser = XML::LibXML->new();
#  print "parsing string:\n\n";
#  my $dom = $parser->parse_string($xmlstring);
#  print "dom->toSring\n\n";
#  print $dom->toString();
#  print "findnodes:\n\n";
#  for my $project ($dom->findnodes('//*')) { # was //PROJECT
#    print $project->nodeName(), ": ", $project->textContent(), "\n";
#    print "\n";
#  }

  my $primary_id;
  my $pubmed=0;
  my @xmlstring = split /\n/, $xmlstring;
  foreach my $line (@xmlstring) {
#    if ($line =~ /<PRIMARY_ID>(\S+)<\/PRIMARY_ID>/) {
#      $primary_id = $1;

      # the following are for PROJECT sections
#    } els
    if ($line =~ /<SECONDARY_ID>(\S+)<\/SECONDARY_ID>/) {
      $primary_id = $1;
      $data->{$primary_id}{secondary_id} = $1;
      # Study ID debug
      if ($data->{$primary_id}{secondary_id} !~ /^(D|S|E)RP/) {print "Have a secondary_id field in find_studies_search() = ",$data->{$primary_id}{secondary_id},"\n"}

    } elsif ($line =~ /<TITLE>(.+?)<\/TITLE>/) {
      $data->{$primary_id}{title} = $1;
    } elsif ($line =~ /<DESCRIPTION>(.+?)<\/DESCRIPTION>/) {
      $data->{$primary_id}{description} = $1;
    } elsif ($pubmed) {
      $line =~ /<ID>(\d+)<\/ID>/;
      $data->{$primary_id}{pubmed} = $1;
      $pubmed = 0;
    } elsif ($line =~ /<DB>PUBMED<\/DB>/) {
      # get next line
      $pubmed=1;

      # the following are for STUDY sections
    } elsif ($line =~ /<STUDY_TITLE>(.+?)<\/STUDY_TITLE>/) {
      $data->{$primary_id}{title} = $1;
    } elsif ($line =~ /<STUDY_ABSTRACT>(.+?)<\/STUDY_ABSTRACT>/) {
      $data->{$primary_id}{description} = $1;
    } 
  }
  return $data;
}


=head2 

    Title   :   find_experiments
    Usage   :   $data = $self->find_experiments
    Function:   Finds all experiments for the current species by querying the ENA - gets basic details of the Study accession and title
    Returns :   ref to hash of hash: primary keys are the experiment accessions, secondary are types of information about each experiment
    Args    :   

=cut


sub find_experiments {
  my ($self) = @_;
  my $taxon = $self->{wormbase}->ncbi_tax_id;
  if (! defined $taxon) {die("RNASeq: taxon ID not given in find_studies()\n")};
  
  my %data;
  my $data = \%data; # hashref
  $self->find_experiments_search($taxon, $data);
  
  # C. briggsae has some libraries with the AF16 strain taxon_id so search for this as well
  if ($taxon == 6238) {$self->find_experiments_search(473542, $data);} 
  
  return $data;
}

=head2 

    Title   :   find_experiments_search
    Usage   :   $data = $self->find_experiments_search($taxon_id, $hashref)
    Function:   Finds all experiments for the given species taxon_id by querying the ENA
    Returns :   ref to hash of hash: primary keys are the experiment accessions, secondary are types of information about each experiment
    Args    :   taxon_id to search for, hashref to populate

=cut


sub find_experiments_search {
  my ($self, $taxon_id, $data) = @_;
  
  my $query = 'http://www.ebi.ac.uk/ena/data/warehouse/search?query=%22tax_eq%28TAXON%29%22&result=read_experiment&display=XML';
  $query =~ s/TAXON/$taxon_id/;
  
  open (DATA, "wget -q -O - '$query' |") || die("RNASeq: Can't get information on SRA entries in find_experiments_search()\n");

  my $primary_id;
  my $study_id=0;
  my $in_study=0;
  my $in_geo=0;
  
  #     <IDENTIFIERS>
  #          <PRIMARY_ID>SRX007173</PRIMARY_ID>
  #          <SUBMITTER_ID namespace="University of Toronto">N2 Mixed Stage</SUBMITTER_ID>
  #     </IDENTIFIERS>
  #     <TITLE>Illumina Genome Analyzer paired end sequencing; A high resolution transcriptome map for both wild-type and NMD defective C. elegans</TITLE>
  #     <STUDY_REF accession="SRP001010" refname="Elegans Transcriptome" refcenter="University of Toronto">
  #          <IDENTIFIERS>
  #               <PRIMARY_ID>SRP001010</PRIMARY_ID>
  #               <SUBMITTER_ID namespace="University of Toronto">Elegans Transcriptome</SUBMITTER_ID>
  #          </IDENTIFIERS>
  #     </STUDY_REF>
  

  #   <EXPERIMENT_ATTRIBUTES>
  #        <EXPERIMENT_ATTRIBUTE>
  #             <TAG>GEO Accession</TAG>
  #             <VALUE>GSM1576710</VALUE>
  #        </EXPERIMENT_ATTRIBUTE>
  
#  my $trace='SRX494916'; # for debugging
#  my $trace_on=0;

  my $line_count=0;
  while (my $line = <DATA>) {

#    if ($trace_on) {print $line}

    if ($line =~ /<PRIMARY_ID>(..X\S+)<\/PRIMARY_ID>/ && !$in_study) {
      $primary_id = $1;
#      if ($primary_id eq $trace) {$trace_on=1;print $line} else {$trace_on=0;}
    } elsif ($line =~ /<STUDY_REF /) { # the study accession should be somewhere in this line
      $in_study=1;
      if ($line =~ /accession=\"(\S+)\"/) {
	$study_id=$1;
	$data->{$primary_id}{study_accession} = $1;
	if ($study_id !~ /^(D|E|S)RP/) {print "Have an unusual STUDY_REF accession field in find_experiments_search() = $study_id\n"}
      }
    } elsif ($line =~ /<\/STUDY_REF>/) {
      $in_study=0;
    } elsif ($line =~ /<TITLE>(.+?)<\/TITLE>/) {
      $data->{$primary_id}{study_title} = $1;
    } elsif ($line =~ /<TAG>GEO Accession<\/TAG>/) {
      $in_geo=1;
    } elsif ($in_geo) {
      $in_geo=0;
      if ($line =~ /<VALUE>(GSM\d+)<\/VALUE>/) {
	$data->{$primary_id}{geo_accession} = $1;
      }
    }
  }
  close(DATA);

  if (!defined $data) {die("RNASeq: Can't get information on SRA entries in find_experiments_search()\n")}
  return $data;
}

=head2 

    Title   :   read_all_studies_config
    Usage   :   $ini = $self->read_all_studies_config()
    Function:   Gets the Studies in the current species from the INI config files
    Returns :   ref to INI object
    Args    :   

=cut


sub read_all_studies_config {
  my ($self) = @_;

  my $path = $self->{wormbase}->misc_dynamic . "/SHORT_READS/" . $self->{wormbase}->{species} . "/Studies.ini";
  if (!-e $path) {system("touch $path")}
  my $ini = Config::IniFiles->new( -file => $path, -allowempty => 1);
  return $ini;
}


=head2 

    Title   :   find_new_studies
    Usage   :   @new_study_accessions = $self->find_new_studies($study_ini, $studies)
    Function:   Gets a list of the accessions of the new studies that are not in the 'studies' INI config file
    Returns :   array of new accession names
    Args    :   ini object of known Studies, ref to hash of all Studies

=cut


sub find_new_studies {
  my ($self, $ini, $studies) = @_;

  my @new_study_accessions;

  foreach my $study (keys %{$studies}) {
    if (! $ini->SectionExists($study)) {
      push @new_study_accessions, $study;
    }
  }

  return @new_study_accessions;
}

=head2 

    Title   :   add_new_studies
    Usage   :   $ini = $self->add_new_studies($study_ini, $studies, @new_study_accessions)
    Function:   Updates the 'study' INI object and the 'study' INI file with details of new studies
    Returns :   updated ini object and ini file
    Args    :   ini object of known Studies, ref to hash of all Studies, array of names of new studies

=cut


sub add_new_studies {
  my ($self, $ini, $studies, @new_study_accessions) = @_;

  foreach my $study_accession (@new_study_accessions) {
    print "\nAdding new Study $study_accession\n";
    $ini->AddSection ($study_accession);
    $ini->RewriteConfig;
    foreach my $param (keys %{$studies->{$study_accession}}) {
      my $value = $studies->{$study_accession}{$param};
      $ini->newval($study_accession, $param, $value);
    }
    $ini->RewriteConfig;
    # now write out the experiments details for each new study
    $self->add_new_experiments_in_study($study_accession, $study_accession);
  }

  return $ini;
}


=head2 

    Title   :   get_study_config_file
    Usage   :   my $ini = $self->get_study_config_file($study_accession);
    Function:   opens (creates if necessary) the config file holding the experiments details of a study
    Returns :   config INI handle
    Args    :   study accession

=cut

sub get_study_config_file {

  my ($self, $study_accession) = @_;

  my $path = $self->{wormbase}->misc_dynamic . "/SHORT_READS/" . $self->{wormbase}->{species} . "/${study_accession}.ini";
  if (!-e $path) {system("touch $path")}
  my $ini = Config::IniFiles->new( -file => $path, -allowempty => 1 );
  return $ini;

}

=head2 

    Title   :   add_new_experiments_in_study
    Usage   :   $self->add_new_experiments_in_study($study_accession)
    Function:   Updates the 'experiments' INI object for a study and writes the study's 'experiments' INI file
    Returns :   updates Experiments ini file
    Args    :   Accession of the new Study (used to name INI file) to search for in ENA to get Experiment data

=cut


sub add_new_experiments_in_study {
  my ($self, $study_accession) = @_;
  print "\n";

  my $taxon = $self->{wormbase}->ncbi_tax_id;
  if (! defined $taxon) {die("RNASeq: taxon ID not given in find_studies()\n")};

  my $runs_in_this_study = $self->read_accession($study_accession);
  my $experiments = $self->get_experiments_from_run_data($runs_in_this_study);

  foreach my $experiment_accession (keys %{$experiments}) {
    # no , don't do the following. We now want to store everything under the Study ID, not the Project ID
    # if we are using a secondary accession as the Study name for the file, then use the primary study accession instead
    #if ($experiments->{$experiment_accession}{study_accession} ne $study_accession) {$study_accession = $experiments->{$experiment_accession}{study_accession}}

    my $ini = $self->get_study_config_file($study_accession);

    # make sure that this experiment was done in this species - some Studies have multiple species in
    if ($experiments->{$experiment_accession}{tax_id} ne $taxon) {
      if ($taxon == 6238) { # C. briggsae has some libraries with the AF16 strain taxon_id so search for this as well
	 if ($experiments->{$experiment_accession}{tax_id} ne '473542') {next}
       } else {
	 next;
       }
    } 

    print "Adding new Experiment $experiment_accession to Study $study_accession\n";
    $ini->AddSection ($experiment_accession);
    $ini->newval($experiment_accession, 'quality', 'phred');   # assume that all new experiments used phred quality format - this can be manually edited in the INI file

    foreach my $param (keys %{$experiments->{$experiment_accession}}) {
      my $value = $experiments->{$experiment_accession}{$param};
      $ini->newval($experiment_accession, $param, $value);
    }
    $ini->RewriteConfig;
  }
}

=head2 

    Title   :   add_one_new_experiment_to_config
    Usage   :   $self->add_one_new_experiment_to_config($study_accession, $experiment_accession, $geo_accession, $study_ini, $pubmed, $wbpaper);
    Function:   Updates the 'experiments' INI object for a study and writes the study's 'experiments' INI file
    Returns :   updates Experiments ini file
    Args    :   Accession of the new Study (used to name INI file), experiment accession

=cut


sub add_one_new_experiment_to_config {
  my ($self, $study_accession, $experiment_accession, $geo_accession, $study_ini, $pubmed, $wbpaper) = @_;

  my $taxon = $self->{wormbase}->ncbi_tax_id;
  if (! defined $taxon) {die("RNASeq: taxon ID not given in find_studies()\n")};

  my $runs_in_this_study = $self->read_accession($study_accession);
  my $experiments = $self->get_experiments_from_run_data($runs_in_this_study);

  my $ini = $self->get_study_config_file($study_accession);

  # make sure that this experiment was done in this species - some Studies have multiple species in
  if ($experiments->{$experiment_accession}{tax_id} ne $taxon) {
    if ($taxon == 6238) { # C. briggsae has some libraries with the AF16 strain taxon_id so search for this as well
      if ($experiments->{$experiment_accession}{tax_id} ne '473542') {return}
    } else {
      return;
    }
  } 

  print "Adding new Experiment $experiment_accession to Study $study_accession\n";
  $ini->AddSection ($experiment_accession);
  $ini->newval($experiment_accession, 'quality', 'phred');   # assume that all new experiments used phred quality format - this can be manually edited in the INI file

  foreach my $param (keys %{$experiments->{$experiment_accession}}) {
    my $value = $experiments->{$experiment_accession}{$param};
    $ini->newval($experiment_accession, $param, $value);
  }

  # we got the GEO_accession from find_experiments() getting the complete list of experiments for this species
  if (defined $geo_accession && $geo_accession ne '') {$ini->newval($experiment_accession, 'geo_accession', $geo_accession)}

  $ini->RewriteConfig;
  

  # now update the Study.ini record based on this experiment information
    if (exists $experiments->{$experiment_accession}{study_alias}) {
      my $study_alias = $experiments->{$experiment_accession}{study_alias};
      my $pubmed = $pubmed->{$study_alias};
      if (defined $pubmed) {
	$study_ini->newval($study_accession, 'pubmed', $pubmed);
	my $wbpaper = $wbpaper->{$pubmed};
	if (defined $wbpaper) {
	  $study_ini->newval($study_accession, 'wbpaper', $wbpaper);
	}
      }
    }
    # library_source=GENOMIC
    if ($experiments->{$experiment_accession}{library_source} eq 'GENOMIC') {$study_ini->newval($study_accession, 'ignore', 'Genomic.')}

    # library_strategy=ChIP-Seq
    if ($experiments->{$experiment_accession}{library_strategy} eq 'ChIP-Seq') {$study_ini->newval($study_accession, 'ignore', 'ChIP-Seq.')}

    # library_selection=size fractionation
    if ($experiments->{$experiment_accession}{library_selection} eq 'size fractionation') {$study_ini->newval($study_accession, 'ignore', 'Small RNA stuff.')}

    $study_ini->RewriteConfig;

}



=head2 

    Title   :   update_experiment_config_record
    Usage   :   $ini = $self->read_experiments_from_study_config($study_accession, $experiment_accession, $geo_accession, $study_ini, \%pubmed, \%wbpaper, %expt_config)
    Function:   Updates an existing experiment config record adding new parameters, unless the config parameter 'locked' is set. Changed config parameters are reported to STDOUT but not altered.
    Returns :   updates Experiments ini file
    Args    :   Accession of the new Study (used to name INI file), Experiment_accession to update, hash of existing experiment parameter/value pairs

=cut

sub update_experiment_config_record {
  my ($self, $study_accession, $experiment_accession, $geo_accession, $study_ini, $pubmed, $wbpaper, %expt_config) = @_;

  my $changed_study_id = 0;
  my $changed_ini=0;
  my $changed_study=0;

  # we just want to add anything that may have been added to the
  # experiment record, unless it has known problems that have been
  # manually corrected, in which case we should have put a keyword
  # 'locked' in the record.

  if (exists $expt_config{locked}) {return}

  # anything that has changed probably needs to be manually looked at,
  # so it is reported to STDOUT and not changed
  
  my $taxon = $self->{wormbase}->ncbi_tax_id;
  if (! defined $taxon) {die("RNASeq: taxon ID not given in find_studies()\n")};

  my $runs_in_this_study = $self->read_accession($study_accession);
  my $experiments = $self->get_experiments_from_run_data($runs_in_this_study);

  # if we are using a primary accession as the Study name for the file, then use the secondary study accession instead
  if (exists $experiments->{$experiment_accession}{secondary_study_accession} ) {
    if ($experiments->{$experiment_accession}{secondary_study_accession} ne $study_accession) {
      print "In update_experiment_config_record() :\nExpt $experiment_accession has different config file secondary study accession (",$experiments->{$experiment_accession}{secondary_study_accession},") to the ENA secondary_study_accession ($study_accession)\n";
      $changed_study_id = 1;
      $study_accession = $experiments->{$experiment_accession}{study_accession}
    }
  }
  my $ini;

  # make sure that this experiment was done in this species - some Studies have multiple species in
  if ($experiments->{$experiment_accession}{tax_id} ne $taxon) {
    if ($taxon == 6238) { # C. briggsae has some libraries with the AF16 strain taxon_id so search for this as well
      if ($experiments->{$experiment_accession}{tax_id} ne '473542') {return}
    } else {
      return;
    }
  } 
  
  if ($changed_study_id) { # the existing experiment record may be using the old primary ID, so write it all out
    foreach my $param (keys %{$experiments->{$experiment_accession}}) {
      my $value = $experiments->{$experiment_accession}{$param};
      if (!$changed_ini) {$ini = $self->get_study_config_file($study_accession);} # only get the details from the INI file if we wish to change it
      $ini->newval($experiment_accession, $param, $value);
      $changed_ini = 1;
    }
    
  } else { # be selective in what we write
    foreach my $param (keys %{$experiments->{$experiment_accession}}) {
      my $value = $experiments->{$experiment_accession}{$param};
      if (!exists $expt_config{$param}) {
	print "Updating Experiment $experiment_accession to Study $study_accession: adding '$param = $value\n";
	if (!$changed_ini) {$ini = $self->get_study_config_file($study_accession);} # only get the details from the INI file if we wish to change it
	$ini->newval($experiment_accession, $param, $value);
	$changed_ini = 1;
      } else {
	my $old_value = $expt_config{$param};
	if ($old_value ne $value) {
	  if ($param eq 'fastq_ftp' && defined $old_value && $old_value ne '') {
            # do nothing
	  } elsif ($param eq 'library_selection' || $param eq 'library_strategy' || $param eq 'fastq_ftp') {
	    print "Changed value '$param = $value' found for Experiment $experiment_accession in Study $study_accession - old value: '$param = $old_value'\n";
	  } else {
	    if (!$changed_ini) {$ini = $self->get_study_config_file($study_accession);} # only get the details from the INI file if we wish to change it
	    $ini->newval($experiment_accession, $param, $value);
	    $changed_ini = 1;
	  }
	}
      }
    }
  }

  # we got the GEO_accession from find_experiments() getting the complete list of experiments for this species
  if (defined $geo_accession && $geo_accession ne '') {
    if (!$changed_ini) {$ini = $self->get_study_config_file($study_accession);} # only get the details from the INI file if we wish to change it
    $ini->newval($experiment_accession, 'geo_accession', $geo_accession);
    $changed_ini = 1;
  }

  if ($changed_ini) {
    $ini->RewriteConfig;
  }

  # now update the Study.ini record based on this experiment information
  my $pubmed_id = undef;
  if (exists  $expt_config{pubmed}) {
    $pubmed_id = $expt_config{pubmed};
  } elsif (exists $experiments->{$experiment_accession}{study_alias}) {
    my $study_alias = $experiments->{$experiment_accession}{study_alias};
    $pubmed_id = $pubmed->{$study_alias};
    if (defined $pubmed_id) {
      $study_ini->newval($study_accession, 'pubmed', $pubmed_id);
      $changed_study = 1;
    }
  }
  if (defined $pubmed_id) {
    if (!exists $expt_config{wbpaper}) {
      my $wbpaper = $wbpaper->{$pubmed_id};
      if (defined $wbpaper) {
	$study_ini->newval($study_accession, 'wbpaper', $wbpaper);
	$changed_study = 1;
      }
    }
  }
  

  if (!exists $expt_config{ignore}) {
    # library_source=GENOMIC
    if ($experiments->{$experiment_accession}{library_source} eq 'GENOMIC') {$changed_study = 1; $study_ini->newval($study_accession, 'ignore', 'Genomic.')}
  
    # library_strategy=ChIP-Seq
    if ($experiments->{$experiment_accession}{library_strategy} eq 'ChIP-Seq') {$changed_study = 1; $study_ini->newval($study_accession, 'ignore', 'ChIP-Seq.')}
  
    # library_selection=size fractionation
    if ($experiments->{$experiment_accession}{library_selection} eq 'size fractionation') {$changed_study = 1; $study_ini->newval($study_accession, 'ignore', 'Small RNA stuff.')}
  }

  if ($changed_study) {
    $study_ini->RewriteConfig;
  }
}

=head2 

    Title   :   read_experiments_from_study_config
    Usage   :   $ini = $self->read_experiments_from_study_config($study_accession)
    Function:   Gets Experiments in a specified Study in the current species
    Returns :   Experiments ini object of experiments in the specified Study
    Args    :   study accession

=cut


sub read_experiments_from_study_config {
  my ($self, $study_accession) = @_;

  my $path = $self->{wormbase}->misc_dynamic . "/SHORT_READS/" . $self->{wormbase}->{species} . "/${study_accession}.ini";
  if (!-e $path) {die "RNASeq: Can't find file $path\n"}
  my $ini = Config::IniFiles->new( -file => $path, -allowempty => 1 );

  return $ini;
}



=head2 

    Title   :   convert_ini_to_hashref
    Usage   :   $data = $self->convert_ini_to_hashref($ini_object)
    Function:   converts and APPENDS an INI object onto a hashref object
    Returns :   hashref of the contents of the INI object
    Args    :   INI object, hashref to append to

=cut


sub convert_ini_to_hashref {
  my ($self, $ini, $hashref) = @_;

  if (defined $ini) {
    my @sections = $ini->Sections;
    foreach my $section (@sections) {
      my @params = $ini->Parameters($section);
      foreach my $param (@params) {
	my $value = $ini->val($section, $param);
	$hashref->{$section}{$param} = $value;
      }
    }
  }
  return $hashref;
}



=head2 

    Title   :   get_all_experiments
    Usage   :   $data = $self->get_all_experiments()
    Function:   gets all of the experiments data for the current species as a hashref of hash
    Returns :   ref to hash of hash
    Args    :   

=cut


sub get_all_experiments {
  my ($self) = @_;

  my %hash;
  my $hashref = \%hash;

  my $studies_ini = $self->read_all_studies_config();
  my @studies = $studies_ini->Sections;
  foreach my $study_accession (@studies) {
    my $experiments_ini = $self->read_experiments_from_study_config($study_accession);
    $self->convert_ini_to_hashref($experiments_ini, $hashref);
  }

  return $hashref;
}


=head2 

    Title   :   update_study_config_data
    Usage   :   @studies_changed = $self->update_study_config_data()
    Function:   queries the ENA and updates the study and experiment config files if it finds new studies
                this is now deprecated as it pulls out too much Project data insstead of Study data
    Returns :   array of accessions of studies changed
    Args    :   

=cut


sub update_study_config_data {
  my ($self) = @_;

  my $studies = $self->find_studies();
  my $ini = $self->read_all_studies_config();
  my @new_studies = $self->find_new_studies($ini, $studies);
  $ini = $self->add_new_studies($ini, $studies, @new_studies);

  return @new_studies;
}


=head2 

    Title   :   update_experiment_config_data
    Usage   :   @studies_changed = $self->update_experiment_config_data()
    Function:   queries the ENA and updates the study and experiment config files if it finds new experiments
    Returns :   array of accessions of experiments changed
    Args    :   

=cut


sub update_experiment_config_data {
  my ($self) = @_;

  my $ena_experiments = $self->find_experiments(); # hash key = all experiment_accessions, hash of hash {study_accession} or {study_title} or {geo_accession}
  my $config_experiments = $self->get_all_experiments();

  my @new_experiment_accessions;
  my @existing_experiment_accessions;

  my %pubmed = $self->get_pubmed(); # hash of pubmed IDs keyed by primary study ID (study_alias=PRJNA259320)
  my %wbpaper = $self->get_WBPaper(); # hash of WBPaper IDs keyed by PubMed ID

  foreach my $experiment (keys %{$ena_experiments}) {
    if (!exists $config_experiments->{$experiment}) {
      push @new_experiment_accessions, $experiment;
    } else {
      push @existing_experiment_accessions, $experiment;
    }
  }

  # create new studies
  my $study_ini = $self->read_all_studies_config();

  foreach my $experiment_accession (@new_experiment_accessions) {
    # store study details in INI file
    my $study_accession = $ena_experiments->{$experiment_accession}{study_accession};
    my $study_title = $ena_experiments->{$experiment_accession}{study_title};
    my $geo_accession = (exists $ena_experiments->{$experiment_accession}{geo_accession}) ? $ena_experiments->{$experiment_accession}{geo_accession} : '';
    $study_ini->AddSection ($study_accession);
    $study_ini->newval($study_accession, 'study_title', $study_title);
    print "Added Study $study_accession to Study.ini\n";

    # get full experiment details from ENA and store experiment details in INI file
    $self->add_one_new_experiment_to_config($study_accession, $experiment_accession, $geo_accession, $study_ini, \%pubmed, \%wbpaper);
  }

  # add new experiments to existing studies
  foreach my $experiment_accession (@existing_experiment_accessions) {
    # we don't want to change the Study INI record much, 
    # we just want to add anything that may have been added or changed
    # to the experiment record, unless it has known problems that have
    # been manually corrected in which case we should have put a
    # keyword 'locked' in the record.

    # But it is useful to add pubmed and wbpaper data to the Study.ini
    my $study_accession = $ena_experiments->{$experiment_accession}{study_accession};
    my $geo_accession = (exists $ena_experiments->{$experiment_accession}{geo_accession}) ? $ena_experiments->{$experiment_accession}{geo_accession} : '';

    # get full experiment details from ENA and store experiment details in INI file
    print "Checking existing expt $experiment_accession in study $study_accession\n";
    $self->update_experiment_config_record($study_accession, $experiment_accession, $geo_accession, $study_ini, \%pubmed, \%wbpaper, %{$config_experiments->{$experiment_accession}});
  }

  $study_ini->RewriteConfig;

  return @new_experiment_accessions;
}


=head2 

    Title   :   add_param_to_experiment_config
    Usage   :   $self->add_param_to_experiment_config( $experiment_accession, $param, $value)
    Function:   adds a new param and value to an experiment config INI file
    Returns :   
    Args    :   experiment_accession, param, value

=cut


sub add_param_to_experiment_config {
  my ($self, $experiment_accession, $param, $value) = @_;

  my $experiments = $self->get_all_experiments();
  my $study_accession = $experiments->{$experiment_accession}{study_accession};
  my $experiment_ini = $self->read_experiments_from_study_config($study_accession);
  $experiment_ini->newval($experiment_accession, $param, $value);
  $experiment_ini->RewriteConfig;
}

=head2 

    Title   :   get_transcribed_long_experiments
    Usage   :   $data = $self->get_transcribed_long_experiments()
    Function:   gets the transcribed and not small RNA experiments data for the current species as a hashref of hash
    Returns :   ref to hash of hash
    Args    :   

=cut


sub get_transcribed_long_experiments {
  my ($self) = @_;

  my %hash;
  my $hashref = \%hash;
  my %transcribed;
  my $transcribedref = \%transcribed;

  my $studies_ini = $self->read_all_studies_config();
  my @studies = $studies_ini->Sections;
  foreach my $study_accession (@studies) {
    my $experiments_ini = $self->read_experiments_from_study_config($study_accession);
    $self->convert_ini_to_hashref($experiments_ini, $hashref);
  }

  foreach my $expt (keys %{$hashref}) {
    if (
	($hashref->{$expt}{library_source} eq 'TRANSCRIPTOMIC') &&
	($hashref->{$expt}{library_strategy} eq 'RNA-Seq' || 
	 $hashref->{$expt}{library_strategy} eq 'FL-cDNA' || # full-length cDNA commonly specified when doing 454 sequencing
	 $hashref->{$expt}{library_strategy} eq 'EST') &&
	($hashref->{$expt}{library_selection} ne 'size fractionation') &&
	(! exists $hashref->{$expt}{ignore})
       ) {
      $transcribedref->{$expt} = $hashref->{$expt};
    }
  }

  return $transcribedref;
}

=head2 

    Title   :   get_transcribed_small_experiments
    Usage   :   $data = $self->get_transcribed_small_experiments()
    Function:   gets the transcribed and small RNA experiments data for the current species as a hashref of hash
    Returns :   ref to hash of hash
    Args    :   

=cut


sub get_transcribed_small_experiments {
  my ($self) = @_;

  my %hash;
  my $hashref = \%hash;
  my %small;
  my $smallref = \%small;

  my $studies_ini = $self->read_all_studies_config();
  my @studies = $studies_ini->Sections;
  foreach my $study_accession (@studies) {
    my $experiments_ini = $self->read_experiments_from_study_config($study_accession);
    $self->convert_ini_to_hashref($experiments_ini, $hashref);
  }

  foreach my $expt (keys %{$hashref}) {
    if (
	($hashref->{$expt}{library_source} eq 'TRANSCRIPTOMIC') &&
	($hashref->{$expt}{library_selection} eq 'size fractionation') &&
	(! exists $hashref->{$expt}{ignore})
       ) {
      $smallref->{$expt} = $hashref->{$expt};
    }
  }

  return $smallref;
}

=head2 

    Title   :   filter_experiments
    Usage   :   $filtered_experiments = $self->filter_experiments($experiments, 'library_layout', 'PAIRED')
    Function:   takes a hasref of hash experiments data and returns only those that match the specified condition
    Returns :   ref to hash of hash experiments data
    Args    :   ref to hash of hash experiments, parameter to filter by, value to filter by

=cut


sub filter_experiments {
  my ($self, $experiments, $param, $value) = @_;

  my $filtered;

  foreach my $expt (keys %{$experiments}) {
    if (exists $experiments->{$expt}{$param} && $experiments->{$expt}{$param} eq $value) {
      $filtered->{$expt} = $experiments->{$expt};
    }
  }

  return $filtered;
}

=head2 

    Title   :   exclude_experiments
    Usage   :   $filtered_experiments = $self->exclude_experiments($experiments, 'library_layout', 'PAIRED')
    Function:   takes a hasref of hash experiments data and returns only those that do not match the specified condition
    Returns :   ref to hash of hash experiments data
    Args    :   ref to hash of hash experiments, parameter to exclude with, value to exclude with

=cut


sub exclude_experiments {
  my ($self, $experiments, $param, $value) = @_;

  my $filtered;

  foreach my $expt (keys %{$experiments}) {
    if (!exists $experiments->{$expt}{$param} || $experiments->{$expt}{$param} ne $value) {
      $filtered->{$expt} = $experiments->{$expt};
    }
  }

  return $filtered;
}

#####################################################################################################
# general mapping methods

=head2 

    Title   :   remove_old_experiment_files
    Usage   :   self->remove_old_experiment_files()
    Function:   deletes old experiment files from the directories on /nfs/nobackup, keeping stuff as specified by $new_genome and $check
    Returns :   
    Args    :   experiments - hashref of all experiments to be aligned 

=cut

sub remove_old_experiment_files  {
  my ($self, $experiments) = @_;

  my $log = $self->{log};
  if (! $self->{check}) {
    foreach my $SRX (keys %{$experiments}) {
      if (-d $self->{RNASeqSRADir}."/$SRX") {
	chdir $self->{RNASeqSRADir}."/$SRX";
	if ($self->{new_genome}) {
	  $self->{'wormbase'}->run_command("rm -rf ".$self->{'alignmentDir'}, $log);
	  $self->{'wormbase'}->run_command("rm -rf Introns", $log);
	  $self->{'wormbase'}->run_command("rm -rf SRR", $log);
	}
	$self->{'wormbase'}->run_command("rm -rf cufflinks", $log); # always remove cufflinks unless -check
      }
    }
  }
}

=head2 

    Title   :   strip_prefix_from_fasta_file
    Usage   :   self->strip_prefix_from_fasta_file($input_file, $output_file)
    Function:   remove the chromosome_prefix from the chromosome names in a fasta file
    Reason  :   GBrowse doesn't have the chromosome_prefix and this makes it easier to import BAM files
    Returns :   
    Args    :   name of input file, name of output file

=cut

sub strip_prefix_from_fasta_file {
  my ($self, $input_file, $output_file) = @_;

  my $prefix = $self->{wormbase}->chromosome_prefix;

  open (IN, "<$input_file");
  open (OUT, ">$output_file");
  while (my $line = <IN>) {
    $line =~ s/>${prefix}(.+)/>$1/;
    print OUT $line;
  }
  close (OUT);
  close (IN);

}

=head2 

    Title   :   get_SRA_files
    Usage   :   $self->get_SRA_files($experiment_accession);
    Function:   gets the SRA fastq files for an experiment and place them in the SRA directory path/SRR dir
    Returns :   
    Args    :   experiment accession

=cut

sub get_SRA_files {
  my ($self, $experiment_accession) = @_;
  
  my $RNASeqSRADir = $self->{RNASeqSRADir};
  my $log = $self->{log};
  my $SRRDir = "$RNASeqSRADir/$experiment_accession/SRR";
  my $srr_done_file = "$SRRDir/fastq.done";
    
  if (!$self->{check} || ! -e $srr_done_file) {
    my $status = $self->{wormbase}->run_command("rm -rf $SRRDir", $log);
    chdir "$RNASeqSRADir";
    $self->get_SRX_file($experiment_accession);
  }
}

=head2 

    Title   :   get_SRX_file
    Usage   :   $self->get_SRX_files($experiment_accession);
    Function:   gets the SRA fastq files for an experiment
    Returns :   success = 0
    Args    :   experiment accession

=cut

sub get_SRX_file {
  my ($self, $experiment_accession) = @_;

  my $RNASeqSRADir = $self->{RNASeqSRADir};
  my $log = $self->{log};
  my $status;

  # in the EBI, the set of files available for download from an
  # experiment can be seen by parsing the page resulting from a FTP
  # query like:
  # http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=ERX200392&result=read_run&fields=experiment_accession,run_accession,fastq_ftp
  
  
  # The EBI specify the preferred method of download:
  # "Aspera is a commercial file transfer protocol that provides better
  # transfer speeds than ftp over long distances. For short distance file
  # transfers we continue to recommend the use of ftp."
  # See: http://www.ebi.ac.uk/ena/about/sra_data_download
  # See: https://www.ebi.ac.uk/ena/about/browser
  #
  # Aspera doesn't work from the LSF farm.

# get page like:
#experiment_accession	run_accession	fastq_ftp
#ERX200392	ERR225732	ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225732/ERR225732_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR225/ERR225732/ERR225732_2.fastq.gz

  my $failed = 0;
  
  if (!-e $experiment_accession) {
    unless(mkdir $experiment_accession, 0777) {
      $log->log_and_die( "ERROR: Could not mkdir subdir $experiment_accession: $!\n");
    }
  }
  chdir $experiment_accession;
  if (!-e "SRR") {
    unless(mkdir "SRR", 0777) {
      $log->log_and_die( "ERROR: Could not mkdir subdir SRR: $!\n");
    }
  }
  chdir "SRR";

  my $count;
  open(ENTRY, "/sw/arch/bin/wget -q -O - 'http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=$experiment_accession&result=read_run&fields=experiment_accession,run_accession,fastq_ftp' |") || $log->log_and_die("Can't get information on SRA entry $experiment_accession\n");
  while (my $line = <ENTRY>) {
    chomp $line;
    if ($line =~ /^experiment_accession/) {next}
    my @line = split /\s+/, $line;
    my $run_accession = $line[1];
    my $ftp = $line[$#line];

    # use Aspera to fetch the file - this doesn't work on the LSF cluster
    # command like: ~gw3/.aspera/connect/bin/ascp -QT -l 300m -i /net/nas17b/vol_homes/homes/gw3/.aspera/connect/etc/asperaweb_id_dsa.putty era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR941/SRR941702/SRR941702.fastq.gz .
    #$address =~ s#ftp://ftp.sra.ebi.ac.uk#era-fasp\@fasp.sra.ebi.ac.uk:#; # change FTP address to Aspera address
    #$status = $self->{wormbase}->run_command("~gw3/.aspera/connect/bin/ascp -QT -l 300m -i /net/nas17b/vol_homes/homes/gw3/.aspera/connect/etc/asperaweb_id_dsa.putty $address .", $log);
    #if ($status != 0) {$log->log_and_die("Aspera fetch of fastq file $file failed for $experiment_accession\n");}

    # split up the FTP address
    my @addresses = split /;/, $ftp;
    foreach my $address (@addresses) {
      my $file = basename($address);

      # if the full URL is not gven, make a plausible one to use
      if ($address !~ /ftp.sra.ebi.ac.uk/) {
	my ($dirbit) = ($file =~ /^(\S\S\S\d\d\d)/);
	my ($base) = ($file =~ /(\S+)[\.\_]/);
	$address = "ftp.sra.ebi.ac.uk/vol1/fastq/$dirbit/$base/$file";
      }

      # use FTP to fetch the file
      $status = $self->{wormbase}->run_command("/sw/arch/bin/wget -q $address", $log);
      if ($status == 0) {
	if (-s $file) {
	  $count++;
	  $status = $self->{wormbase}->run_command("gunzip -f $file", $log);
	  if ($status != 0) {$log->log_and_die("gunzip of fastq file $file failed for $experiment_accession\n");}
	}
      } else {
	$log->write_to("FTP fetch of fastq file $file failed for $experiment_accession:\n/sw/arch/bin/wget -q $address\n");
	
	# Attempt to get it from the NCBI
      	$log->write_to("Trying FTP from NCBI...\n");
	my ($dirbit) = ($run_accession =~ /^(\S\S\S\d\d\d)/);
	$address = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/$dirbit/${run_accession}/${run_accession}.sra";
	$status = $self->{wormbase}->run_command("/sw/arch/bin/wget -q $address", $log);
	if ($status != 0) {$log->log_and_die("FTP fetch of fastq file $file from NCBI failed for ${experiment_accession}:\n/sw/arch/bin/wget -q $address\n");}
	my $cmd = "$ENV{WORM_PACKAGES}/sratoolkit/bin/fastq-dump --origfmt ${run_accession}.sra --outdir .";
	my $status = $self->{wormbase}->run_command($cmd, $log);
	if ($status != 0) {$log->log_and_die("unpack of .sra file ${run_accession}.sra from NCBI failed for ${experiment_accession}:\n/sw/arch/bin/wget -q $address\n");}
	if (-s "${run_accession}.fastq") {
	  $status = $self->{wormbase}->run_command("rm -rf ${run_accession}.sra", $log);
	  $count++;
	}
      }
    }
  }
  close(ENTRY);

  if ($count) {
    $status = $self->{wormbase}->run_command("touch fastq.done", $log); # set flag to indicate we have finished this
    if ($status != 0) {$log->log_and_die("'touch fastq.done' failed for $experiment_accession\n");}
  } else {
    $log->log_and_die("No fastq files were found for $experiment_accession\n");
  }

  return $failed;

# the following is for downloading the files using Aspera from the NCBI
# This doesnt work on the EBI LSF farm
#  # get the name of the subdirectory in the SRA
#  my ($dirbit) = ($experiment_accession =~ /SRX(\d\d\d)/);
#
#
#  # the Aspera file transfer system is error prone and (in Sanger) only works when run on the farm2-login head node
#  $log->write_to("Get the SRA files for $experiment_accession\n");
#  my $failed = 1;
#  while ($failed) {
#    my $cmd = "~gw3/.aspera/connect/bin/ascp -i /net/nas17b/vol_homes/homes/gw3/.aspera/connect/etc/asperaweb_id_dsa.openssh -k 1 -T -l200m anonftp\@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/SRR150/SRR1503645/SRR1503645.sra  /."
# used to be this but doesn't seem to work any more    my $cmd = "~gw3/.aspera/connect/bin/ascp -k 1 -l 300M -QTr -i /net/nas17b/vol_homes/homes/gw3/.aspera/connect/etc/asperaweb_id_dsa.putty anonftp\@ftp-private.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByExp/litesra/SRX/SRX${dirbit}/$experiment_accession ./";
#    print "cmd = $cmd\n";
#    my $return_status = system($cmd);
#    if ( ( $return_status >> 8 ) == 0 )  {
#       $failed=0;
#       $cmd = "$WORM_PACKAGES/sratoolkit/bin/fastq-dump --origfmt SRR304976.sra --outdir /tmp"
#       my $return_status = system($cmd);
#    }
#  }
}


=head2 

    Title   :   trimmomatic
    Usage   :   $self->trimmomatic($experimenta_accession)
    Function:   trims adapter and poor sequence regions from the reads using Trinity's Trimmomatic program
    Returns :   success = 0
    Args    :   experiment_accession - string

=cut

sub trimmomatic {
  my ($self, $experiment_accession) = @_;

  my $log = $self->{log};
  my $status;
  my $RNASeqSRADir = $self->{RNASeqSRADir};
  my $trimmomatic_done_file = "$RNASeqSRADir/$experiment_accession/SRR/trimmomatic.done";
  chdir "$RNASeqSRADir/$experiment_accession/SRR";
  
  my @files = glob("*.fastq");
  $log->write_to("Have fastq files: @files\n");
  
  # do we have paired reads?
  my $have_paired_reads = 0;
  my @files1 = sort glob("*_1.fastq"); # sort to ensure the two sets of files are in the same order
  my @files2 = sort glob("*_2.fastq");
  if ((@files1 == @files2) && @files2 > 0) {
    $log->write_to("Have paired-end files.\n");
    $have_paired_reads = 1;
    @files = @files1;
  }
  
  my @SRR_ids; # the SRR IDs processed
  
  # run alignment on each fastq file (or paired-read files) in turn
  foreach my $fastq_file (@files) {

    # get the SRR ID of the file
    my $srr_name = $fastq_file;
    $srr_name =~ s/.fastq$//;
    $srr_name =~ s/_1$//;
    push @SRR_ids, $srr_name;
    
    my $cmd='';
    if ($have_paired_reads) {
      my $next_paired_file = shift @files2;
      #$fastq_file .= " $next_paired_file"; # space character between the paired file names
      my $out_paired_one = $fastq_file;
      $out_paired_one =~ s/\.fastq/\.P\.qtrim/;
      my $out_paired_two = $next_paired_file;
      $out_paired_two =~ s/\.fastq/\.P\.qtrim/;
      my $out_unpaired_one = $fastq_file;
      $out_unpaired_one =~ s/\.fastq/\.U\.qtrim/;
      my $out_unpaired_two = $next_paired_file;
      $out_unpaired_two =~ s/\.fastq/\.U\.qtrim/;

      # For paired-ended data:
      $cmd = "$ENV{WORM_PACKAGES}/java7/bin/java -jar $ENV{WORM_PACKAGES}/trinity/trinity-plugins/Trimmomatic/trimmomatic.jar PE -threads 2 -phred33 $fastq_file $next_paired_file $out_paired_one $out_unpaired_one $out_paired_two $out_unpaired_two ILLUMINACLIP:$ENV{WORM_PACKAGES}/trinity/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25";

    } else {
      my $output_file = $fastq_file;
      $output_file =~ s/\.fastq/\.P\.qtrim/;

      # For single-ended data
      $cmd = "$ENV{WORM_PACKAGES}/java7/bin/java -jar $ENV{WORM_PACKAGES}/trinity/trinity-plugins/Trimmomatic/trimmomatic-0.30.jar SE -phred33 $fastq_file $output_file ILLUMINACLIP:$ENV{WORM_PACKAGES}/trinity/trinity-plugins/Trimmomatic/adapters/TruSeq3-SE:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:5 MINLEN:25";
    }

    $log->write_to("Trimmomatic run $cmd\n");
    
    $status = $self->{wormbase}->run_command($cmd, $log);
    if ($status != 0) { 
      $log->write_to("Didn't run Trimmomatic successfully\n");
    } else {
      $status = $self->{wormbase}->run_command("touch $trimmomatic_done_file", $log);
    }
  }

}


=head2 

    Title   :   sam_to_bam
    Usage   :   $self->sam_to_bam($samfile, $bamfile)
    Function:   converts a SAM file to a BAM file and deletes the SAM file
    Returns :   success = 0
    Args    :   SAM file name, BAM file name

=cut

sub sam_to_bam {
  my ($self, $samfile, $bamfile) = @_;
  
  my $failed = 0;
  my $result = system("rm -f $bamfile");
  $result = system($self->{Software}."/samtools/samtools view -Sbh $samfile > $bamfile");
  if ( ( $result >> 8 ) != 0 )  {$failed=1} 
  $result = system("rm -f $samfile");
  return $failed;
}

=head2 

    Title   :   sort_bam_file
    Usage   :   $self->sort_bam_file($bamfile);
    Function:   sorts a BAM file
    Returns :   success = 0
    Args    :   BAM file name

=cut

sub sort_bam_file {
  my ($self, $bamfile) = @_;
  my $log = $self->{log};
  
  my $failed = 0;
  my $status = $self->{wormbase}->run_command($self->{Software}."/samtools/samtools sort -m 3500000000 $bamfile -o $bamfile.sorted", $log); # 3.5Gb
  $failed = $status;
  $status = $self->{wormbase}->run_command("mv $bamfile.sorted $bamfile", $log);
  return ($failed & $status);
}

=head2 

    Title   :   merge_bam_files
    Usage   :   $self->merge_bam_files(\@list_of_files, $output_filename)
    Function:   merges BAM files
    Returns :   
    Args    :   ref to list of Run directories, output BAM file name

=cut

sub merge_bam_files {
  my ($self, $bamfiles, $output_filename) = @_;

  my @SRR_ids = @{$bamfiles};
  my $log = $self->{log};

  if (scalar @SRR_ids > 1) {
    my $bam_files =  join ' ', map {$_ . "/Aligned.out.bam"} @SRR_ids; # join into one string
    $self->{wormbase}->run_command($self->{Software}."/samtools/samtools merge $output_filename $bam_files", $log);
  } else {
    # have just the one bam file, so simply link to it
    my $file = $SRR_ids[0]."/Aligned.out.bam";
    $self->{wormbase}->run_command("ln -s $file $output_filename", $log);
  }
}

=head2 

    Title   :   index_bam_file
    Usage   :   $self->index_bam_file($name_of_file)
    Function:   indexes a BAM file
    Returns :   
    Args    :   name of the BAM file

=cut

sub index_bam_file {
  my ($self, $filename) = @_;

  my $log = $self->{log};

  $self->{wormbase}->run_command($self->{Software}."/samtools/samtools index $filename", $log);
}

=head2 

    Title   :   merge_splice_junction_files
    Usage   :   $self->merge_splice_junction_files(\@SRR_ids, $output_filename);
    Function:   merges splice junction files
    Returns :   
    Args    :   ref to list of Run file directories, output filename

=cut

sub  merge_splice_junction_files {
  my ($self, $bamfiles, $output_filename) = @_;

  my @SRR_ids = @{$bamfiles};

  my %gff; # gff hash holding junction sites and read numbers

  foreach my $srr_name (@SRR_ids) {
    $self->read_junction_file(\%gff, "${srr_name}/SJ.out.tab"); # add junctions to %gff
  }
  $self->write_junction_gff_file(\%gff, $output_filename);
}

=head2 

    Title   :   read_junction_file
    Usage   :   $self->read_junction_file(\%gff, $junctions_file)
    Function:   reads a splice junction files into a hash structure
    Returns :   
    Args    :   ref to hash to receive data, name of junction file

=cut

sub  read_junction_file {
  my ($self, $gff_href, $junction_file) = @_;
  my $log = $self->{log};

  # SJ.out.tab - high confidence collapsed splice junctions in
  # tab-delimited format. Only junctions supported by uniquely mapping
  # reads are reported.
  
  # Column 1: chromosome
  # Column 2: first base of the intron (1-based)
  # Column 3: last base of the intron (1-based)
  # Column 4: strand
  # Column 5: intron motif: 0: non-canonical; 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
  # Column 6: 0: unannotated, 1: annotated (only if splice junctions database is used)
  # Column 7: number of uniquely mapping reads crossing the junction
  # Column 8: number of multi-mapping reads crossing the junction
  # Column 9: maximum spliced alignment overhang

  # the SJ.out.tab file looks like:
  #I       15161   16472   1       1       0       59      0       25
  #I       17959   18005   2       2       0       2       0       12
  #I       19242   20270   2       2       0       9       0       24
  #I       111066  111160  1       1       0       40      0       22
  #I       111261  111509  1       1       0       549     0       25
  #I       111708  111754  1       1       0       712     0       25
  #I       111972  112018  1       1       0       2104    0       25


  open(JUNCT, "<$junction_file") || $log->log_and_die("Can't open the file $junction_file\n");
  while (my $line = <JUNCT>) {
    if ($line =~ /^track/) {next}
    my @cols = split /\s+/, $line;
    my $chrom = $cols[0];
    my $start = $cols[1];
    my $end = $cols[2];
    my $sense = $cols[3];
    my $reads = $cols[6] + $cols[7]; # add the unique and multimapped reads

    $chrom = $self->{wormbase}->chromosome_prefix . $chrom;

    $gff_href->{$chrom}{$start}{$end}{$sense} += $reads;

  }
  close(JUNCT);
 
}

=head2 

    Title   :   write_junction_gff_file
    Usage   :   $self->write_junction_gff_file(\%gff, $junctions_file)
    Function:   writes a splice junction hash into a GFF file
    Returns :   
    Args    :   ref to hash to receive data, name of junction GFF file

=cut

sub  write_junction_gff_file {
  my ($self, $gff_href, $junction_file) = @_;
  my $log = $self->{log};


  my $strand;
  open (GFF, ">$junction_file") || $log->log_and_die("Can't open $junction_file\n");
  foreach my $chrom (sort {$a cmp $b} keys %{$gff_href} ) {
    foreach my $start (sort {$a <=> $b} keys %{$gff_href->{$chrom}}) {
      foreach my $end (sort {$a <=> $b} keys %{$gff_href->{$chrom}{$start}}) {
	foreach my $sense (sort {$a <=> $b} keys %{$gff_href->{$chrom}{$start}{$end}}) {
	  my $reads = $gff_href->{$chrom}{$start}{$end}{$sense};
	  if ($sense == 1) {$strand='+'} elsif ($sense == 2) {$strand='-'} else {$strand = '.'} # sense=0 is used by STAR when there is a non-canonical splice
	  print GFF "$chrom\tintron\tintron\t$start\t$end\t.\t$strand\t$reads\n";
	}
      }
    }
  }
  close (GFF);

}

=head2

    Title   :   run_make_gtf_transcript
    Usage   :   $self->run_make_gtf_transcript($database)
    Function:   runs make_gtf_transcript
    Returns :   success = 0
    Args    :   path of acedb database to query (e.g. autoace)

=cut

sub run_make_gtf_transcript {
  my ($self, $database) = @_;

  my $RNASeqSRADir = $self->{RNASeqSRADir};
  my $RNASeqGenomeDir = $self->{RNASeqGenomeDir};
  my $gtf_file = "$RNASeqGenomeDir/transcripts.gtf";
  my $species = $self->{wormbase}->{species};
  my $log = $self->{log};
  my $status;

  if ($self->new_genome) {system("rm -f $gtf_file")}
  if ($self->check && -s $gtf_file) {
    $log->write_to("GTF file already exists - not making them again\n");
    $status = 0;
  } else {

    my $gtf_file = " $RNASeqGenomeDir/transcripts.gtf";
    my $cmd = "make_GTF_transcript.pl -database $database -out $gtf_file -species $species -noprefix"; # -noprefix to strip the chromosome prefix out
    
    $log->write_to("run $cmd\n");
    
    $status = $self->{wormbase}->run_script($cmd, $log);
    if ($status != 0) { $log->log_and_die("Didn't run make_GTF_transcript.pl successfully\n"); }
    sleep(60); # wait a minute to let the NFS system catch up with all the LSF nodes
    #if (!-e $gtf_file || -z $gtf_file) { $log->log_and_die("Didn't run make_GTF_transcript.pl successfully\n"); }
  }
  return $status;
}

=head2

    Title   :   run_cufflinks
    Usage   :   $self->run_cufflinks($experiment_accession);
    Function:   runs cufflinks on the accepted_hits.bam file
    Returns :   success = 0
    Args    :   experiment accession, strandedness

=cut

sub run_cufflinks {
  my ($self, $experiment_accession, $strandedness, $library_type) = @_;
  
  my $RNASeqSRADir = $self->{RNASeqSRADir};
  my $RNASeqGenomeDir = $self->{RNASeqGenomeDir};
  my $log = $self->{log};
  my $done_file = "$RNASeqSRADir/$experiment_accession/cufflinks/genes.fpkm_tracking.done"; # made when cufflinks is run successfully
  my $Software = $self->{Software};
  my $alignmentDir = $self->{'alignmentDir'};
  my $status;

  $log->write_to("run cufflinks\n");
  chdir "$RNASeqSRADir/$experiment_accession";
  mkdir "cufflinks", 0777;
  chdir "cufflinks";

  if ($self->{check} && -e $done_file) {
    $log->write_to("Cufflinks files already exist - not making them again\n");
    $status = 0;
  } else {

    # For strand-specific data, you do not need any extra parameters for
    # STAR runs, but you need to use --library-type option for
    # Cufflinks. For example, for the "standard" dUTP protocol you need to
    # use --library-type fr-firststrand in Cufflinks.
    my $strand_option='';
    #my $strand_option='--library-type fr-firststrand';
    # default is to assume not stranded - this is correct for the Hillier modENCODE data
    if (defined $strandedness && $strandedness eq 'stranded' && defined $library_type && $library_type eq 'fr') {
      $strand_option = '--library-type fr-firststrand'
    }

    my $gtf = "--GTF $RNASeqGenomeDir/transcripts.gtf";
    $status = $self->{wormbase}->run_command("$Software/cufflinks/cufflinks $gtf $strand_option ../$alignmentDir/accepted_hits.bam", $log);
    if ($status != 0) {  $log->log_and_die("Didn't run cufflinks to get the isoform/gene expression successfully\n"); }
    if (-s "genes.fpkm_tracking" > 10000 && -s "isoforms.fpkm_tracking" > 10000) {
      $self->{wormbase}->run_command("touch $done_file", $log); # set flag to indicate we have finished this
    } else {
      $log->log_and_die("Didn't run cufflinks to get the isoform/gene expression successfully\n");
    }
  }

  # make the cufflinks assembly file for Michael Paulini to model into operons :-)
  my $species = $self->{wormbase}->{species};
  if ($species ne 'elegans' && !-e "assembly/transcripts.gtf") {
    mkdir "assembly", 0777;
    chdir "assembly";
    my $strand_option = "--min-intron-length 25 --max-intron-length 30000";
    $status = $self->{wormbase}->run_command("$Software/cufflinks/cufflinks $strand_option ../../$alignmentDir/accepted_hits.bam", $log);
    if ($status != 0) {  $log->log_and_die("Didn't run cufflinks to get the cufflinks gene structures successfully\n"); }    
  }

  return $status;
}


  

=head2

    Title   :   get_introns
    Usage   :   $self->get_introns($experiment_accession);
    Function:   gets the introns from the alignment output
    Returns :   success = 0
    Args    :   experiment accession

=cut

sub get_introns {
  my ($self, $experiment_accession) = @_;

  my $RNASeqSRADir = $self->{RNASeqSRADir};
  my $RNASeqGenomeDir = $self->{RNASeqGenomeDir};
  my $log = $self->{log};
  my $done_file = "$RNASeqSRADir/$experiment_accession/Introns/Intron.ace.done";
  my $Software = $self->{Software};
  my $alignmentDir = $self->{'alignmentDir'};
  my $status;
  my $database = $self->{wormbase}->{autoace};

  $log->write_to("get introns\n");
  chdir "$RNASeqSRADir/$experiment_accession";
  mkdir "Introns", 0777;
  chdir "Introns";
  

  my %seqlength;
  
  my $coords = Coords_converter->invoke($database, 0, $self->{wormbase});
  
  my $output = "Intron.ace";
  my $junctions = "../$alignmentDir/junctions.gff";
  my $virtual;
  my $method = "RNASeq_splice";
  my $text = "\"RNASeq intron\"";


  open (ACE, ">$output") || $log->log_and_die("Can't open the file $output\n");
  
  open(GFF, "<$junctions") || $log->log_and_die("Can't open the file $junctions\n");

  my $vfile = "virtual_objects." . $self->{wormbase}->{species} . ".RNASeq.ace";
  open(VIRT, ">$vfile") or $log->log_and_die("Could not open $vfile for writing\n");

  my (@tiles, @whole_chromosome);


  my $sequence = '';
  my $sequence_len = 0;

  while (my $line = <GFF>) {

    if ($line =~ /^#/) {next}
    my @cols = split /\s+/, $line;
    my $chrom = $cols[0];
    my $start = $cols[3];
    my $end   = $cols[4];
    my $sense = $cols[6];
    my $reads = $cols[7];
    
    if ($sense eq '-') {
      ($end, $start) = ($start, $end);
    }
    
    if ($chrom ne $sequence) { # new sequence
 
      $self->write_tiles(\@tiles, \@whole_chromosome, $virtual, $sequence, $sequence_len); #  write the old data


      $sequence = $chrom;
      $virtual = "${sequence}:Confirmed_intron_RNASeq";
      $sequence_len = $self->initialise_tiles($sequence, \@tiles, $coords);

    }

    $self->store_feature_in_tile(\@tiles, \@whole_chromosome, $method, $start, $end, $reads, $text);

  }

  # write the last sequence
  $self->write_tiles(\@tiles, \@whole_chromosome, $virtual, $sequence, $sequence_len); #  write the old data


  close(GFF);
  close(ACE);
  close(VIRT);
  
  $self->{wormbase}->run_command("touch $done_file", $log); # set flag to indicate we have finished this
  
  return $status;
}



##########################################
# find the tile to store the Feature in
sub store_feature_in_tile {
  my ($self, $tiles_aref, $whole_chromosome_aref, $method, $start, $end, $reads, $text) = @_;

  my $found = 0;
  for( my $tile_idx = 1; $tile_idx <= @{$tiles_aref}; $tile_idx++) {
    my $tile = $tiles_aref->[$tile_idx-1];
    if ($start < $end) {
      if ($start > $tile->{start} && $end <= $tile->{end}) { # find the tile containing this forward Feature
	push @{$tile->{segs}}, [$method, $start - $tile->{start} + 1, $end - $tile->{start} + 1, $reads, $text];
	$found = 1;
      }
    } else {
      if ($end > $tile->{start} && $start <= $tile->{end}) { # find the tile containing this reverse Feature
	push @{$tile->{segs}}, [$method, $start - $tile->{start} + 1, $end - $tile->{start} + 1, $reads, $text];
	$found = 1;
      }
    }
  }
  if (!$found) { # it falls between two tiles, so place it on the top-level Sequence
    push @{$whole_chromosome_aref}, [$method, $start, $end, $reads, $text];
  }
}

##########################################
sub write_tiles {
  my ($self, $tiles_aref, $whole_chromosome_aref, $virtual, $sequence, $sequence_len) = @_;


  # output the new Sequence lines

  my @sequence_out;
  my @feature_out;

  if (scalar @{$tiles_aref}) {
    push @sequence_out, "\nSequence : \"${sequence}\"\n";

    for(my $tile_idx = 1; $tile_idx <= @{$tiles_aref}; $tile_idx++) {
      my $tile = $tiles_aref->[$tile_idx-1];
      
      my $vseq = "${virtual}:$tile_idx";
      
      if (@{$tile->{segs}}) {
	push @sequence_out, "S_Child Feature_data ". $vseq ." ". $tile->{start} ." ". $tile->{end} ."\n";
	
	push @feature_out, "\nFeature_data : \"$vseq\"\n";
	foreach my $seg (@{$tile->{segs}}) {
	  push @feature_out, "Feature @$seg\n";
	}
      }
    }
  }


  if (scalar @{$whole_chromosome_aref}) {
    push @sequence_out, "\nSequence : \"${sequence}\"\n";
    push @sequence_out, "S_Child Feature_data ${virtual} 1 $sequence_len\n";
    push @feature_out, "\nFeature_data : ${virtual}\n";
    foreach my $seg (@{$whole_chromosome_aref}) {
      push @feature_out,  "Feature @$seg\n";
    }
  }

  print VIRT @sequence_out;
  print VIRT "\n"; # acezip.pl concatenates another line to the last line if this is not blank

  print ACE @feature_out;
  print ACE "\n";

  @{$tiles_aref} = ();
  @{$whole_chromosome_aref} = ();
}


##########################################
sub initialise_tiles {
  my ($self, $sequence, $tiles_aref, $coords) = @_;

  my $log = $self->{log};

  my $chr_len = $coords->Superlink_length($sequence);
  if (!defined $chr_len) {$log->log_and_die("Can't find the length of the Sequence $sequence\n")}

  for(my $i=0; $i < $chr_len; $i += 300000) {
    my $chr_start = $i + 1;
    my $chr_end = $chr_start + 300000 - 1;
    $chr_end = $chr_len if $chr_end > $chr_len;
    push @{$tiles_aref}, {
		  start => $chr_start, 
		  end   => $chr_end,
		  segs  => [],
    }
  }
  return $chr_len;
}
##########################################

=head2

    Title   :   check_all_done
    Usage   :   $self->check_all_done($experiment_accession);
    Function:   checks that all of the 'done' files exist for an experiment
    Returns :   1 if all the files exist
    Args    :   experiment accession, not in build int (so don't want to run cufflinks)

=cut

sub check_all_done {
  my ($self, $experiment_accession, $notbuild) = @_;

  my $RNASeqSRADir = $self->{RNASeqSRADir};
  my $alignmentDir = $self->{'alignmentDir'};
  my $intron_done_file    = "$RNASeqSRADir/$experiment_accession/Introns/Intron.ace.done"; 
  my $bam_done_file       = "$RNASeqSRADir/$experiment_accession/".$self->{'alignmentDir'}."/accepted_hits.bam.done";
  my $hits_file           = "$RNASeqSRADir/$experiment_accession/".$self->{'alignmentDir'}."/hits.tmp";
  my $stranded_hits_file  = "$RNASeqSRADir/$experiment_accession/".$self->{'alignmentDir'}."/hits.stranded";
  my $cufflinks_done_file = "$RNASeqSRADir/$experiment_accession/cufflinks/genes.fpkm_tracking.done";
  my $status=0;
  if (-e $intron_done_file &&
      -e $bam_done_file &&
      -e $hits_file &&
      -e $stranded_hits_file &&
      (-e $cufflinks_done_file || $notbuild)) {
    $status=1;
  }

  return $status;
}

=head2

    Title   :   make_hits
    Usage   :   $self->make_hits($experiment_accession);
    Function:   make the hits file giving the count of hits at each base and the sense of the alignment (NB not the sense of the transcript)
    Returns :   
    Args    :   experiment accession

=cut

sub make_hits {
  my ($self, $experiment_accession) = @_;
  
  my $RNASeqSRADir = $self->{RNASeqSRADir};
  my $alignmentDir = $self->{'alignmentDir'};
  my $exptDir = "$RNASeqSRADir/$experiment_accession/$alignmentDir";
  my $bamfile = "$exptDir/accepted_hits.bam";
  my $hitsfile = "$exptDir/hits.tmp";
  my $Software = $self->{Software};
  
  if (-e $bamfile) {
    # get the hits with a count
    print "writing $hitsfile\n";
    unlink "$hitsfile";
    system(qq#/$Software/BEDTools/bin/bamToBed -split -i $bamfile | awk '{OFS=\"\t\"; print \$1,\$2,\$3,\$6 }' | sort -k1,1 -k2,3n | uniq -c > $hitsfile#);
  }
}

=head2

    Title   :   make_stranded_hits
    Usage   :   $self->make_stranded_hits($experiment_accession);
    Function:   make the hits file giving the count of hits at each base and the sense of the transcript that they came from
    Returns :   
    Args    :   experiment hashref

=cut
  
sub make_stranded_hits {
  my ($self, $experiment) = @_;

  my $experiment_accession = $experiment->{experiment_accession};
  my $strandedness = $experiment->{strandedness};
  my $library_layout = $experiment->{library_layout};
  my $library_type = $experiment->{library_type};

  my $log = $self->{log};
  my $RNASeqSRADir = $self->{RNASeqSRADir};
  my $alignmentDir = $self->{'alignmentDir'};
  my $exptDir = "$RNASeqSRADir/$experiment_accession/$alignmentDir";
  my $bamfile = "$exptDir/accepted_hits.bam";
  my $stranded_hitsfile = "$exptDir/hits.stranded";
  my $Software = $self->{Software};

  $self->{wormbase}->run_command("rm -f $stranded_hitsfile", $log);

  # just make an empty hits file if the reads are not stranded
  if ($strandedness eq 'unknown' || $strandedness eq 'unstranded') {
    $self->{wormbase}->run_command("touch $stranded_hitsfile", $log);
  } else {


    if (-e $bamfile) {
      print "writing $stranded_hitsfile\n";
      my $bamtobed = "$Software/BEDTools/bin/bamToBed -split -i $bamfile";
      open (HITS, ">${stranded_hitsfile}.tmp") || $log->log_and_die("Can't open ${stranded_hitsfile}.tmp in make_stranded_hits()\n");
      open (BED, "$bamtobed |") || $log->log_and_die("Can't run $bamtobed in make_stranded_hits()\n");
      while (my $line = <BED>) {
	#      IV      4248069 4248094 SRR548309.19790789/1    255     +
	#      IV      4247931 4247972 SRR548309.15239765/2    255     -
	
	my ($bed_chrom, $bed_start, $bed_end, $bed_name, $bed_sense) = ($line =~ /^(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+\d+\s+(\S+)/);

	if ($library_layout eq 'PAIRED' && $library_type ne 'unknown') {
	  my $mate='';
	  if ($bed_name =~ /\/1$/) {$mate='first'}
	  if ($bed_name =~ /\/2$/) {$mate='second'}

	  if ($strandedness eq 'stranded') {
	    if ($library_type eq 'fr') {
	      if ($mate eq 'second') {$bed_sense = ($bed_sense eq '+') ? '-' : '+' }
	      
	    } elsif ($library_type eq 'ff') { 
	      # no action - both mates are already in forward sense of transcript
	    }
	    
	  } elsif ($strandedness eq 'reverse_stranded') { # not sure if this occurs in our data, but deal with it if it does
	    if ($library_type eq 'rr') {
	      if ($mate eq 'first')  {$bed_sense = ($bed_sense eq '+') ? '-' : '+' }
	      if ($mate eq 'second') {$bed_sense = ($bed_sense eq '+') ? '-' : '+' }
	      
	    } elsif ($library_type eq 'rf') {
	      if ($mate eq 'first') {$bed_sense = ($bed_sense eq '+') ? '-' : '+' }
	    }
	  }
	  print HITS "$bed_chrom\t$bed_start\t$bed_end\t$bed_sense\n";
	  
	} elsif ($library_layout eq 'SINGLE') {
	  if ($strandedness eq 'reverse_stranded') { # e.g. if it is ABI SOLID
	    $bed_sense = ($bed_sense eq '+') ? '-' : '+';
	  }
	  print HITS "$bed_chrom\t$bed_start\t$bed_end\t$bed_sense\n";
	}
      }
      close(BED);
      close(HITS);
      
      # now sort the file
      $self->{wormbase}->run_command("sort -k1,1 -k2,3n ${stranded_hitsfile}.tmp | uniq -c > $stranded_hitsfile", $log);
      unlink "${stranded_hitsfile}.tmp";
    }
  }


}
#####################################################################################################
# set up a genome for Tophat

=head2 

    Title   :   setup_genome_for_tophat
    Usage   :   $self->setup_genome_for_tophat()
    Function:   sets up and indexes a fresh copy of the genome for aligning with Tophat
    Returns :   
    Args    :   masked - use a masked genome if true

=cut

sub setup_genome_for_tophat {
  my ($self, $masked) = @_;
 
}


#####################################################################################################
# set up a genome and run alignments for STAR

=head2 

    Title   :   setup_genome_for_star
    Usage   :   self->setup_genome_for_star($masked)
    Function:   sets up and indexes a fresh copy of the genome for aligning with STAR
    Returns :   
    Args    :   

=cut

sub setup_genome_for_star {
  my ($self) = @_;
  
  my $species = $self->{'wormbase'}->species;
  my $RNASeqGenomeDir = $self->{RNASeqGenomeDir};
  my $database = $self->{wormbase}->autoace;
  my $chrom_file = "$species.genome.fa"; # not the masked file
  my $species_genome = "$RNASeqGenomeDir/STAR/${species}.fa";
  my $done_file = "$RNASeqGenomeDir/STAR/${species}.fa.done";
  my $status;
  my $log = $self->{log};

  # only create the genome sequence index if the assembly changed or the genome has not been set up
  if ($self->{new_genome} || !-e $done_file) {
    mkdir $RNASeqGenomeDir, 0777;
    mkdir $RNASeqGenomeDir."/STAR", 0777;
    mkdir $RNASeqGenomeDir."/STAR/reference-indexes", 0777; # this looks bogus/unused - remove?
    chdir $RNASeqGenomeDir."/STAR";
    unlink $done_file;
    unlink glob("chr*");
    unlink glob("*.txt");
    unlink glob("*.fa");
    unlink ("Genome");
    my $source_file = "${database}/SEQUENCES/${chrom_file}";
    if (-e $source_file) {
      my $copy_cmd = "cp $source_file .";
      $status = $self->{'wormbase'}->run_command($copy_cmd, $log);
      
    } elsif (-e "${source_file}.gz") {
      my $gzip_cmd = "gunzip -c ${source_file}.gz > ./$chrom_file";
      $status = $self->{'wormbase'}->run_command($gzip_cmd, $log);
      
    } else {
      $log->log_and_die("RNASeq: Can't locate the genome file: $source_file\n");
    }
    
    $self->strip_prefix_from_fasta_file($chrom_file, $species_genome);
    

    # run in four threads - takes 5 mins
    # uses 3.3 Gb memory for elegans
    # uses 5.7 Gb memory for remanei

    # For parallel jobs, EBI recommend you use 4 CPUs by default, as
    # this will give you access to the most number of machines. Tell
    # LSF how many CPUs you want with -n, and tell the software the
    # same number.
    my $star_cmd = $self->{Software}."/star/STAR --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles ${species_genome} --runThreadN 4";
      $status = $self->{'wormbase'}->run_command($star_cmd, $log);
    if ($status != 0) {  $log->log_and_die("RNASeq: Didn't create the STAR indexes $RNASeqGenomeDir/STAR\n"); }
    $self->{wormbase}->run_command("touch $done_file", $log); # set flag to indicate we have finished this

  } else {
    $log->write_to("Genome database files already exist - not making them again\n");
    $status = 0;
  }
}



=head2 

    Title   :   align_star
    Usage   :   self->align_star($experiment_accession)
    Function:   gets the SRR fastq files, runs a STAR session to align the an experiment then runs cufflinks
    Returns :   
    Args    :   experiment hashref, number of threads for STAR to use int, not in build int (so don't run cufflinks)

=cut

sub align_star {
  my ($self, $experiment, $threads, $notbuild) = @_;

  my $log = $self->{log};
  my $RNASeqSRADir = $self->{RNASeqSRADir};
  my $experiment_accession = $experiment->{experiment_accession};
  my $bam_done_file = "$RNASeqSRADir/$experiment_accession/".$self->{'alignmentDir'}."/accepted_hits.bam.done";
  my $srr_done_file = "$RNASeqSRADir/$experiment_accession/SRR/fastq.done";
  my $trimmomatic_done_file = "$RNASeqSRADir/$experiment_accession/SRR/trimmomatic.done";
  my $solid_done_file = "$RNASeqSRADir/$experiment_accession/SRR/solid.done";
  my $alignmentDir = $self->{'alignmentDir'};
  my $exptDir = "$RNASeqSRADir/$experiment_accession/$alignmentDir";
  my $bam_hits = "$exptDir/hits.tmp";
  my $bam_stranded_hits = "$exptDir/hits.stranded";
  my $introns_done_file = "$RNASeqSRADir/$experiment_accession/Introns/Intron.ace.done";

  if ($self->{new_genome} || !-e $bam_done_file || !-e $srr_done_file) { 
#  if ($self->{new_genome} || !-e $bam_done_file || !-e $srr_done_file || !-e $trimmomatic_done_file) { 

    if ($self->{new_genome} || !-e $srr_done_file) { 
      $self->get_SRA_files($experiment_accession);
    }

    # convert any ABI_SOLID files to fastq
    if ($experiment->{instrument_platform} eq 'ABI_SOLID' && !-e $solid_done_file) {
      chdir "$RNASeqSRADir/$experiment_accession/SRR";
      my @files = glob("*.fastq");
      foreach my $file (@files) {
	my $cmd = "mv $file $file.solid";
	$self->{wormbase}->run_command($cmd, $log);
	$cmd = "perl ".$self->{Software}."/SOLID2std/SOLID2std.pl -fastq $file.solid -o $file";
	$self->{wormbase}->run_command($cmd, $log);
	$cmd = "rm -f $file.solid";
	$self->{wormbase}->run_command($cmd, $log);
      }
      my $status = $self->{wormbase}->run_command("touch $solid_done_file", $log); # set flag to indicate we have finished this
    }

#    if (!-e $trimmomatic_done_file) {$self->trimmomatic($experiment_accession)}

    if ($self->{new_genome} || !-e $bam_done_file) {
      chdir "$RNASeqSRADir/$experiment_accession"; # get into the correct experiment directory
      $self->{wormbase}->run_command("rm -rf ".$self->{'alignmentDir'}, $log); #delete the old tophat_out or star_out
      $self->run_star_alignment($experiment, $threads);
    }
  }

  if (!-e $bam_hits)              {$self->make_hits($experiment_accession)}
  if (!-e $bam_stranded_hits)     {$self->make_stranded_hits($experiment)}
  if (!-e $introns_done_file)     {$self->get_introns($experiment_accession)}

  if (! $notbuild) {
    my $strandedness = $experiment->{strandedness};
    my $library_type = $experiment->{library_type};
    $self->run_cufflinks($experiment_accession, $strandedness, $library_type);
  }
}

=head2 

    Title   :   run_star_alignment
    Usage   :   self->run_star_alignment
    Function:   cd to the experiment directory, identify the SRR files to use, runs the STAR alignment command on each SRR file
    Returns :   
    Args    :   experiment hashref, number of threads to use

=cut

sub run_star_alignment {
  my ($self, $experiment, $threads) = @_;
  
  my $status;
  my $log = $self->{log};
  my $RNASeqGenomeDir = $self->{RNASeqGenomeDir};
  my $RNASeqSRADir = $self->{RNASeqSRADir};
  my $genomeDir = $RNASeqGenomeDir."/STAR";
  my $experiment_accession = $experiment->{experiment_accession};
  my $strandedness = $experiment->{strandedness};
  my $alignmentDir = $self->{'alignmentDir'}; # tophat_out or star_out
  
  if (!defined $threads || $threads !~ /^\d+$/) {$threads = 1}

  mkdir "$RNASeqSRADir/$experiment_accession/$alignmentDir", 0777;
  chdir "$RNASeqSRADir/$experiment_accession/$alignmentDir"; # get into the alignment directory
  
  my @files = glob("../SRR/*.fastq");
  $log->write_to("Have fastq files: @files\n");
  
  # do we have paired reads?
  my $have_paired_reads = 0;
  my @files1 = sort glob("../SRR/*_1.fastq"); # sort to ensure the two sets of files are in the same order
  my @files2 = sort glob("../SRR/*_2.fastq");
  if ((@files1 == @files2) && @files2 > 0) {
    $log->write_to("Have paired-end files.\n");
    $have_paired_reads = 1;
    @files = @files1;
  }
  
  my @SRR_ids; # the SRR IDs processed
  
  # run alignment on each fastq file (or paired-read files) in turn
  foreach my $fastq_file (@files) {

    # get the SRR ID of the file
    my $srr_name = $fastq_file;
    $srr_name =~ s/^\.\.\/SRR\///;
    $srr_name =~ s/.fastq$//;
    $srr_name =~ s/_1$//;
    push @SRR_ids, $srr_name;
    
    # as we are now chdir'd a directory level down, we must add '../' to the SRR path
    $fastq_file = '../'.$fastq_file;

    if ($have_paired_reads) {
      my $next_paired_file = shift @files2;
      $next_paired_file = '../'.$next_paired_file; # as we are no chdir'd a directory level dow, we must add '../' to the SRR path
      $fastq_file .= " $next_paired_file"; # space character between the paired file names
    }
    
    # make a subdirectory for each run file's output data
    mkdir $srr_name, 0777;
    chdir $srr_name;

    # useful descriptions of parameters are found in the program's STAR/parametersDefault file

    #/pathToStarDir/STAR --genomeDir /path/to/GenomeDir --readFilesIn /path/
    #to/read1 [/path/to/read2] --runThreadN <n> --<inputParameterName> <input
    #parameter value(s)>
    
    #my $options = "";
    # outReadsUnmapped Fastx   : output in separate fasta/fastq files, Unmapped.out.mate1/2
    # alignIntronMin 25
    # alignIntronMax maximum intron size, if 0, max intron size will be determined by (2^winBinNbits)*winAnchorDistNbins (500Kb)
    # alignMatesGapMax maximum gap between two mates, if 0, max intron gap will be determined by (2^winBinNbits)*winAnchorDistNbins (500Kb)
    # outFilterMultimapNmax 10 int: read alignments will be output only if the read maps fewer than this value, otherwise no alignments will be output (1=Only uniquely mapping reads)
    # outFilterMismatchNoverLmax 0.05 float: alignment will be output only if its ratio of mismatches to mapped length is less than this value
    # chimSegmentMin int>0: minimum length of chimeric segment length, if ==0, no chimeric output
    # chimJunctionOverhangMin int>0: minimum overhang for a chimeric junction
    # outSAMstrandField None string: Cufflinks-like strand field flag None : not used intronMotif : strand derived from the intron motif. Reads with inconsistent and/or non-canonical introns are filtered out.

      # uses 3.2 Gb memory for single-end elegans - 4 threads takes 15 mins (tophat takes 2 hours)
      # uses 4.3 Gb memory for paired-end remanei - 4 threads takes 90 mins - 3670 contigs

    my $options = "--alignIntronMin 25 --outReadsUnmapped Fastx --alignIntronMax 15000 --alignMatesGapMax 50000 --outFilterMultimapNmax 2 --outFilterMismatchNoverLmax 0.02 --chimSegmentMin 15 --chimJunctionOverhangMin 15";

    # For non-strand-specific data, you need to use STAR option
    # --outSAMstrandField intronMotif which will add the XS attribute
    # to all canonically spliced alignments using their introns'
    # motifs - that's exactly what Cufflinks needs.

    # For strand-specific data, you do not need any extra parameters for
    # STAR runs, but you need to use --library-type option for
    # Cufflinks. For example, for the "standard" dUTP protocol you need to
    # use --library-type fr-firststrand in Cufflinks.
    
    # default is to assume not stranded - this is correct for the Hillier modENCODE data
    ##if (!defined $strandedness || $strandedness eq 'unstranded') { 
      $options .= " --outSAMstrandField intronMotif";
    ##}
    
    my $STARcmd = $self->{Software}."/star/STAR --genomeDir $genomeDir --readFilesIn $fastq_file $options --runThreadN $threads";
    $log->write_to("Running STAR command:\n$STARcmd\n\n");
    $status = $self->{wormbase}->run_command($STARcmd, $log);
    if ($status != 0) {$log->log_and_die("Didn't run STAR to do the alignment successfully for $srr_name in $experiment_accession\n");}

    $status = $self->sam_to_bam('Aligned.out.sam', 'Aligned.out.bam');
    if ($status != 0) {$log->log_and_die("sam_to_bam() failed for $srr_name in $experiment_accession\n");}
    $status = $self->sort_bam_file('Aligned.out.bam');
    if ($status != 0) {$log->log_and_die("sort_bam_file() failed for $srr_name in $experiment_accession\n");}
    if (-e glob("Aligned.out.bam.sorted*")) {$log->log_and_die("sort_bam_file() failed for $srr_name in $experiment_accession - a BAM sort file is still existant\n");}
    chdir ".."; # back up out of the $srr_name run directory
  }

  $status = $self->merge_bam_files(\@SRR_ids, 'accepted_hits.bam');
  if ($status != 0) {$log->log_and_die("merge_bam_files() failed for $experiment_accession\n");}
  $status = $self->index_bam_file('accepted_hits.bam');
  if ($status != 0) {$log->log_and_die("index_bam_file() failed for $experiment_accession\n");}
  $self->merge_splice_junction_files(\@SRR_ids, 'junctions.gff');
  $status = $self->{wormbase}->run_command("touch accepted_hits.bam.done", $log); # set flag to indicate we have finished this
  if ($status != 0) {$log->log_and_die("'touch accepted_hits.bam.done' failed for $experiment_accession\n");}
}


1;

__END__


Comments about STAR
-------------------

alexdobin

Long intron "spurious" alignments: as @pbluescript wrote, it's all
about filtering. If you are interested in high confidence junctions
only, I would recommend keeping only canonical GT/AG introns, unique
mappers, and supported by 2 or more reads junctions. Annotated
junctions can also be considered trustworthy. Another useful filtering
option is to remove junctions which are supported by a very short
overhang, STAR is using this filtering to create the collapsed
junction file. If your species does not have long introns, you can
make STAR (newest version) filter them out with --alignIntronMax and
--alignMatesGapMax options.

Q

We are currently using the STAR aligner software (version 2.3.0.1) to
align RNA sequences to a reference genome. The majority of inserts in
our sequencing libraries have short fragment lengths (< 60bp) and
consequently STAR returns an output file stating that we have a large
number of reads that are too short to be mapped to unique regions of
the reference genome. We would like to increase the number of mappable
reads, but are unsure of the default read length threshold used by
STAR and more importantly how we can lower this using the STAR command
line. We were wondering if anyone could assist with these queries?

A Chipper

It is explained in the google group somewhere, you can try
--outFilterScoreMinOverLread 0.25 --outFilterMatchNminOverLread 0.25
to allow 50 bp alignments from 2x100 bp runs. Default is 0.66 which
means alignments < 2*100*0.66 are discarded.

A alexdobin

Another option would be to clip the adapter sequence that must appear
on 3' of you reads if your fragments are so short. You can do it with
an external of trimmer software, or use --clip3pAdapterSeq within
STAR.



Q 

1. I used STAR for mapping RNA seq data with default parameters only
changing --outFilterMismatchNmax 5 (allowing 5 mismatches, default was
10). I am now confused if I should have gone for 3 or less? My read
length varies between 10-50. What is the ideal cutoff one should allow
for mismatch parameter?

2. Also, one should go for multi-mapping reads or not? because around
24% of reads in my case is showing multi-mapping, is it significant? i
don't think so.. should I switch it to zero multi-mapping?

3. I also want to do some statistics to see the correlation between my
replicates, and some post-mapping statistical analysis like length
counts, log(n) etc..can anybody refer a good tool for this..?

A alexdobin

1. The max number of mismatches depends on many factors, the most
important being read length (more mismatches for longer reads),
sequencing quality (more mismatches for poorer quality), sequence
divergence of your sample with respect to the reference.

Since you have reads of widely varying length (are these short RNA?),
I would recommend setting the max number of mismatches relative to the
read length using --outFilterMismatchNoverLmax 0.05 This is what we
typically use for ENCODE small RNA-seq. This means that for reads <20b
no mismatches are allowed, 20-39b: 1 mismatch, 40-59b 2 mismatches and
so on.

2. Multi-mapping reads are quite common especially when you map short
reads. 24% is a decent number for short RNA-seq. What to do with
multi-mappers is a complicated issue. Typically they are used to
represent expression of whole classes/families of RNA (repeats
(e.g. transposons), gene families etc). If you are looking at your
data in a browser, it's useful to have two separate tracks - for
unique mappers and multi-mappers.

3. We typically check correlation between bio-reps for read counts on
annotated small RNAs. A good tool to do it is 'coverageBed' from
BEDtools package.  Another nice tool is HTseq.





Q 

I am using STAR for my single-end small RNAseq data (15-46 read
length) and its working great. thanks! However, since its new, I can't
find much help on public forums, there are couple of things I need to
know, and who can tell it better than you..

I used STAR with the following parameters:
--genomeDir </path/genome/hg19>
--genomeLoad NoSharedMemory
--readFilesIn <path/inputReadsFile/min15_max_46length_reads>
--runThreadN 4
--outFilterMismatchNoverLmax 0.05
--outSAMattributes All
--outSJfilterOverhangMin 20 10 10 10
--outFileNamePrefix <output_location>
--sjdbFileChrStartEnd <introns_ucsc_hg19>
--sjdbOverhang 45


Now, when I created contigs / or visualize it on genome browser I am
getting large chunks of regions spanning 100,000bp and more in the
genome, and there are many such occurrences in multiple files. which
therefore is giving large set of contigs as well. Is this a mapping
error?

Also, when I identified those regions from my bam files, I found
bitwise flag '64' in the second column, and '255' in 5th (for uniquely
mapped).


A alexdobin


I am not sure how you make contigs, but I suspect that these long
clusters are made of alignments that are spliced. Some small RNA can
be spliced (e.g. short mRNA decay products that happen to cross
exon-exon junction). Unless you are looking for them specifically, I
would suggest prohibiting splicing with --alignIntronMax 10.  If you
do not want splicing, you also do not --sjdbFileChrStartEnd
<introns_ucsc_hg19> --sjdbOverhang 45, and in any case, these
parameters should be used at the genome generation step, not mapping
step.

Note, that you can also filter out spliced alignments after the
mapping - they will contain "N" in the CIGAR string.

We are routinely using STAR to map "small RNA" (~<200b) data within
the ENCODE project with the following parameters:
--outFilterMismatchNoverLmax 0.05 --outFilterMatchNmin 16
--outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0
--alignIntronMax 1 (>=16b matched to the genome, number of mismatches
<= 5% of mapped length, i.e. 0MM for 16-19b, 1MM for 20-39b etc,
splicing switched off).

You can clip 3' adapter before feeding the reads to STAR, or you can
use simple built-in clipper --clip3pAdapterSeq TGGAATTCTC
--clip3pAdapterMMp 0.1 (second parameter is the proportion of
mismatches in the matched adapter length).

You would also likely want to filter out reads that STAR "genomically"
trims at the 5'.  This simple awk one-liner will filter out all
alignments that are trimmed by more than 1 base from the 5':

awk '{S=0; split($6,C,/[0-9]*/); n=split($6,L,/[NMSID]/); if
(and($2,0x10)>0 && C[n]=="S") {S=L[n-1]} else if (and($2,0x10)==0 &&
C[2]=="S") {S=L[1]}; if (S<=1) print }' Aligned.out.sam >
Aligned.filtered.sam





Q 

I was wondering if anyone knows how to integrate cufflinks with star.

alexdobin

For non-strand-specific data, you need to use STAR option
--outSAMstrandField intronMotif which will add the XS attribute to all
canonically spliced alignments using their introns' motifs - that's
exactly what Cufflinks needs.

For strand-specific data, you do not need any extra parameters for
STAR runs, but you need to use --library-type option for
Cufflinks. For example, for the "standard" dUTP protocol you need to
use --library-type fr-firststrand in Cufflinks.

Do not use the option --outSAMattributes All, as Cufflinks (as well as
many other tools) have compatibility problems with the latest SAM
format specifications.

For poly-A+ selected RNAs, all other default parameters generally work
well in our ENCODE experience.  As always I would highly recommend
using annotations for mapping.  Also, --outFilterType BySJout could
help clean up rare junctions, and --outFilterIntronMotifs
RemoveNoncanonicalUnannotated will clean up non-canonical junctions.

Total RNA data is much harder to assemble with Cufflinks since it
contains a substantial intronic signal. Please check this post
(https://groups.google.com/d/msg/rna-star/X8mjUc7nm1U/dl7Jy12BNAIJ),
you may need to mask "complex" loci.





Notes on strandedness and cufflinks
-----------------------------------

Looking at the stranded library SRX028191.

STAR with --outSAMstrandField intronMotif = star_out_outSAMstrandField/
  cufflinks with --library-type default     = cufflinks_outSAMstrandField/
    - single exon gene WBGene00023193 at the start of Chrom I has FPKM = 0 because no reads in BAM file for it
    - probably because we are restricting reads to 2 unique hits
    - otherwise broadly in line with Tophat results

  cufflinks with --library-type fr-firststrand = cufflinks_outSAMstrandField_library-type/
    - same single-exon with FPKM = 0 because no reads in BAM file for it
    - probably because we are restricting reads to 2 unique hits
    - many genes with FPKM = 0 <= PROBLEM!!!!!!!!
    - only doing the ones with a splice junction <= PROBLEM!!!!!!!!
    - and the FPKM values are very low <= PROBLEM!!!!!!!!

STAR with default options (no --outSAMstrandField intronMotif) = star_out_default/
  cufflinks with --library-type default = cufflinks_default_default/
    - cufflinks complains long and loud about "BAM record error: found spliced alignment without XS attribute" <= PROBLEM!!!!!!!!

  cufflinks with --library-type fr-firststrand = cufflinks_default_library-type/
        - same single-exon with FPKM = 0 because no reads in BAM file for it
    - probably because we are restricting reads to 2 unique hits
    - many genes with FPKM = 0 <= PROBLEM!!!!!!!!
    - only doing the ones with a splice junction <= PROBLEM!!!!!!!!
    - and the FPKM values are very low <= PROBLEM!!!!!!!!



Notes on sorting using samtools
-------------------------------

From SEQAnswers:
http://seqanswers.com/forums/showthread.php?t=29652

This merged bam file is rather large, 200GB in size 

Originally, I was using samtools 0.1.18, and with the sort command I
set the -m flag to 16000000000. This took ~12 hrs to finish.

After updating to samtools 0.1.19 and using the sort command with -m
2G and -@ 32, the file sorted in about 4hrs...

Useful options for cufflinks described in
-----------------------------------------
http://seqanswers.com/forums/archive/index.php/t-20702.html


--max-bundle-frags 999999999
--no-effective-length-correction
--min-isoform-fraction 0
--pre-mrna-fraction 0.05
--junc-alpha 0.05
--max-bundle-length 5500000
-b <genome.fa>
-G <my_annotation.gtf>


I've set --max-bundle-frags and --max-bundle-length to these values in
order to force it to evaluate all of the bundles it finds in the mouse
data that I usually work with. So you would want to tweak those
depending on your depth of sequencing and the genome you're working
with. Just watch for "warning skipping ...." in the cufflinks output
it generates while running.

Per Cole in this same thread I've disabled effective length correction
(--no-effective-length-correction) because I don't like that this
correction extrapolates the raw data. It also greatly exaggerates the
expression estimation of single exon genes that are in the gene
annotations. Basically disabling this feature makes it so the
counts/expression seem to match up with the raw coverage a little
better.

Bias correction (-b) makes a big difference in the quality of the
quantifications. It takes longer to run but it seems to really help
cufflinks produce logical output.

--min-isoform-fraction and --pre-mrna-fraction are both settings that
  allow cufflinks to discard information based on some arbirary
  thresholding. --min-isoform-fraction filters out low-count isoforms
  and --pre-mrna-fraction filters out intronic alignment data. I'm not
  sure what a difference --pre-mrna-fraction makes but when I was
  picking out these options I was looking for any options that made
  any filtering less strict. My and the researchers I work with would
  rather do the filtering afterwards. It's more useful for us to see
  as much of the raw data as possible.

--junc-alpha is one that I haven't tested with different values but I
  have set it to be less strict. The default is 0.001. I plan to mess
  around with that one a little more to see what sort of impact it
  has. I'm not sure if it even has any impact at all on quantification
  (maybe it only applies to assembly).

I believe --min-isoform-fraction has an impact on what is reported
when you use cufflinks to do assembly. For example if you aren't
prepared to trust cufflinks' estimations of expression quantification
but you want to use it to build the most robust assembly possible then
by setting this value to 0 you ensure that cufflinks doesn't throw out
assembled isoforms it thinks are low expressed. I assume --junc-alpha
will have some kind of impact as well but I haven't tested it out
much.

I used the above settings to quantify some real data against the mm9
known gene annotation (from UCSC) and compared those quantifications
to those I got using eXpress. As I posted before - it seemed like
cufflinks was doing a much better job of making specific decisions
about which isoforms really needed to be expressed to fully explain
the alignments. For example if all of the coverage and junction
information in a locus can be explained by a single isoform why should
these programs report that multiple isoforms are expressed? If they do
then I think that decreases the sensitivity of differential splicing
analysis. Maybe in another sample there's new junction and coverage
information in that locus that DOES justify expression of a second
isoform. While eXpress would have givin you expression of both
isoforms in both cases in my tests cufflinks would more likely report
that second isoform as something that was activated in the second
sample. That translates to me that there was sufficient splicing or
coverage evidence of that new isoform and not just that some
proportion of reads are being assigned to it in both cases because
they share exons.

Useful description of strandedness of mate pairs
------------------------------------------------

http://www.cureffi.org/2012/12/19/forward-and-reverse-reads-in-paired-end-sequencing/

Therefore when you open your FASTQ files and look at a pair of reads,
the sequences you see are, conceptually, pointing towards each other
on opposite strands.  When you align them to the genome, one read
should align to the forward strand, and the other should align to the
reverse strand, at a higher base pair position than the first one so
that they are pointed towards one another.  This is known as an “FR”
read – forward/reverse, in that order.

Useful description of strandedness of mate pairs
------------------------------------------------

http://onetipperday.blogspot.co.uk/2012/07/how-to-tell-which-library-type-to-use.html

Comment on strandedness and cufflinks
-------------------------------------

http://seqanswers.com/forums/showpost.php?p=41795&postcount=9

I have asked one of the developers of Tophat about this library type option:

"If you don't specify --library-type, TopHat just treats your reads as
unstranded. (The default is *unstranded*). 

Actually, library-type is only intended for paired-end, so if you
specify the library-type fr-unstranded option, should be the same as
non-specified. "
