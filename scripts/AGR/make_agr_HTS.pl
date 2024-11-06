#!/bin/env perl
#
# script to create the JSON for 
#   https://github.com/alliance-genome/agr_schemas/tree/release-1.0.1.1/ingest/htp
# based on: 
#   https://docs.google.com/spreadsheets/d/1H7bzMJVj6KevZHDGp61tCCdRo8ubHHzEKWaYu_7gn-w

use strict;
use Storable;
use Getopt::Long;
use Ace;
use JSON;
use Storable qw(dclone);

use lib $ENV{CVS_DIR};
use Wormbase;
use Modules::AGR;

# MMO:0000659 - RNA-seq assay
# MMO:0000649 - microarray
# MMO:0000648 - transcript array
# MMO:0000664 - proteomic profiling
# MMO:0000666 - high throughput proteomic profiling
# MMO:0000000 - measurement method



# OBI:0000895 - total RNA extract
# OBI:0000869 - polyA RNA extract
# OBI:0000880 - RNA extract
# OBI:0000862 - nuclear RNA extract
# OBI:0000423 - extract
# OBI:0000876 - cytoplasmic RNA extract
# OBI:0002573 - polyA depleted RNA extract

my ($debug, $test, $verbose, $store, $acedbpath, $outfile,$ws_version,$outFdata,$wb_to_uberon_file);
GetOptions (
  'debug=s'     => \$debug,
  'test'        => \$test,
  'verbose'     => \$verbose,
  'store:s'     => \$store,
  'database:s'  => \$acedbpath,
  'samplefile:s'=> \$outfile,
  'datasetfile:s'=> \$outFdata,
  'wsversion=s' => \$ws_version,
  'wb2uberon=s' => \$wb_to_uberon_file,
)||die("unknown command line option: $@\n");

my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
      );
}

my $tace = $wormbase->tace;
my $date = AGR::get_rfc_date();

$acedbpath  ||= $wormbase->autoace;
$ws_version ||= $wormbase->get_wormbase_version_name;

my $db = Ace->connect(-path => $acedbpath, -program => $tace) or die('Connection failure: '. Ace->error);

my @htps; # high troughput samples
my @htpd; # high thruput datasets

my $assembly = fetch_assembly($db,$wormbase);


my %uberon=%{read_uberron($wb_to_uberon_file)};

my $it = $db->fetch_many(-query => 'find Analysis;Species_in_analysis="Caenorhabditis elegans"') || die (Ace->error);
while (my $analysis = $it->next) {    
    if ($analysis->name =~ /^RNASeq_Study/ or $analysis->Microarray_experiment) {
	my $datasetId = 'WB:' . $analysis->name;
       
	my @p = $analysis->Reference;
	my @papers = map {get_paper($_)} @p;
	
	$db->timestamps(1);
	my $timestamp = $analysis->timestamp;
	$timestamp=~s/_/T/;
	$timestamp=~s/_\w+/\+01:00/;
	$db->timestamps(0);

	my %dataset_json;
	$dataset_json{datasetId}= {primaryId => $datasetId,
				   preferredCrossReference => ,{id =>$datasetId,pages => ['htp/dataset']}
	}; # required
	$dataset_json{publications} = \@papers if @papers;
	$dataset_json{dateAssigned} = $timestamp; #required
	$dataset_json{title} = $analysis->Title ? $analysis->Title->name : "$analysis"; #required
	$dataset_json{summary} = $analysis->Description->name if $analysis->Description;
	if ($analysis->Category) {
	    my @tags = map {$_->name} $analysis->Category;
	    $dataset_json{categoryTags} = \@tags;
	}
	else {
	    $dataset_json{categoryTags} = ['unclassified'];
	}

	my %num_channels;
	if ($analysis->name =~ /^RNASeq_Study/) {
	    foreach my $subproject ($analysis->Subproject){
		# RNASeq sample
		next unless $subproject->Species_in_analysis eq 'Caenorhabditis elegans';
		my %sample_json;
		my $sample=$subproject->Sample;
		
		$sample_json{sampleId} = {
		    primaryId => "WB:$subproject",
		}; # required
		$sample_json{sampleTitle} = $subproject->Title->name;
		$sample_json{sampleType} = 'OBI:0000895'; # total RNA extract
		$sample_json{assayType} = 'MMO:0000659'; # RNA-seq assay
		$sample_json{taxonId} = 'NCBITaxon:' . $wormbase->ncbi_tax_id;
		$sample_json{datasetIds} = [$datasetId];
		$sample_json{sex}=lc($sample->Sex->name) if $sample->Sex;

		sampleAge(\%sample_json,$sample) if $sample->Life_stage;
		sampleLocation(\%sample_json,$sample) if $sample->Tissue;

		$sample_json{assayType} = 'MMO:0000659'; # RNA-seq assay
		$sample_json{assemblyVersions} = [$assembly];
		push @htps,\%sample_json;
	    }
	}
	else {
	    # Microarray sample
	    for my $array ($analysis->Microarray_experiment) {
		if ($array->Microarray_sample) {
		    # 1-channel microarray samples
		    $num_channels{1} = 1;
		    
		    my %sample_json;
		    my $sample=$array->Microarray_sample;
		    $sample_json{sampleId} = {
			primaryId => "WB:$sample",
			secondaryId => ["WB:${\$array->name}"],
		    };
		    $sample_json{sampleTitle} = "WB:$sample";
		    $sample_json{sampleType} = 'OBI:0000880'; # RNA extract
		    
		    
		    $sample_json{taxonId} = 'NCBITaxon:' . $wormbase->ncbi_tax_id;
		    $sample_json{datasetIds} = ['WB:' . $analysis->name];
		    $sample_json{sex}=lc($sample->Sex->name) if $sample->Sex;
		    
		    sampleAge(\%sample_json,$sample) if ($sample->Life_stage);
		    sampleLocation(\%sample_json,$sample) if ($sample->Tissue);
		    
		    $sample_json{assayType}='MMO:0000649'; # micro array
		    $sample_json{assemblyVersions}=[$assembly];
		    $sample_json{microarraySampleDetails} = {channelId => "WB:$sample", channelNum => 1};
		    push @htps,\%sample_json;
		}
		elsif ($array->Sample_A and $array->Sample_B) {
		    # 2-channel microarray samples
		    $num_channels{2} = 1;

		    my $n=1;
		    foreach my $s($array->Sample_A , $array->Sample_B){
			my %sample_json;
			my $suffix = ($n == 1) ? 'A' : 'B'; # for the datasets

			$sample_json{sampleTitle} = "WB:$s:$suffix";
			$sample_json{sampleType} = 'OBI:0000880'; # RNA extract
			$sample_json{sampleId} = {
			    primaryId => "WB:$s",
			    secondaryId => ["WB:${\$array->name}"],
			};

			$sample_json{taxonId} = 'NCBITaxon:' . $wormbase->ncbi_tax_id;
			$sample_json{datasetIds} = ['WB:' . $analysis->name];
			$sample_json{sex}=$s->Sex->name if $s->Sex;

			sampleAge(\%sample_json,$s) if $s->Life_stage;
			sampleLocation(\%sample_json,$s) if $s->Tissue;

			$sample_json{assayType}='MMO:0000649'; # micro array
			$sample_json{assemblyVersions}=[$assembly];
			$sample_json{microarraySampleDetails} = {channelId => "WB:$s:$suffix", channelNum => $n};
			push @htps,\%sample_json;
			
			$n++;
		    }
		}
		else {
		    # Unknown microarray details
		    $num_channels{'?'} = 1;
		}
	    }
	}

	my @channels = keys %num_channels;
	if (scalar @channels == 1 and $channels[0] ne '?') {
	    $dataset_json{numChannels} = int($channels[0]);
	}
	push @htpd,\%dataset_json;
    }
}



################################################

print_json(\@htps,$outfile);
print_json(\@htpd,$outFdata);

$db->close;


####################################
# helper functions

# will change the hashref, so no return

sub sampleAge{
	my ($jsonRef,$sample)=@_;

	my $ls=$sample->Life_stage;

	return unless $ls && $ls->Public_name;

	$jsonRef->{sampleAge}={stage => { stageTermId => $ls->name,
		                          stageName => $ls->Public_name->name ,
					  stageUberonSlimTerm => {uberonTerm => ((keys %{$uberon{$ls->name}})[0]||'post embryonic, pre-adult')},
				        }
	};
}

sub sampleLocation{
	my ($jsonRef,$sample)=@_;
	$jsonRef->{sampleLocations}=[];
	map {push @{$jsonRef->{sampleLocations}},
				{anatomicalStructureTermId => "$_",
				 anatomicalStructureUberonSlimTermIds => [{uberonTerm =>((keys %{$uberon{$_->name}})[0] || 'Other')}],
				 whereExpressedStatement => ($_->Term->name || 'C.elegans Cell and Anatomy'),
			        }
  	} $sample->Tissue;
}

sub read_uberron{
	my ($file) = @_;
        my %wb2u;
	my $fh = IO::File->new($file);
	while (<$fh>){
          if (/^(\S+)\s+(\S+)/){
    		my $uberon = $1;
		my @list = split(/,/, $2);
    		foreach my $item (@list) {
      		  $wb2u{$item}->{"$uberon"} = 1;
		}
	  }
	}
	$wb2u{'WBbt:0000100'}->{Other} = 1;
	return \%wb2u;
}

sub print_json{
	my ($data,$file) = @_;
	print STDERR "printing $file\n";
	my $completeData = {
	  metaData => AGR::get_file_metadata_json( (defined $ws_version) ? $ws_version : $wormbase->get_wormbase_version_name(), $date ),
	  data     => $data,
	};
        my $fh;
	if ($file) {
	  open $fh, ">$file" or die "Could not open $outfile for writing\n";
	} else {
	  $fh = \*STDOUT;
	}

	my $sample_json = JSON->new;
	print $fh $sample_json->allow_nonref->canonical->pretty->encode($completeData);
}

###########
sub get_paper {
  my ($wb_paper) = @_;

  my $json_paper = {};

  my $pmid;
  foreach my $db ($wb_paper->Database) {
    if ($db->name eq 'MEDLINE') {
      $pmid = $db->right->right->name;
      last;
    }
  }
  $json_paper->{publicationId} = $pmid ? "PMID:$pmid" : "WB:$wb_paper";
  $json_paper->{crossReference} = {id =>"WB:$wb_paper",pages => ['reference']};

  return $json_paper;
}
###########

sub fetch_assembly{
	my ($db,$wormbase) = @_;
	# get the assembly name for the canonical bioproject
	my $assembly_name;
  
	my $species_obj = $db->fetch(-class => 'Species', -name => $wormbase->full_name);
	my @seq_col = $species_obj->at('Assembly');
        
	foreach my $seq_col_name (@seq_col) {
		my $bioproj;
		my $seq_col = $seq_col_name->fetch;
		my $this_assembly_name = $seq_col->Name;

		my @db = $seq_col->at('Origin.DB_info.Database');
		foreach my $db (@db) {
       			$bioproj = $db->right->right->name if ($db->name eq 'NCBI_BioProject');
		}
    
		if (defined $bioproj and $wormbase->ncbi_bioproject eq $bioproj) {
			$assembly_name = $this_assembly_name->name;
			last;
		}
        }
    
	if (not defined $assembly_name) {
	    die "Could not find name of current assembly for " . $wormbase->species() . "\n";
	}
        return $assembly_name;
}
