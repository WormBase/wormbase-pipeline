#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;
use Path::Class;
use Const::Fast;

use lib $ENV{CVS_DIR};
use Wormbase;
use Modules::AGR;

my ($debug, $test, $verbose, $store, $wormbase, $curation_test, $limit, $schema);
my ($outfile, $acedbpath, $ws_version, $out_fh, $bgi_json,$disease_file);

const my $LINKML_SCHEMA => 'v1.7.0';

GetOptions (
    'debug=s'       => \$debug,
    'test'          => \$test,
    'verbose'       => \$verbose,
    'store:s'       => \$store,
    'database:s'    => \$acedbpath,
    'outfile:s'     => \$outfile,
    'curationtest' => \$curation_test,  
    'wsversion=s'   => \$ws_version,
    'limit:i'       => \$limit
)||die();

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
      );
}

my $tace = $wormbase->tace;
my $date = AGR::get_rfc_date();
my $alt_date = join("/", $date =~ /^(\d{4})(\d{2})(\d{2})/);
my $full_name = $wormbase->full_name;

$acedbpath  ||= $wormbase->autoace;
$ws_version ||= $wormbase->get_wormbase_version_name;
if (!$outfile) {
    if ($curation_test) {
	$outfile = "./wormbase.alleles.${ws_version}.test_data.json";
    } else {
	$outfile = "./wormbase.alleles.${ws_version}.json";
    }
}

my $db = Ace->connect(-path => $acedbpath, -program => $tace) or die("Connection failure: ". Ace->error);

my $out_fh  = file($outfile)->openw;
$out_fh->print("{\n   \"linkml_version\" : \"" . $LINKML_SCHEMA . "\",\n   \"allele_ingest_set\" : [");

my $alleles = process_variations($db, $out_fh);

$alleles = process_transgenes($db, $out_fh);
$out_fh->print("\n   ]\n}");

$db->close;

sub process_variations {
    my ($db, $out_fh) = @_;

    my @alleles = $db->fetch(-query => "Find Variation WHERE species = \"Caenorhabditis elegans\"");
    
    my $var_count = 0;
    for my $obj (@alleles) {
	$var_count++;
	next unless $obj->isObject();
	unless ($obj->Species) {
	    print "No species for $obj - skipping\n";
	    next;
	}

	my $dp_xref_dto_json = {
	    referenced_curie => 'WB:' . $obj->name,
	    page_area => 'allele',
	    display_name => $obj->name,
	    prefix => 'WB',
	    internal => JSON::false,
	    obsolete => JSON::false
	};
	
	my $data_provider_dto_json = {
	    source_organization_abbreviation => 'WB',
	    cross_reference_dto => $dp_xref_dto_json,
	    internal => JSON::false,
	    obsolete => JSON::false
	}; 
	
	my $json_obj = {
	    curie             => "WB:" . $obj->name,
	    allele_symbol_dto => get_symbol_dto($obj),
	    taxon_curie       => "NCBITaxon:" . $obj->Species->NCBITaxonomyID->name,
	    internal          => JSON::false,
	    obsolete          => $obj->Status && $obj->Status->name eq 'Live' ? JSON::false : JSON::true,
	    created_by_curie  => 'WB:curator',
	    updated_by_curie  => 'WB:curator',
	    data_provider_dto => $data_provider_dto_json
	};
	if ($obj->Method) {
	    my $collection = $obj->Method->name;
	    unless ($collection eq 'Allele'
		    || $collection eq 'CGH_allele'
		    || $collection eq 'Deletion_allele'
		    || $collection eq 'Deletion_and_insertion_allele'
		    || $collection eq 'Deletion_polymorphism'
		    || $collection eq 'Engineered_allele'
		    || $collection eq 'Insertion_allele'
		    || $collection eq 'Insertion_polymorphism'
		    || $collection eq 'Mos_insertion'
		    || $collection eq 'SNP'
		    || $collection eq 'Substitution_allele'
		    || $collection eq 'Transposon_insertion'
		){
		$collection =~ s/ /_/g;
		$collection .= '_project' if $collection eq 'Million_mutation';
		$json_obj->{in_collection_name} = $collection;
	    }
	}

	my $mutation_types = get_mutation_types($obj);
	$json_obj->{allele_mutation_type_dtos} = $mutation_types if scalar @$mutation_types > 0;

	my @synonym_dtos;
	if ($obj->Other_name) {
	    for my $other_name ($obj->Other_name) {
		next if $other_name->name =~ /^[^:]+:[pcg]\./;
		my $synonym_dto = {
		    display_text   => $other_name->name,
		    format_text    => $other_name->name,
		    name_type_name => 'unspecified'
		};
		my $synonym_evidence = get_evidence_curies($other_name);
		$synonym_dto->{evidence_curies} = $synonym_evidence if @$synonym_evidence;
		push @synonym_dtos, $synonym_dto;
	    }
	}
	$json_obj->{allele_synonym_dtos} = \@synonym_dtos if @synonym_dtos;
	
	if ($obj->Reference) {
	    $json_obj->{reference_curies} = get_papers($obj);
	}

	# Stick some random data in for curation DB test files
	if ($curation_test) {
	    my @inheritance_modes = qw(dominant recessive semi-dominant unknown);
	    $json_obj->{is_extinct} = $var_count % 2 == 0 ? JSON::true : JSON::false;
	    $json_obj->{date_updated} = get_random_datetime(10);
	    $json_obj->{date_created} = get_random_datetime(1);
	}

	my $json = JSON->new;
	my $string = $json->allow_nonref->canonical->pretty->encode($json_obj);
	$out_fh->print(',') unless $var_count == 1;
	$out_fh->print("\n" . '      ' . $string);
	
	last if $limit && $limit <= $var_count * 2;
    }

    return;
}

sub process_transgenes {
    my ($db, $out_fh) = @_;


    my @transgenes = $db->fetch(-query => 'find Transgene WHERE Species = "Caenorhabditis elegans"');

    my $var_count = 0;
    for my $obj(@transgenes) {
	$var_count++;
	next unless $obj->isObject();
	unless ($obj->Species) {
	    print "No species for $obj - skipping\n";
	    next;
	}
	
	my $dp_xref_dto_json = {
	    referenced_curie => 'WB:' . $obj->name,
	    page_area => 'transgene',
	    display_name => $obj->name,
	    prefix => 'WB',
	    internal => JSON::false,
	    obsolete => JSON::false
	};
	
	my $data_provider_dto_json = {
	    source_organization_abbreviation => 'WB',
	    cross_reference_dto => $dp_xref_dto_json,
	    internal => JSON::false,
	    obsolete => JSON::false
	}; 
	
	my $json_obj = {
	    curie         => "WB:" . $obj->name, 
	    allele_symbol_dto => get_symbol_dto($obj),
	    taxon_curie   => "NCBITaxon:" . $obj->Species->NCBITaxonomyID->name,
	    internal      => JSON::false,
	    obsolete      => JSON::false,
	    data_provider_dto => $data_provider_dto_json
	};	

	my @synonym_dtos;
	if ($obj->Synonym) {
	    for my $synonym ($obj->Synonym) {
		my $synonym_dto = {
		    display_text   => $synonym->name,
		    format_text    => $synonym->name,
		    name_type_name => 'unspecified',
		    internal      => JSON::false
		};
		push @synonym_dtos, $synonym_dto;
	    }
	}
	$json_obj->{allele_synonym_dtos} = \@synonym_dtos if @synonym_dtos;
	


	if ($obj->Reference) {
	    $json_obj->{reference_curies} = get_papers($obj);
	}

	
	# Stick some random data in for curation DB test files
	if ($curation_test) {
	    $json_obj->{is_extinct} = $var_count % 2 == 0 ? JSON::true : JSON::false;
	    $json_obj->{date_updated} = get_random_datetime(10);
	    $json_obj->{date_created} = get_random_datetime(1);
	}
	
	my $json = JSON->new;
	my $string = $json->allow_nonref->canonical->pretty->encode($json_obj);
	$out_fh->print(",\n" . '      ' . $string);
	
	last if $limit && $limit <= $var_count * 2;
    }

    return;
}


sub get_papers {
    my $obj = shift;

    my @references;
    for my $ref ($obj->Reference) {
	next unless defined $ref;
	my $publication_id = get_paper($ref);
	next unless defined $publication_id;
	push @references, $publication_id;
    }

    return \@references;
}

sub get_paper {
    my $ref = shift;

    if (!$ref->Status) {
	return;
    }

    my $level = 1;
    while ($ref->Merged_into && $level < 6) {
	$level++;
	$ref = $ref->Merged_into;
    }
    return if $ref->Status->name eq 'Invalid';

    my $pmid;
    for my $db ($ref->Database) {
	if ($db->name eq 'MEDLINE') {
	    $pmid = $db->right->right->name;
	    last;
	}
    }
    my $publication_id = $pmid ? "PMID:$pmid" : 'WB:' . $ref->name;
    if ($publication_id eq 'WB:WBPaper000045183') {
	$publication_id eq 'WB:WBPaper00045183';
    }

    return $publication_id;
}

sub get_random_datetime {
    my $add = shift;
    
    my $year = int(rand(10)) + $add;
    my $month = int(rand(11)) + 1;
    my $day = int(rand(27)) + 1;

    my $year_string = $year > 9 ? $year : '0' . $year;
    my $month_string = $month > 9 ? $month : '0' . $month;
    my $day_string = $day > 9 ? $day : '0' . $day;
    
    return '20' . $year_string . '-' . $month_string . '-' . $day_string . 'T00:00:00+00:00';
}

sub get_mutation_types {
    my $obj = shift;

    my @mutation_types;
    if ($obj->Substitution) {
	my $substitution = {
	   mutation_type_curies => ["SO:1000002"],
	   internal => JSON::false
	};
	if ($obj->Substitution->right) {
	    my @substitution_paper_evidence;
	    for my $evi ($obj->Substitution->right->col) {
		next unless $evi->name eq 'Paper_evidence';
		for my $paper ($evi->col) {
		    my $paper_id = get_paper($paper);
		    push @substitution_paper_evidence, $paper_id if defined $paper_id;
		}
	    }
	    $substitution->{evidence_curies} = \@substitution_paper_evidence if @substitution_paper_evidence > 0;
	}
	push @mutation_types, $substitution;
    }

    if ($obj->Insertion) {
	my $insertion = {
	    mutation_type_curies => ["SO:0000667"],
	    internal => JSON::false
	};
	my @insertion_paper_evidence;
	for my $evi ($obj->Insertion->col) {
	    next unless $evi->name eq 'Paper_evidence';
	    for my $paper ($evi->col) {
		my $paper_id = get_paper($paper);
		push @insertion_paper_evidence, $paper_id if defined $paper_id;
	    }
	}
	$insertion->{evidence_curies} = \@insertion_paper_evidence if @insertion_paper_evidence > 0;
	push @mutation_types, $insertion;
    }

    if ($obj->Deletion) {
	my $deletion = {
	    mutation_type_curies => ["SO:0000159"],
	    internal => JSON::false
	};
	my @deletion_paper_evidence;
	for my $evi ($obj->Deletion->col) {
	    next unless $evi->name eq 'Paper_evidence';
	    for my $paper ($evi->col) {
		my $paper_id = get_paper($paper);
		push @deletion_paper_evidence, $paper_id if defined $paper_id;
	    }
	}
	$deletion->{evidence_curies} = \@deletion_paper_evidence if @deletion_paper_evidence > 0;
	push @mutation_types, $deletion;
    }

    if ($obj->Tandem_duplication) {
	my $duplication = {
	    mutation_type_curies => ["SO:1000173"],
	    internal => JSON::false
	};
	my @duplication_paper_evidence;
	for my $evi ($obj->Tandem_duplication->col) {
	    next unless $evi->name eq 'Paper_evidence';
	    for my $paper ($evi->col) {
		my $paper_id = get_paper($paper);
		push @duplication_paper_evidence, $paper_id if defined $paper_id;
	    }
	}
	$duplication->{evidence_curies} = \@duplication_paper_evidence if @duplication_paper_evidence > 0;
	push @mutation_types, $duplication
    };

    return \@mutation_types;
}

sub get_symbol_dto {
    my $obj = shift;

    my $symbol_dto;
    if ($obj->Public_name) {
	$symbol_dto = {
	    display_text       => $obj->Public_name->name,
	    format_text        => $obj->Public_name->name,
	    name_type_name     => 'nomenclature_symbol',
	    synonym_scope_name => 'exact',
	    internal      => JSON::false
	};
	my $symbol_evidence = get_evidence_curies($obj->Public_name);
	$symbol_dto->{evidence_curies} = $symbol_evidence if @$symbol_evidence;
    } else {
	$symbol_dto = {
	    display_text       => $obj->name,
	    format_text        => $obj->name,
	    name_type_name     => 'systematic_name',
	    synonym_scope_name => 'exact',
	    internal      => JSON::false
	};
    }

    return $symbol_dto;
}

sub get_evidence_curies {
    my $name_obj = shift;

    my @paper_evidence;
    for my $evi ($name_obj->col) {
	next unless $evi->name eq 'Paper_evidence';
	for my $paper ($evi->col) {
	    my $paper_id = get_paper($paper);
	    push @paper_evidence, $paper_id if defined $paper_id;
	}
    }

    return \@paper_evidence;
}

    
