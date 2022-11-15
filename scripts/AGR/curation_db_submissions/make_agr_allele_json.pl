#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;
use Path::Class;

use lib $ENV{CVS_DIR};
use Wormbase;
use Modules::AGR;

my ($debug, $test, $verbose, $store, $wormbase, $curation_test, $limit, $schema);
my ($outfile, $acedbpath, $ws_version, $out_fh, $bgi_json,$disease_file);

GetOptions (
    'debug=s'       => \$debug,
    'test'          => \$test,
    'verbose'       => \$verbose,
    'store:s'       => \$store,
    'database:s'    => \$acedbpath,
    'outfile:s'     => \$outfile,
    'curationtest' => \$curation_test,  
    'wsversion=s'   => \$ws_version,
    'limit:i'       => \$limit,
    'schema=s'      => \$schema
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
$out_fh->print("{\n   \"linkml_version\" : \"" . $schema . "\",\n   \"allele_ingest_set\" : [");

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
	
	my $json_obj = {
	    curie             => "WB:" . $obj->name,
	    symbol            => $obj->Public_name ? $obj->Public_name->name : $obj->name,
	    taxon_curie       => "NCBITaxon:" . $obj->Species->NCBITaxonomyID->name,
	    internal          => JSON::false,
	    obsolete          => $obj->Status && $obj->Status->name eq 'Live' ? JSON::false : JSON::true,
	    created_by_curie  => 'WB:curator',
	    updated_by_curie  => 'WB:curator'
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
	
	#if ($obj->Other_name) {
	#   my @synonyms_to_submit;
	#    my %synonyms = map {$_->name => 1} $obj->Other_name;
	#    for my $syn (keys %synonyms) {
	#	push @synonyms_to_submit, {name => $syn} unless $syn =~ /^[^:]+:[pcg]\./; # remove HGVS identifiers (already generated in Alliance)
	#    }
	#	
	#    $json_obj->{synonyms} = \@synonyms_to_submit if @synonyms_to_submit > 0;
	#}
	if ($obj->Reference) {
	    $json_obj->{reference_curies} = get_papers($obj);
	}

	# Stick some random data in for curation DB test files
	if ($curation_test) {
	    my @inheritance_modes = qw(dominant recessive semi-dominant unknown);
	    $json_obj->{inheritance_mode_name} = $inheritance_modes[rand @inheritance_modes];
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
	
	my $json_obj = {
	    curie         => "WB:" . $obj->name, 
	    symbol        => $obj->Public_name ? $obj->Public_name->name : $obj->name,
	    taxon_curie   => "NCBITaxon:" . $obj->Species->NCBITaxonomyID->name,
	    internal      => JSON::false,
	    obsolete      => JSON::false
	};	
	#if ($obj->Synonym) {
	#    my %synonyms = map {$_->name => 1}$obj->Synonym;
	#    if (scalar keys %synonyms > 0) {
	#	my @synonyms_to_submit;
	#	for my $syn (keys %synonyms) {
	#	    push @synonyms_to_submit, {name => $syn};
	#	}
	#	$json_obj->{synonyms} = \@synonyms_to_submit;
	#    }
	#}
	if ($obj->Reference) {
	    $json_obj->{reference_curies} = get_papers($obj);
	}

	
	# Stick some random data in for curation DB test files
	if ($curation_test) {
	    my @inheritance_modes = qw(dominant recessive semi-dominant unknown);
	    $json_obj->{inheritance_mode_name} = $inheritance_modes[rand @inheritance_modes];
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
	my $pmid;
	for my $db ($ref->Database) {
	    if ($db->name eq 'MEDLINE') {
		$pmid = $db->right->right->name;
		last;
	    }
	}
	my $publication_id = $pmid ? "PMID:$pmid" : 'WB:' . $ref->name;
	push @references, $publication_id;
    }

    return \@references;
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
		    push @substitution_paper_evidence, 'WB:' . $paper->name;
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
		push @insertion_paper_evidence, 'WB:' . $paper->name;
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
		push @deletion_paper_evidence, 'WB:' . $paper->name;
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
		push @duplication_paper_evidence, 'WB:' . $paper->name;
	    }
	}
	$duplication->{evidence_curies} = \@duplication_paper_evidence if @duplication_paper_evidence > 0;
	push @mutation_types, $duplication
    };

    return \@mutation_types;
}
