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

const my $LINKML_SCHEMA => 'v1.9.0';

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
	$outfile = "./wormbase.constructs.${ws_version}.test_data.json";
    } else {
	$outfile = "./wormbase.constructs.${ws_version}.${LINKML_SCHEMA}.json";
    }
}

my $db = Ace->connect(-path => $acedbpath, -program => $tace) or die("Connection failure: ". Ace->error);

my $out_fh  = file($outfile)->openw;

process_constructs($db, $out_fh);

$db->close;

sub process_constructs {
    my ($db, $out_fh) = @_;

    my @constructs = $db->fetch(-query => "find Construct WHERE (Public_name OR Summary OR Driven_by_gene OR Gene)");

    my $all_annots = {linkml_version => $LINKML_SCHEMA};
    $all_annots->{construct_ingest_set => []};
    for my $obj (@constructs) {
	next unless $obj->isObject();
	my $dp_xref_dto_json = {
	    referenced_curie => 'WB:' . $obj->name,
	    page_area => 'construct',
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

	my $is_obsolete = ($obj->Status && $obj->Status->name eq 'Dead') ? JSON::true : JSON::false;
	    
	my $json_obj = {
	    mod_entity_id              => "WB:" . $obj->name,
	    name                       => $obj->Public_name ? "${\$obj->Public_name}" : "$obj",
	    internal                   => JSON::false,
	    obsolete                   => $is_obsolete,
	    created_by_curie           => 'WB:curator',
	    updated_by_curie           => 'WB:curator',
	    data_provider_dto          => $data_provider_dto_json
	};

	if ($obj->Gene || $obj->Driven_by_gene) {
	    $json_obj->{construct_components}=[];
	    my $driven_by_gene_components = get_components($obj->Driven_by_gene);
	    my $gene_components = get_components($obj->Gene);
	    push @{$json_obj->{construct_components}}, @$gene_components;
	    push @{$json_obj->{construct_components}}, @$driven_by_gene_components;
	}
	
	if ($curation_test) {
	    $json_obj->{taxon_curie} = "NCBITaxon:6239";
	    $json_obj->{date_updated} = get_random_datetime(10);
	    $json_obj->{date_created} = get_random_datetime(1);  
	} 
	push @{$all_annots->{construct_ingest_set}}, $json_obj;
    }	

    my $json_obj = JSON->new;
    my $string = $json_obj->allow_nonref->canonical->pretty->encode($all_annots);
    $out_fh->print($string);
    
    
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
	$publication_id = 'WB:WBPaper00045183';
    }
    if ($publication_id eq 'WB:WBPaper000042571') {
	$publication_id = 'WB:WBPaper00042571';
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

sub get_components {
    my $obj = shift;

    my @components;
    for my $g($obj) {
	my $component_json = {
	    componentSymbol => $g->Public_name->name;
	};

	my @component_paper_evidence;
	for my $evi ($g->right->col) {
	    next unless $evi->name eq 'Paper_evidence';
	    for my $paper ($evi->col) {
		my $paper_id = get_paper($paper);
		push @component_paper_evidence, $paper_id if defined $paper_id;
	    }
	}
	$component_json->{evidence_curies} = \@component_paper_evidence if @component_paper_evidence > 0;
	
	my @note_dtos;
	my $cmp_note_json = {
	    free_text => $g->Public_name->name . ' test note',
	    note_type_name => 'construct_component_note',
	    internal => JSON::false,
	    obsolete => JSON::false
	};
	push @note_dtos, $cmp_note_json;
	$component_json->{related_notes} = \@note_dtos if $curation_test;
	$component_json->{taxon_curie} = "NCBITaxon:6239" if $curation_test;
	$component_json->{taxon_text} = "Caenorhabditis elegans" if $curation_test;

	push @components, $component_json;
    }
    return \@components;
}

    
