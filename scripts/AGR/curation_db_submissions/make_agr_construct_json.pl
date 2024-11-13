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
my ($outdir, $acedbpath, $ws_version, $out_fh, $bgi_json,$disease_file);

const my $LINKML_SCHEMA => 'v2.8.1';

GetOptions (
    'debug=s'       => \$debug,
    'test'          => \$test,
    'verbose'       => \$verbose,
    'store:s'       => \$store,
    'database:s'    => \$acedbpath,
    'outdir:s'     => \$outdir,
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
my ($outfile, $assoc_outfile);
if (!$outdir) {
    $outdir = '.';
}
if ($curation_test) {
    $outfile = $outdir . "/wormbase.constructs.${ws_version}.${LINKML_SCHEMA}.test_data.json";
    $assoc_outfile = $outdir . "/wormbase.construct_associations.${ws_version}.${LINKML_SCHEMA}.test_data.json";
} else {
    $outfile = $outdir . "/wormbase.constructs.${ws_version}.${LINKML_SCHEMA}.json";
    $assoc_outfile = $outdir . "/wormbase.construct_associations.${ws_version}.${LINKML_SCHEMA}.json";
}

my $db = Ace->connect(-path => $acedbpath, -program => $tace) or die("Connection failure: ". Ace->error);

my $out_fh = file($outfile)->openw;
my $assoc_out_fh = file ($assoc_outfile)->openw;

process_constructs($db, $out_fh);

$db->close;

sub process_constructs {
    my ($db, $out_fh) = @_;
    
    my @constructs = $db->fetch(-query => "find Construct");
    
    my $all_annots = {linkml_version => $LINKML_SCHEMA, alliance_member_release_version => $ws_version};
    $all_annots->{construct_ingest_set => []};
    my $cnst_count = 0;
    my $associations = {linkml_version => $LINKML_SCHEMA, alliance_member_release_version => $ws_version};
    $associations->{construct_genomic_entity_association_ingest_set => []};
    
    for my $obj (@constructs) {
	$cnst_count++;
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
	
	my $json_obj = {
	    mod_entity_id              => "WB:" . $obj->name,
	    construct_symbol_dto       => get_symbol_dto($obj),
	    internal                   => JSON::false,
	    obsolete                   => JSON::false,
	    created_by_curie           => 'WB:curator',
	    updated_by_curie           => 'WB:curator',
	    data_provider_dto          => $data_provider_dto_json
	};
	
	$json_obj->{construct_full_name_dto} = get_full_name_dto($obj) if $obj->Summary;

	my $synonyms = get_synonym_dtos($obj);
	if (@$synonyms) {
	    $json_obj->{construct_synonym_dtos} = $synonyms;
	}
	
	if ($obj->Reference) {
	    $json_obj->{reference_curies} = get_papers($obj);
	}
	
	if ($obj->Gene || $obj->Driven_by_gene || $obj->Sequence_feature) {
	    $json_obj->{construct_component_dtos}=[];
	    my @notes;
#	    if ($obj->Fusion_reporter) {
#		push @notes, 'Fusion reporter: ' . $obj->Fusion_reporter->name;
#	    }
#	    if ($obj->Other_reporter) {
#		push @notes, 'Reporter: ' . $obj->Other_reporter->name;
#	    }
#	    if ($obj->Purification_tag) {
#		push @notes, 'Purification tag: ' . $obj->Purification_tag->name;
#	    }
#	    if ($obj->Recombination_site) {
#		push @notes, 'Recombination site: ' . $obj->Recombination_site->name;
#	    }
#	    if ($obj->Type_of_construct) {
#		push @notes, 'Construct type: ' . $obj->Type_of_construct->name;
#	    }
#	    if ($obj->Selection_marker) {
#		push @notes, 'Selection marker: ' . $obj->Selection_marker->name;
#	    }
#	    if ($obj->Construction_summary) {
#		push @notes, 'Construction sumary: ' . $obj->Construction_summary->name;
#	    }
#	    if ($obj->DNA_text) {
#		push @notes, 'DNA: ' . $obj->DNA_text->name;
#	    }
	    
	    if (scalar @notes == 0 && $curation_test) {
		push @notes, $json_obj->{name} . ' test note';
	    }
	    
	    if ($obj->Driven_by_gene) {
		for my $gene ($obj->Driven_by_gene) {
		    push @{$associations->{construct_genomic_entity_association_ingest_set}}, get_gene_component('is_regulated_by', $gene, \@notes, $obj) if $gene->Species and $gene->Species->name eq 'Caenorhabditis elegans';
		}
	    }
	    if ($obj->Gene) {
		for my $gene ($obj->Gene) {
		    push @{$associations->{construct_genomic_entity_association_ingest_set}}, get_gene_component('expresses', $gene, \@notes, $obj) if $gene->Species and $gene->Species->name eq 'Caenorhabditis elegans';
		}
	    }
	    if ($obj->UTR_3) {
		for my $gene ($obj->UTR_3) {
		    push @{$associations->{construct_genomic_entity_association_ingest_set}}, get_gene_component("expresses_3'_UTR_of", $gene, \@notes, $obj) if $gene->Species and $gene->Species->name eq 'Caenorhabditis elegans';
		}
	    }
	    # Sequence features will go as own slot annotation in future
	    #if ($obj->Sequence_feature) {
	    #     for my $feature ($obj->Sequence_feature) {
	    #         push @{$json_obj->{construct_component_dtos}}, get_feature_component($feature, \@notes, $obj);
	    #     }
	    #}
	}
	
	if ($curation_test) {
	    $json_obj->{date_updated} = get_random_datetime(10);
	    $json_obj->{date_created} = get_random_datetime(1);
	    $json_obj->{secondary_identifiers} = [$obj->name . '_secondary_test_two', $obj->name . '_secondary_test_one'];
	} 
	push @{$all_annots->{construct_ingest_set}}, $json_obj;
    }	
    
    my $json = JSON->new;
    $json->allow_blessed(1);
    my $string = $json->allow_nonref->canonical->pretty->encode($all_annots);
    $out_fh->print($string);
    
    my $assoc_json = JSON->new;
    $assoc_json->allow_blessed(1);
    my $assoc_string = $json->allow_nonref->canonical->pretty->encode($associations);
    $assoc_out_fh->print($assoc_string);
    
    return;
}

sub get_symbol_dto {
    my $obj = shift;

    my $symbol_dto;
    if ($obj->Public_name) {
	$symbol_dto = {
	    display_text       => $obj->Public_name->name,
	    format_text        => $obj->Public_name->name,
	    name_type_name     => 'nomenclature_symbol',
	    internal      => JSON::false
	};
    } else {
	$symbol_dto = {
	    display_text       => $obj->name,
	    format_text        => $obj->name,
	    name_type_name     => 'systematic_name',
	    internal      => JSON::false
	};
    }
    if ($curation_test) {
	if ($obj->Reference) {
	    $symbol_dto->{evidence_curies} = get_papers($obj);
	}
	$symbol_dto->{synonym_url} = "http://symboltest.org/" . $obj->name;
	$symbol_dto->{synonym_scope_name} = "exact";
    }
    

    return $symbol_dto;
}

sub get_synonym_dtos {
    my $obj = shift;
    my @synonym_dtos;
    
    if ($obj->Other_name) {
	for my $on ($obj->Other_name) {
	    my $synonym_dto = {
		display_text => $on->name,
		format_text => $on->name,
		name_type_name => 'unspecified'
	    };
	    push @synonym_dtos, $synonym_dto;
	}
    }
		    
    return \@synonym_dtos;
}

sub get_full_name_dto {
    my $obj = shift;

    my $fn_dto;
    if ($obj->Summary) {
	$fn_dto = {
	    display_text => $obj->Summary->name,
	    format_text => $obj->Summary->name,
	    name_type_name => 'full_name',
	    internal => JSON::false
	};
	my @evidence;
	for my $evi ($obj->Summary->col) {
	    next unless $evi->name eq 'Paper_evidence';
	    for my $paper ($evi->col) {
		my $paper_id = get_paper($paper);
		push @evidence, $paper_id if defined $paper_id;
	    }
	}
	$fn_dto->{evidence_curies} if @evidence > 0;
    }

    return $fn_dto;
}


sub get_papers {
    my ($obj) = @_;
    
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

sub get_feature_component {
    my ($feature, $notes, $obj) = @_;

    my $component_json = {
	relation_name => 'expresses',
	component_symbol => $feature->Public_name ? $feature->Public_name->name : $feature->name,
	internal => JSON::false,
	obsolete => JSON::false
    };

    my @note_dtos;
    for my $note (@$notes) {
	my $cmp_note_dto = {
	    free_text => $note,
	    note_type_name => 'comment',
	    internal => JSON::false,
	    obsolete => JSON::false
	};
	if ($curation_test && $obj->Reference) {
	    $cmp_note_dto->{evidence_curies} = get_papers($obj);
	}
	push @note_dtos, $cmp_note_dto;
    }
    $component_json->{note_dtos} = \@note_dtos if scalar @note_dtos > 0;
   

    return $component_json;
}

sub get_gene_component {
    my ($relation, $gene, $notes, $obj) = @_;
    
    my $component_json = {
	construct_identifier => 'WB:' . $obj->name,
	genomic_entity_relation_name => $relation,
	genomic_entity_identifier  => 'WB:' . $gene->name,
	internal                   => JSON::false,
	obsolete                   => JSON::false,
	created_by_curie           => 'WB:curator',
	updated_by_curie           => 'WB:curator',
    };
	
    my @component_paper_evidence;
    for my $evi ($gene->col) {
	next unless $evi->name eq 'Paper_evidence';
	for my $paper ($evi->col) {
	    my $paper_id = get_paper($paper);
	    push @component_paper_evidence, $paper_id if defined $paper_id;
	}
    }
    if ($curation_test && scalar @component_paper_evidence == 0) {
	push @component_paper_evidence, "WB:WBPaper00042571";
    }
    $component_json->{evidence_curies} = \@component_paper_evidence if scalar @component_paper_evidence > 0;
	
    my @note_dtos;
    for my $note (@$notes) {
	my $cmp_note_dto = {
	    free_text => $note,
	    note_type_name => 'comment',
	    internal => JSON::false,
	    obsolete => JSON::false
	};
	if ($curation_test && $obj->Reference) {
	    $cmp_note_dto->{evidence_curies} = get_papers($obj);
	}
	push @note_dtos, $cmp_note_dto;
    }
    $component_json->{note_dtos} = \@note_dtos if scalar @note_dtos > 0;
   	
    return $component_json;
}


