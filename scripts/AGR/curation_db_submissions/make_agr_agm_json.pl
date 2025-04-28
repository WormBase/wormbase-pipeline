#!/usr/bin/env perl
use strict;
use Getopt::Long;
use JSON;
use Ace;
use Wormbase;
use Const::Fast;

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($outfile, $acedbpath, $ws_version, $out_fh);

const my $LINKML_SCHEMA => 'v2.12.0';

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "test"        => \$test,
	    "verbose"     => \$verbose,
	    "store:s"     => \$store,
	    "database:s"  => \$acedbpath,
	    "outfile:s"   => \$outfile,
            "wsversion=s" => \$ws_version
	    );

if ( $store ) {
    $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
    $wormbase = Wormbase->new(-debug => $debug);
}

my $tace = $wormbase->tace;

$acedbpath = $wormbase->autoace unless $acedbpath;

my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or die("Connection failure: ". Ace->error);

my $query = 'FIND Strain WHERE Species = "Caenorhabditis elegans"';
my $it = $db->fetch_many(-query => $query);

my $suffix = $test ? '.test.json' : '.json';

$outfile = "./wormbase.agms.${ws_version}.${LINKML_SCHEMA}" . $suffix unless defined $outfile;

my @agms;

while (my $obj = $it->next) {
    unless ($obj->Species) {
	print "No species for $obj - skipping\n";
	next;
    }

    my $dp_xref_dto_json = {
	referenced_curie => 'WB:' . $obj->name,
	page_area => 'strain',
	display_name => $obj->name,
	prefix => 'WB',
	inteternal => JSON::false,
	obsolete => JSON::false
    };
	
    my $data_provider_dto_json = {
	source_organization_abbreviation => 'WB',
	cross_reference_dto => $dp_xref_dto_json,
	internal => JSON::false,
	obsolete => JSON::false
    };
	    
    my $strain = {
	primary_external_id     => "WB:$obj",
	subtype_name      => 'strain',
	taxon_curie       => 'NCBITaxon:' . $obj->Species->NCBITaxonomyID,
	internal          => JSON::false,
	obsolete          => JSON::false,
	data_provider_dto => $data_provider_dto_json
    };

    my ($full_name, $synonyms) = get_name_slot_annotations($obj);
    $strain->{agm_full_name_dto} = $full_name if $full_name;
    $strain->{agm_synonym_dtos} = $synonyms if @$synonyms;
    push @agms, $strain;
}


my $query = 'FIND Genotype WHERE Species = "Caenorhabditis elegans"';
my $it2 = $db->fetch_many(-query => $query);

while (my $obj = $it2->next) {
    unless ($obj->Species) {
	print "No species for $obj - skipping\n";
	next;
    }

    my $dp_xref_dto_json = {
	referenced_curie => 'WB:' . $obj->name,
	page_area => 'genotype',
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

    
    my $genotype = {
	primary_external_id => "WB:$obj",
	subtype_name        => 'genotype',
	taxon_curie         => 'NCBITaxon:' . $obj->Species->NCBITaxonomyID,
	internal            => JSON::false,
	obsolete            => JSON::false,
	data_provider_dto   => $data_provider_dto_json
    };

    my ($full_name, $synonyms) = get_genotype_name_slot_annotations($obj);
    $genotype->{agm_full_name_dto} = $full_name if $full_name;
    $genotype->{agm_synonym_dtos} = $synonyms if @$synonyms;

    push @agms, $genotype;
}
my $data = {
    linkml_version => $LINKML_SCHEMA,
    alliance_member_release_version => $ws_version,
    agm_ingest_set => \@agms,
};

my $json_obj = JSON->new;
my $string   = $json_obj->allow_nonref->canonical->pretty->encode($data);
my $out_fh;
open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";

print $out_fh $string;
$db->close;

sub get_name_slot_annotations {
    my $obj = shift;

    my $full_name;
    if ($obj->Public_name) {
	$full_name = {
	    display_text => $obj->Public_name->name,
	    format_text => $obj->Public_name->name,
	    name_type_name => 'full_name',
	    internal => JSON::false,
	    obsolete => JSON::false
	};
    }
    
    my @synonym_objs = $obj->Other_name;

    my @synonyms;
    if (@synonym_objs) {
	for my $synonym_obj (@synonym_objs) {
	    my $symbol_evidence = get_evidence_curies($synonym_obj);
	    my $synonym = {
		format_text => $synonym_obj->name,
		display_text => $synonym_obj->name,
		name_type_name => "unspecified",
		internal => JSON::false,
		obsolete => JSON::false
	    };
	    $synonym->{evidence_curies} = $symbol_evidence if @$symbol_evidence;
	    push @synonyms, $synonym;
	}
    }
	
    return ($full_name, \@synonyms);
}

sub get_genotype_name_slot_annotations {
    my $obj = shift;

    my $full_name;
    if ($obj->Genotype_name) {
	$full_name = {
	    display_text => $obj->Public_name->name,
	    format_text => $obj->Public_name->name,
	    name_type_name => 'full_name',
	    internal => JSON::false,
	    obsolete => JSON::false
	};
    }
    
    my @synonym_objs = $obj->Genotype_synonym;

    my @synonyms;
    if (@synonym_objs) {
	for my $synonym_obj (@synonym_objs) {
	    my $synonym = {
		format_text => $synonym_obj->name,
		display_text => $synonym_obj->name,
		name_type_name => "unspecified",
		internal => JSON::false,
		obsolete => JSON::false
	    };
	    push @synonyms, $synonym;
	}
    }
	
    return ($full_name, \@synonyms);
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
