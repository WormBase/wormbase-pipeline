#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;
use Const::Fast;
use Wormbase;

my ($help, $debug, $test, $verbose, $store, $wormbase, $schema);
my ($outfile, $acedbpath, $ws_version, $out_fh);

const my $LINKML_SCHEMA => 'v1.7.2';

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
    $wormbase = Wormbase->new(-debug => $debug,-test => $test);
}

my $tace = $wormbase->tace;

$acedbpath = $wormbase->autoace unless $acedbpath;

my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or die("Connection failure: ". Ace->error);

my $query = 'FIND Gene WHERE Species = "Caenorhabditis elegans"';

$outfile = "./wormbase.genes.${ws_version}.json" unless defined $outfile;

my @genes;

my $it = $db->fetch_many(-query => $query);

while (my $obj = $it->next) {
    next unless $obj->isObject();
    next unless $obj->name =~ /^WBGene/;
    unless ($obj->Species) {
	print "No species for $obj - skipping\n";
	next;
    }

    my ($symbol, $full_name, $systematic_name, $synonyms) = get_name_slot_annotations($obj);
   
    if (!defined $symbol && $obj->Status && $obj->Status->name eq 'Live') {
	print "Live gene with no symbol $obj - skipping\n";
	next;
    }

    # Secondary ids
    my @secondary_ids;
    if ($obj->Acquires_merge) {
    	foreach my $g ($obj->Acquires_merge) {
	    my $sid_json = {
		secondary_id => 'WB:' . $g->name,
		internal => JSON::false,
		obsolete => JSON::false
	    };
    	    push @secondary_ids, $sid_json;
    	}
    }
    
    my $dp_xref_dto_json = {
	referenced_curie => 'WB:' . $obj->name,
	page_area => 'gene',
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

    my $is_obsolete = $obj->Status && $obj->Status->name eq 'Dead' ? JSON::true : JSON::false;
    my $gene = {
	curie    => 'WB:' . $obj->name,
	gene_symbol_dto   => $symbol,
	taxon_curie    => 'NCBITaxon:' . $obj->Species->NCBITaxonomyID,
	obsolete => $is_obsolete,
	internal => JSON::false,
	created_by_curie => 'WB:curator',
	updated_by_curie => 'WB:curator',
	data_provider_dto => $data_provider_dto_json
    };


    $gene->{gene_full_name_dto} = $full_name if $full_name;
    $gene->{gene_systematic_name_dto} = $systematic_name if $systematic_name;
    $gene->{gene_synonym_dtos} = $synonyms if @$synonyms;
    $gene->{gene_secondary_id_dtos} = \@secondary_ids if @secondary_ids;
    
    push @genes, $gene;
}

my $data = {
    linkml_version => $LINKML_SCHEMA,
    gene_ingest_set => \@genes,
};

open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";

my $json_obj = JSON->new;
my $string = $json_obj->allow_nonref->canonical->pretty->encode($data);
print $out_fh $string;

$db->close;

exit(0);

sub get_name_slot_annotations {
    my $obj = shift;

    my ($symbol_obj, $systematic_name_obj);
    my $symbol_type = "nomenclature_symbol";
    my @synonym_objs = $obj->Other_name;
    if ($obj->CGC_name) {
	$symbol_obj = $obj->CGC_name;
    }
    if ($obj->Sequence_name) {
	unless (defined $symbol_obj) {
	    $symbol_obj = $obj->Sequence_name;
	    $symbol_type = "systematic_name";
	}
	$systematic_name_obj = $obj->Sequence_name;
    }
    if ($obj->Public_name) {
	if (defined $symbol_obj) {
	    push @synonym_objs, $obj->Public_name unless $symbol_obj->name eq $obj->Public_name->name;
	} else {
	    $symbol_obj = $obj->Public_name;
	}
    }

    my $symbol;
    if ($symbol_obj) {
	my $symbol_evidence = get_evidence_curies($symbol_obj);
	$symbol = {
	    display_text => $symbol_obj->name,
	    format_text  => $symbol_obj->name,
	    name_type_name => $symbol_type,
	    internal => JSON::false
	};
	$symbol->{evidence_curies} = $symbol_evidence if @$symbol_evidence;
    } else {
	$symbol = {
	    display_text   => $obj->name,
	    format_text    => $obj->name,
	    name_type_name => 'systematic_name',
	    internal => JSON::false
	};
    }

    my $systematic_name;
    if ($systematic_name_obj) {
	my $systematic_name_evidence = get_evidence_curies($systematic_name_obj);
	$systematic_name = {
	    display_text => $systematic_name_obj->name,
	    format_text => $systematic_name_obj->name,
	    name_type_name => "systematic_name",
	    internal => JSON::false
	};
	$systematic_name->{evidence_curies} = $systematic_name_evidence if @$systematic_name_evidence;
    }

    my @synonyms;
    if (@synonym_objs) {
	for my $synonym_obj (@synonym_objs) {
	    my $symbol_evidence = get_evidence_curies($synonym_obj);
	    my $synonym = {
		format_text => $synonym_obj->name,
		display_text => $synonym_obj->name,
		name_type_name => "unspecified",
		internal => JSON::false
	    };
	    $synonym->{evidence_curies} = $symbol_evidence if @$symbol_evidence;
	    push @synonyms, $synonym;
	}
    }

    my ($full_name, $full_name_evidence);
        
    if ($obj->Gene_class and $obj->Gene_class->Description) {
	my $suffix = '';
	if ($symbol_obj) {
	    if ($symbol_obj->name =~ /-(\S+)$/) {
		$suffix = ' ' . $1;
	    }
	}
	$full_name_evidence = get_evidence_curies($obj->Gene_class->Description);
	$full_name = {
	    display_text => $obj->Gene_class->Description->name . $suffix,
	    format_text => $obj->Gene_class->Description->name . $suffix,
	    name_type_name => "full_name",
	    internal => JSON::false
	};
	$full_name->{evidence_curies} = $full_name_evidence if @$full_name_evidence;
    } elsif ($obj->Corresponding_CDS && $obj->Corresponding_CDS->Brief_identification) {
	$full_name_evidence = get_evidence_curies($obj->Corresponding_CDS->Brief_identification);
	$full_name = {
	    display_text => $obj->Corresponding_CDS->Brief_identification->name,
	    format_text => $obj->Corresponding_CDS->Brief_identification->name,
	    name_type_name => "full_name",
	    internal => JSON::false
	};
	$full_name->{evidence_curies} = $full_name_evidence if @$full_name_evidence;
    }
	
    return ($symbol, $full_name, $systematic_name, \@synonyms);
}

sub get_evidence_curies {
    my $name_obj = shift;

    my @paper_evidence;
    for my $evi ($name_obj->col) {
	next unless $evi->name eq 'Paper_evidence';
	for my $paper ($evi->col) {
	    push @paper_evidence, 'WB:' . $paper->name;
	}
    }

    return \@paper_evidence;
}
