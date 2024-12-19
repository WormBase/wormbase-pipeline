#!/usr/bin/env perl
use strict;
use Getopt::Long;
use JSON;
use Ace;
use Wormbase;
use Const::Fast;

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($outfile, $acedbpath, $ws_version, $out_fh, $association_outfile);

const my $LINKML_SCHEMA => 'v2.8.1';

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "test"        => \$test,
	    "verbose"     => \$verbose,
	    "store:s"     => \$store,
	    "database:s"  => \$acedbpath,
	    "outfile:s"   => \$outfile,
	    "assoc:s"     => \$association_outfile,
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

my $query = 'FIND Strain WHERE Species = "Caenorhabditis elegans"';
my $it = $db->fetch_many(-query => $query);

$outfile = "./wormbase.agms.${ws_version}.${LINKML_SCHEMA}.json" unless defined $outfile;
$association_outfile = "./wormbase.agm_associations.${ws_version}.${LINKML_SCHEMA}.json" unless defined $association_outfile;

my @agms;
my @agm_components;

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
	internal => JSON::false,
	obsolete => JSON::false
    };
	
    my $data_provider_dto_json = {
	source_organization_abbreviation => 'WB',
	cross_reference_dto => $dp_xref_dto_json,
	internal => JSON::false,
	obsolete => JSON::false
    };
	    
    my $strain = {
	mod_entity_id     => "WB:$obj",
	name              => ($obj->Public_name ? "${\$obj->Public_name}" : "$obj"),
	subtype_name      => 'strain',
	taxon_curie       => 'NCBITaxon:' . $obj->Species->NCBITaxonomyID,
	internal          => JSON::false,
	obsolete          => JSON::false,
	data_provider_dto => $data_provider_dto_json
    };
    push @agms, $strain;

    if ($obj->Variation) {
	my @component_variations =  $obj->at('Contains.Variation');
	for my $cv (@component_variations) {
	    my $cv_json = {
		agm_subject_identifier => "WB:$obj",
		zygosity_curie => "GENO:0000137",
		allele_identifier => "WB:" . $cv->name,
		relation_name => "contains"
	    };
	    push @agm_components, $cv_json;
	}
    }
    if ($obj->Transgene) {
	my @component_transgenes =  $obj->at('Contains.Transgene');
	for my $ct (@component_transgenes) {
	    my $ct_json = {
		agm_subject_identifier => "WB:$obj",
		zygosity_curie => "GENO:0000137",
		allele_identifier => "WB:" . $ct->name,
		relation_name => "contains"
	    };
	    push @agm_components, $ct_json;
	}
    }
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
	mod_entity_id     => "WB:$obj",
	name              => ($obj->Genotype_name ? "${\$obj->Genotype_name}" : "$obj"),
	subtype_name      => 'genotype',
	taxon_curie       => 'NCBITaxon:' . $obj->Species->NCBITaxonomyID,
	internal          => JSON::false,
	obsolete          => JSON::false,
	data_provider_dto => $data_provider_dto_json
    };
    push @agms, $genotype;

    if ($obj->Variation) {
	my @component_variations =  $obj->at('Genotype_component.Variation');
	for my $cv (@component_variations) {
	    my $cv_json = {
		agm_subject_identifier => "WB:$obj",
		zygosity_curie => "GENO:0000137",
		allele_identifier => "WB:" . $cv->name,
		relation_name => "contains"
	    };
	    push @agm_components, $cv_json;
	}
    }
    if ($obj->Transgene) {
	my @component_transgenes =  $obj->at('Genotype_component.Transgene');
	for my $ct (@component_transgenes) {
	    my $ct_json = {
		agm_subject_identifier => "WB:$obj",
		zygosity_curie => "GENO:0000137",
		allele_identifier => "WB:" . $ct->name,
		relation_name => "contains"
	    };
	    push @agm_components, $ct_json;
	}
    }
}

my $data = {
    linkml_version => $LINKML_SCHEMA,
    alliance_member_release_version => $ws_version,
    agm_ingest_set => \@agms,
};

my $assoc_data = {
    linkml_version => $LINKML_SCHEMA,
    alliance_member_release_version => $ws_version,
    agm_allele_association_ingest_set => \@agm_components
}

my $json_obj = JSON->new;
my $string   = $json_obj->allow_nonref->canonical->pretty->encode($data);
my $out_fh;
open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";

print $out_fh $string;

my $assoc_json_obj = JSON->new;
my $assoc_string = $assoc_json_obj->allow_nonref->canonical->pretty->encode($assoc_data);
my $assoc_out_fh;
open $assoc_out_fh, ">$association_outfile" or die "Could not open $association_outfile for writing\n";

print $assoc_out_fh $assoc_string;
$db->close;

