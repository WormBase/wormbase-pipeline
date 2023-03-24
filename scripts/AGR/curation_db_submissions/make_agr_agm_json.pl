#!/usr/bin/env perl
use strict;
use Getopt::Long;
use JSON;
use Ace;
use Wormbase;
use Const::Fast;

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($outfile, $acedbpath, $ws_version, $out_fh, $schema);

const my $LINKML_SCHEMA => 'v1.7.0';

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

my $query = 'FIND Strain WHERE Species = "Caenorhabditis elegans"';
my $it = $db->fetch_many(-query => $query);

$outfile = "./wormbase.agms.${ws_version}.${LINKML_SCHEMA}.json" unless defined $outfile;

my @agms;


while (my $obj = $it->next) {
    unless ($obj->Species) {
	print "No species for $obj - skipping\n";
	next;
    }

    my $dp_xref_dto_json = {
	referenced_curie => $obj->name,
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
	curie             => "WB:$obj",
	name              => ($obj->Public_name ? "${\$obj->Public_name}" : "$obj"),
	subtype_name      => 'strain',
	taxon_curie       => 'NCBITaxon:' . $obj->Species->NCBITaxonomyID,
	internal          => JSON::false,
	obsolete          => JSON::false,
	data_provider_dto => $data_provider_dto_json
    };
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
	referenced_curie => $obj->name,
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
	curie             => "WB:$obj",
	name              => ($obj->Genotype_name ? "${\$obj->Genotype_name}" : "$obj"),
	subtype_name      => 'genotype',
	taxon_curie       => 'NCBITaxon:' . $obj->Species->NCBITaxonomyID,
	internal          => JSON::false,
	obsolete          => JSON::false,
	data_provider_dto => $data_provider_dto_json
    };
    push @agms, $genotype;
}
my $data = {
    linkml_version => $LINKML_SCHEMA,
    agm_ingest_set => \@agms,
};

my $json_obj = JSON->new;
my $string   = $json_obj->allow_nonref->canonical->pretty->encode($data);
my $out_fh;
open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";

print $out_fh $string;
$db->close;

