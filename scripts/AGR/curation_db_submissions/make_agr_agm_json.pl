#!/usr/bin/env perl
use strict;
use Getopt::Long;
use JSON;
use Ace;
use Wormbase;

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($outfile, $acedbpath, $ws_version, $out_fh);

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "test"        => \$test,
	    "verbose"     => \$verbose,
	    "store:s"     => \$store,
	    "database:s"  => \$acedbpath,
	    "outfile:s"   => \$outfile,
            "wsversion=s" => \$ws_version,
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

$outfile = "./wormbase.agms.${ws_version}.json" unless defined $outfile;

my @agms;
while (my $obj = $it->next) {
    unless ($obj->Species) {
	print "No species for $obj - skipping\n";
	next;
    }
    my $strain = {
	curie     => "WB:$obj",
	name      => ($obj->Public_name ? "${\$obj->Public_name}" : "$obj"),
	subtype   => 'strain',
	taxonId   => 'NCBITaxon:' . $obj->Species->NCBITaxonomyID,
	internal  => JSON::false
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
    my $genotype = {
	curie => "WB:$obj",
	name      => ($obj->Genotype_name ? "${\$obj->Genotype_name}" : "$obj"),
	subtype   => 'genotype',
	taxonId   => 'NCBITaxon:' . $obj->Species->NCBITaxonomyID,
	internal  => JSON::false
    };
    push @agms, $genotype;
}
 my $data = {
     agm_ingest_set => \@agms,
};

my $json_obj = JSON->new;
my $string   = $json_obj->allow_nonref->canonical->pretty->encode($data);
my $out_fh;
open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";

print $out_fh $string;
$db->close;

