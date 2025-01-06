#!/usr/bin/env perl
use strict;
use Getopt::Long;
use JSON;
use Ace;
use Wormbase;
use Const::Fast;
use Path::Class;
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

my $query = 'FIND Variation WHERE Species = "Caenorhabditis elegans"';

$outfile = "./wormbase.agm_associations.${ws_version}.${LINKML_SCHEMA}.json" unless defined $association_outfile;

open(ASSOC_OUT, '>',$outfile) or die "Could not open $outfile for writing\n";
print ASSOC_OUT "{\n   \"linkml_version\" : \"" . $LINKML_SCHEMA . "\",\n    \"alliance_member_release_version\" : \"" . $ws_version . "\",\n";
print ASSOC_OUT "    \"agm_allele_association_ingest_set\" : [";

print "FETCHING VARIANTS\n";
my @agms;
my $it = $db->fetch_many(-query =>$query);
my $allele_assoc_count = 0;
print "PROCESSING VARIANTS\n";
while (my $obj = $it->next) {

    if ($obj->Strain) {
	for my $strain ($obj->Strain) {
	    if ($allele_assoc_count > 0) {
		print ASSOC_OUT ',';
	    }
	    $allele_assoc_count++;
	    print ASSOC_OUT '{"agm_subject_identifier": "WB:' . $strain . '", "zygosity_curie": "GENO:0000137", "allele_identifier": "WB:' . $obj . '", "relation_name": "contains"}' . "\n";
	}
    }
	    
    if ($obj->Component_of_genotype) {
	for my $gc ($obj->Component_of_genotype) {
	    if ($allele_assoc_count > 0) {
		print ASSOC_OUT ',';
	    }
	    $allele_assoc_count++;
	    print ASSOC_OUT '{"agm_subject_identifier": "WB:' . $gc . '", "zygosity_curie": "GENO:0000137", "allele_identifier": "WB:' . $obj . '", "relation_name": "contains"}' . "\n"; 
	}
    }
}


my $query = 'FIND Transgene WHERE Species = "Caenorhabditis elegans"';
print "FETCHING TRANSGENES\n";
my $git = $db->fetch_many(-query => $query);
print "PROCESSING TRANSGENES\n";
while (my $obj = $git->next) {

    if ($obj->Strain) {
	for my $strain ($obj->Strain) {
	    if ($allele_assoc_count > 0) {
		print ASSOC_OUT ',';
	    }
	    $allele_assoc_count++;
	    print ASSOC_OUT '{"agm_subject_identifier": "WB:' . $strain . '", "zygosity_curie": "GENO:0000137", "allele_identifier": "WB:' . $obj . '", "relation_name": "contains"}' . "\n";
	}
    }
	    
    if ($obj->Component_of_genotype) {
	for my $gc ($obj->Component_of_genotype) {
	    if ($allele_assoc_count > 0) {
		print ASSOC_OUT ',';
	    }
	    $allele_assoc_count++;
	    print ASSOC_OUT '{"agm_subject_identifier": "WB:' . $gc . '", "zygosity_curie": "GENO:0000137", "allele_identifier": "WB:' . $obj . '", "relation_name": "contains"}' . "\n"; 
	}
    }
}
print ASSOC_OUT "\n    ]";
print ASSOC_OUT "\n}\n";
print "DONE\n";

close(ASSOC_OUT);
$db->close;
