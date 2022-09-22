#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;

use lib $ENV{CVS_DIR};
use Wormbase;
use Modules::AGR;

my ($debug, $test, $verbose, $store, $wormbase);
my ($outfile, $acedbpath, $ws_version, $out_fh, $bgi_json,$disease_file);

GetOptions (
  'debug=s'       => \$debug,
  'test'          => \$test,
  'verbose'       => \$verbose,
  'store:s'       => \$store,
  'database:s'    => \$acedbpath,
  'outfile:s'     => \$outfile,
  'wsversion=s'   => \$ws_version,
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
$outfile    ||= "./wormbase.alleles.${ws_version}.json";

my (@alleles, @annots);
my %objectsIds;

my $db = Ace->connect(-path => $acedbpath, -program => $tace) or die("Connection failure: ". Ace->error);

my @variations = $db->fetch(-query => "Find Variation WHERE species = \"Caenorhabditis elegans\"");
process(\@variations);

my $it = $db->fetch_many(-query => 'find Transgene WHERE Species = "Caenorhabditis elegans"');
process_transgenes($it);

sub process{
    my ($variations) = @_;

    for my $obj (@$variations) {
	next unless $obj->isObject();
	unless ($obj->Species) {
	    print "No species for $obj - skipping\n";
	    next;
	}
	
	my $json_obj = {
	    curie         => "WB:$obj", 
	    taxon       => "NCBITaxon:" . $obj->Species->NCBITaxonomyID,
	    internal      => JSON::false,
	    obsolete   => $obj->Live ? JSON::false : JSON::true
	};
	$json_obj->{symbol} = $obj->Public_name->name if $obj->Public_name;
	if ($obj->Other_name) {
	    my %synonyms = map {$_->name => 1}$obj->Other_name;
	    $json_obj->{synonyms} = [keys \%synonyms];
	}
	push @alleles, $json_obj;
    }
}

sub process_transgenes{
    my ($it) = @_;
    
    while (my $obj = $it->next) {
	next unless $obj->isObject();
	unless ($obj->Species) {
	    print "No species for $obj - skipping\n";
	    next;
	}
	
	my $json_obj = {
	    curie         => "WB:$obj", 
	    taxon         => "NCBITaxon:" . $obj->Species->NCBITaxonomyID,
	    internal      => JSON::false,
	    obsolete      => JSON::false
	};	
	$json_obj->{symbol} = $obj->Public_name->name if $obj->Public_name;
	if ($obj->Synonym) {
	    my %synonyms = map {$_->name => 1}$obj->Synonym;
	    $json_obj->{synonyms} = [keys \%synonyms];
	}
	push @alleles, $json_obj;
    }
}

my $data = {
    allele_ingest_set => \@alleles,
};

open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";

my $json_obj = JSON->new;
my $string = $json_obj->allow_nonref->canonical->pretty->encode($data);
print $out_fh $string;

$db->close;
