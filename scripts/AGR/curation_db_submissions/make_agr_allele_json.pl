#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;

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

my (@alleles, @annots);
my %objectsIds;

my $db = Ace->connect(-path => $acedbpath, -program => $tace) or die("Connection failure: ". Ace->error);

my $it = $db->fetch_many(-query => "Find Variation WHERE species = \"Caenorhabditis elegans\"");
process($it);

my $tgit = $db->fetch_many(-query => 'find Transgene WHERE Species = "Caenorhabditis elegans"');
process_transgenes($tgit);


my $data = {
    linkml_version    => $schema,
    allele_ingest_set => \@alleles,
};

open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";

my $json_obj = JSON->new;
my $string = $json_obj->allow_nonref->canonical->pretty->encode($data);
print $out_fh $string;

$db->close;

sub process {
    my ($it) = @_;

    my $var_count = 0;
    while (my $obj = $it->next) {
	$var_count++;
	next unless $obj->isObject();
	unless ($obj->Species) {
	    print "No species for $obj - skipping\n";
	    next;
	}
	
	my $json_obj = {
	    curie    => "WB:" . $obj->name,
	    name     => $obj->Public_name ? $obj->Public_name->name : $obj->name,
	    taxon       => "NCBITaxon:" . $obj->Species->NCBITaxonomyID->name,
	    internal      => JSON::false,
	    obsolete   => $obj->Live ? JSON::false : JSON::true,
	    created_by => 'WB:curator',
	    updated_by => 'WB:curator'
	};
	$json_obj->{symbol} = $obj->Public_name->name if $obj->Public_name;
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
		$collection =~ s/_/ /g;
		$collection .= ' project' if $collection eq 'Million mutation';
		$json_obj->{in_collection} = $collection;
	    }
	}
	$json_obj->{sequencing_status} = lc $obj->SeqStatus->name if $obj->SeqStatus;
	if ($obj->Other_name) {
	    my @synonyms_to_submit;
	    my %synonyms = map {$_->name => 1} $obj->Other_name;
	    for my $syn (keys %synonyms) {
		push @synonyms_to_submit, {name => $syn} unless $syn =~ /^[^:]+:[pcg]\./; # remove HGVS identifiers (already generated in Alliance)
	    }
		
	    $json_obj->{synonyms} = \@synonyms_to_submit if @synonyms_to_submit > 0;
	}
	if ($obj->Reference) {
	    $json_obj->{references} = get_papers($obj);
	}

	# Stick some random data in for curation DB test files
	if ($curation_test) {
	    my @inheritance_modes = qw(dominant recessive semi-dominant unknown);
	    $json_obj->{inheritance_mode} = $inheritance_modes[rand @inheritance_modes];
	    $json_obj->{is_extinct} = $var_count % 2 == 0 ? JSON::true : JSON::false;
	    $json_obj->{date_updated} = get_random_datetime(10);
	    $json_obj->{date_created} = get_random_datetime(1);
	}
	
	push @alleles, $json_obj;
	last if $limit && $limit == $var_count;
    }
}

sub process_transgenes {
    my ($it) = @_;

    my $var_count = 0;
    while (my $obj = $it->next) {
	$var_count++;
	next unless $obj->isObject();
	unless ($obj->Species) {
	    print "No species for $obj - skipping\n";
	    next;
	}
	
	my $json_obj = {
	    curie         => "WB:" . $obj->name, 
	    name          => $obj->Public_name ? $obj->Public_name->name : $obj->name,
	    taxon         => "NCBITaxon:" . $obj->Species->NCBITaxonomyID->name,
	    internal      => JSON::false,
	    obsolete      => JSON::false
	};	
	$json_obj->{symbol} = $obj->Public_name->name if $obj->Public_name;
	if ($obj->Synonym) {
	    my %synonyms = map {$_->name => 1}$obj->Synonym;
	    if (scalar keys %synonyms > 0) {
		my @synonyms_to_submit;
		for my $syn (keys %synonyms) {
		    push @synonyms_to_submit, {name => $syn};
		}
		$json_obj->{synonyms} = \@synonyms_to_submit;
	    }
	}
	if ($obj->Reference) {
	    $json_obj->{references} = get_papers($obj);
	}

	
	# Stick some random data in for curation DB test files
	if ($curation_test) {
	    my @inheritance_modes = qw(dominant recessive semi-dominant unknown);
	    $json_obj->{inheritance_mode} = $inheritance_modes[rand @inheritance_modes];
	    $json_obj->{is_extinct} = $var_count % 2 == 0 ? JSON::true : JSON::false;
	    $json_obj->{date_updated} = get_random_datetime(10);
	    $json_obj->{date_created} = get_random_datetime(1);
	}
	
	push @alleles, $json_obj;
    }
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
