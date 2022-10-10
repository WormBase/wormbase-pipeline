#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;

use Wormbase;

my ($help, $debug, $test, $verbose, $store, $wormbase, $schema);
my ($outfile, $acedbpath, $ws_version, $out_fh);

GetOptions ("help"        => \$help,
            "debug=s"     => \$debug,
	    "test"        => \$test,
	    "verbose"     => \$verbose,
	    "store:s"     => \$store,
	    "database:s"  => \$acedbpath,
	    "outfile:s"   => \$outfile,
            "wsversion=s" => \$ws_version,
	    "schema=s"    => \$schema
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

    unless ($obj->Species) {
	print "No species for $obj - skipping\n";
	next;
    }

    # Symbol and synonym
    my ($symbol, %synonyms);
    if ($obj->Sequence_name) {
	my $seq_name = $obj->Sequence_name->name;
	if ($obj->CGC_name) {
	    $symbol = $obj->CGC_name->name;
	    $synonyms{$seq_name} = 1;
	} else {
	    $symbol = $seq_name;
	}
    } elsif ($obj->CGC_name) {
	$symbol = $obj->CGC_name->name;
    } else {
	print "No symbol for $obj - skipping\n";
	next;
    }
    map { $synonyms{$_->name} = 1 } $obj->Other_name;
    my @synonym_strings = keys %synonyms;

    # Name (from gene class)
    my ($gene_class_name, $brief_id_name);
    if ($obj->Gene_class and $obj->Gene_class->Description) {
	$gene_class_name = $obj->Gene_class->Description->name;
    }
    if ($obj->Corresponding_CDS and $obj->Corresponding_CDS->Brief_identification) {
	$brief_id_name = $obj->Corresponding_CDS->Brief_identification->name;
    }
    my $name;
    if ($gene_class_name) {
	my ($suff) = $symbol =~ /-(\S+)/; 
	$name = "$gene_class_name $suff";
    } elsif ($brief_id_name) {
	$name = $brief_id_name;
    }

    # Secondary ids
    my @secondary_ids;
    if ($obj->Acquires_merge) {
	foreach my $g ($obj->Acquires_merge) {
	    push @secondary_ids, 'WB:'.$g->name;
	}
    }
    
    my $gene = {
	curie    => 'WB:' . $obj->name,
	symbol   => $symbol,
	taxon    => 'NCBITaxon:' . $obj->Species->NCBITaxonomyID,
	obsolete => $obj->Live ? JSON::false : JSON::true,
	internal => JSON::true
	
    };

    $gene->{name} = $name if $name;
    $gene->{synonyms} = \@synonym_strings if @synonym_strings;
    $gene->{secondary_identifiers} = \@secondary_ids if @secondary_ids;
    
    push @genes, $gene;
}

my $data = {
    linkml_version => $schema,
    gene_ingest_set => \@genes,
};

open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";

my $json_obj = JSON->new;
my $string = $json_obj->allow_nonref->canonical->pretty->encode($data);
print $out_fh $string;

$db->close;

exit(0);
