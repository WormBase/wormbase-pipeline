#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Const::Fast;
use Path::Class;
use XML::LibXML;

const my $PURL => 'http://purl.obolibrary.org/obo/chebi.owl';

my ($outfile);

GetOptions (
    "outfile|o=s" => \$outfile
    ) || die ("Unknown comand line option: $@\n");

if (!$outfile) {
    $outfile = './chebi_id_to_name_map.txt';
}

my $out_fh = file($outfile)->openw;


my $chebi_dom = XML::LibXML->load_xml(location => $PURL);
for my $class ($chebi_dom->findnodes('//owl:Class')) {
    my $id = $class->findvalue('./oboInOwl:id');
    my $name = $class->findvalue('./rdfs:label');
    my $deprecated = $class->findvalue('./owl:deprecated');
    if ($deprecated ne 'true') {
	if (defined $id and defined $name) {
	    $out_fh->print("$id\t$name\n");
	}
	else {
	    warn ("Could not parse name corresponding to $id\n");
	}
    }
}



    
