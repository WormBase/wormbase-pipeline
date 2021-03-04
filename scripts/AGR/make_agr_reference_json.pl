#!/usr/bin/env perl

use strict;
use Storable;
use Getopt::Long;
use Ace;
use JSON;
use Wormbase;
use Modules::AGR;

my %WB_AGR_TYPE_MAP = (
    Journal_article           => 'Research Article',
    Review                    => 'Review Article',
    Comment                   => 'Other',
    News                      => 'Other',
    Letter                    => 'Other',
    Editorial                 => 'Other',
    Congresses                => 'Other',
    Historical_article        => 'Other',
    Biography                 => 'Other',
    Interview                 => 'Other',
    Lectures                  => 'Other',
    Interactive_tutorial      => 'Other',
    Retraction_of_publication => 'Retraction',
    Retracted_publication     => 'Other',
    Technical_report          => 'Research Article',
    Directory                 => 'Other',
    Monograph                 => 'Book',
    Published_erratum         => 'Other',
    Meeting_abstract          => 'Conference Publication',
    Gazette_article           => 'Other',
    Book_chapter              => 'Book',
    Book                      => 'Book',
    Email                     => 'Personal Communication',
    WormBook                  => 'Book',
    Other                     => 'Other',
    Micropublication          => 'Research Article',
    Method                    => 'Research Article',
    );

my ($debug, $test, $verbose, $store, $wormbase, $build);
my ($outfile, $re_outfile, $acedbpath, $ws_version);

GetOptions (
    'debug=s'     => \$debug,
    'test'        => \$test,
    'verbose'     => \$verbose,
    'store:s'     => \$store,
    'database:s'  => \$acedbpath,
    'ref:s'       => \$outfile,
    'refex:s'     => \$re_outfile,
    'wsversion=s' => \$ws_version,
)||die("unknown command line option: $@\n");

if ( $store ) {
    $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
    $wormbase = Wormbase->new( -debug   => $debug,
			       -test    => $test,
	);
}

my $tace = $wormbase->tace;
$acedbpath = $wormbase->autoace unless $acedbpath;
$ws_version = $wormbase->get_wormbase_version_name unless $ws_version;
$outfile = "./reference.${ws_version}.json" unless defined $outfile;
$re_outfile = "./reference_exchange.${ws_version}.json" unless defined $re_outfile;


my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or die("Connection failure: ". Ace->error);

my $it = $db->fetch_many(-query => 'find Paper WHERE Valid');
my @references;
my @exchange_refs;
while (my $obj = $it->next) {
    next unless $obj->Species eq 'Caenorhabditis elegans';

    my $ref_id = "WB:$obj";
    if ($obj->Database) {
	for my $db ($obj->Database) {
	    $ref_id = $db->right->name . ':' . $db->right->right->name;
	    last if $db->right->name eq 'MEDLINE';
	}
    }
  
    # Get book chapter publication dates from book
    my $publication_date;
    if ($obj->Publication_date) {
	$publication_date = $obj->Publication_date->name;
    }
    elsif ($obj->Contained_in and $obj->Contained_in->Publication_date) {
	$publication_date = $obj->Contained_in->Publication_date->name;
    }
    else {
	next;
    }

    next unless ($obj->Title and $obj->Brief_citation);

    my $is_pubmed = $ref_id =~ /^PMID:/ ? 1 : 0;
    my ($reference, $exchange_ref);

    my $author_list = get_author_list($obj, $ref_id);
    my $resource_title = $obj->Journal ? $obj->Journal->name : $obj->Title->name;
    $reference = {
	primaryId            => $ref_id,
	authors              => $author_list,
	title                => $obj->Title->name,
	datePublished        => $publication_date,
	citation             => $obj->Brief_citation->name,
	tags                 => [{referenceId => $ref_id, tagName => 'inCorpus', tagSource => 'WB'}],
	resourceAbbreviation => $resource_title,
    };

    $exchange_ref = {
	pubMedId  => $ref_id,
	modId     => 'WB:' . $obj->name,
	tags      => [{referenceId => $ref_id, tagName => 'inCorpus', tagSource => 'WB'}],
    } if $is_pubmed;
	
    
    if ($obj->Type) {
	my $wb_article_type = $obj->Type->name;
	$reference->{allianceCategory} = $WB_AGR_TYPE_MAP{$wb_article_type};
	$exchange_ref->{allianceCategory} = $WB_AGR_TYPE_MAP{$wb_article_type} if $is_pubmed;
	$reference->{MODReferenceTypes} = [{referenceType => $wb_article_type, source => 'WB'}];
	$exchange_ref->{MODReferenceTypes} = [{referenceType => $wb_article_type, source => 'WB'}] if $is_pubmed;
    }
    else {
	$reference->{allianceCategory} = 'Unknown';
	$exchange_ref->{allianceCategory} = 'Unknown' if $is_pubmed;
    }

    my @keywords;
    for my $kw ($obj->at('Keyword')) {
	push @keywords, $kw->name;
    };
    $reference->{keywords}  = \@keywords if @keywords;
    $reference->{volume}    = $obj->Volume->name if $obj->Volume;
    $reference->{pages}     = $obj->Page->name if $obj->Page;
    $reference->{publisher} = $obj->Publisher->name if $obj->Publisher;
    
    my @xrefs;
    push @xrefs, {id => "WB:$obj", pages => ['reference']};

    for my $name ($obj->Name) {
	next unless $name =~ /^doi(10\..+)$/;
	push @xrefs, {id => 'DOI:' . $1};
        last;
    }

    $reference->{crossReferences} = \@xrefs if @xrefs;

    push @references, $reference;
    push @exchange_refs, $exchange_ref if $is_pubmed;
}

open (OUT, '>', $outfile) or die "Cannot open $outfile to write: $!\n";
my $data = {
    metaData => AGR::get_file_metadata_json($ws_version),
    data     => \@references,
};

my $json_obj = JSON->new;
print OUT $json_obj->allow_nonref->canonical->pretty->encode($data);
close(OUT);

open (RE_OUT, '>', $re_outfile) or die "Cannot open $outfile to write: $!\n";
my $re_data = {
    metaData => AGR::get_file_metadata_json($ws_version),
    data     => \@exchange_refs,
};

my $re_json_obj = JSON->new;
print RE_OUT $re_json_obj->allow_nonref->canonical->pretty->encode($re_data);
close(RE_OUT);


exit(0);


sub get_author_list {
    my ($obj, $ref_id) = @_;

    my $author_data = [];
    my ($first_name, $middle_name, $last_name);

    for my $author ($obj->Author) {
        my $nr_people = 1;
        if ($author->right) {
	    my $author_str = $author->name;
	    $author_str =~ s/\./\\\./g;
            my @people = $obj->at('Author')->at($author_str.'[2]');
            $nr_people = @people;
        }
       
	for (1 .. $nr_people) {
	    push @$author_data, {name => $author->name, referenceId => $ref_id};
	}
    }

    return $author_data;
}
