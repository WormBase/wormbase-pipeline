#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;
use lib $ENV{CVS_DIR};
use Wormbase;
use Modules::AGR;
use Path::Class;
use Const::Fast;

const my $LINKML_SCHEMA => 'v2.13.0';

my ($debug, $test, $verbose, $store, $wormbase, $acedbpath, $ws_version, $outfile, $schema, $dates_file);

GetOptions (
    'debug=s'      => \$debug,
    'test'         => \$test,
    'verbose'      => \$verbose,
    'store:s'      => \$store,
    'database:s'   => \$acedbpath,
    'outfile:s'    => \$outfile,
    'wsversion=s'  => \$ws_version
    )||die("unknown command line option: $@\n");

if ( $store ) {
    $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
    $wormbase = Wormbase->new( -debug   => $debug,
			       -test    => $test
	);
}

my $tace = $wormbase->tace;

$acedbpath = $wormbase->autoace unless $acedbpath;
$ws_version = $wormbase->get_wormbase_version_name unless $ws_version;


$outfile = "./wormbase.phenotype_annotations.${ws_version}.${LINKML_SCHEMA}.json" unless defined $outfile;


my %zeco = (
    'experimental conditions' => 'ZECO:0000104',
    'chemical treatment'      => 'ZECO:0000111',
    'temperature exposure'    => 'ZECO:0000160',
    );

my $db = Ace->connect(-path => $acedbpath, -program => $tace) or die('Connection failure: '. Ace->error);

my $data_provider_dto_json = {
    source_organization_abbreviation => 'WB',
    internal => JSON::false,
    obsolete => JSON::false
};

my (@allele_annotations, @agm_annotations);

my $it = $db->fetch_many(-query => 'find Variation WHERE Live AND COUNT(Gene) < 2 AND (Phenotype OR Phenotype_not_observed)  AND NOT Natural_variant');
my $negated;
while(my $obj = $it->next) {
    $negated = 0;
    if ($obj->Phenotype) {
	for my $phenotype ($obj->Phenotype) {
	    process_phenotypes($obj, $phenotype, $negated);
	}
    }
    $negated = 1;
    if ($obj->Phenotype_not_observed) {
	for my $phenotype_not_observed ($obj->Phenotype_not_observed) {
	    process_phenotypes($obj, $phenotype_not_observed, $negated);
	}
    }
}

$it = $db->fetch_many(-query => 'find Transgene WHERE Phenotype OR Phenotype_not_observed');
while(my $obj = $it->next) {
    $negated = 0;
    if ($obj->Phenotype) {
	for my $phenotype ($obj->Phenotype) {
	    process_phenotypes($obj, $phenotype, $negated);
	}
    }
    $negated = 1;
    if ($obj->Phenotype_not_observed) {
	for my $phenotype_not_observed ($obj->Phenotype_not_observed) {
	    process_phenotypes($obj, $phenotype_not_observed, $negated);
	}
    }
}

$it = $db->fetch_many(-query => 'find Strain WHERE Phenotype OR Phenotype_not_observed');
while(my $obj = $it->next) {
    $negated = 0;
    if ($obj->Phenotype) {
	for my $phenotype ($obj->Phenotype) {
	    process_phenotypes($obj, $phenotype, $negated);
	}
    }
    $negated = 1;
    if ($obj->Phenotype_not_observed) {
	for my $phenotype_not_observed ($obj->Phenotype_not_observed) {
	    process_phenotypes($obj, $phenotype_not_observed, $negated);
	}
    }
}

$db->close;

my $all_annots = {
    linkml_version => $LINKML_SCHEMA,
    alliance_member_release_version => $ws_version
};
$all_annots->{phenotype_allele_ingest_set} = \@allele_annotations;
$all_annots->{phenotype_agm_ingest_set} = \@agm_annotations;

print_json($outfile, $all_annots);
    
exit(0);


sub process_phenotypes {

    my ($obj, $phenotype, $negated) = @_;
    
    unless ($phenotype->Primary_name) {
	print STDERR "No primary name for $phenotype - skipping\n";
	next;
    }
	    
    my %papers;
    my @inferred_genes;
    foreach my $evi ($phenotype->col()) {
	if ($evi->name eq 'Paper_evidence') {
	    my $ref_count = 0;
	    foreach my $wb_paper ($evi->col ) {
		$papers{$wb_paper} = get_paper_id($wb_paper);
	    }
	} elsif($obj->name =~ /WBTransgene/ && $evi->name eq 'Caused_by_gene'){
	    foreach my $g ($evi->col){
		push @inferred_genes, "WB:$g";
	    }
	}
    }
    
    if ($obj->name =~ /WBVar/ && $obj->Gene) {
	for my $gene ($obj->Gene) {
	    push @inferred_genes, "WB:" . $gene->name;
	}
    }
    
    if (scalar @inferred_genes > 1) {
	print STDERR "Multiple inferred genes for $obj - using first\n";
    }
    
    for my $paper (keys %papers) {
	my $annot = {
	    mod_internal_id            => "WB:$obj|WB:" . $phenotype->name . "|" . $papers{$paper},
	    data_provider_dto          => $data_provider_dto_json,
	    phenotype_statement        => $phenotype->Primary_name->name,
	    phenotype_term_curies      => [$phenotype->name],
	    evidence_curie             => $papers{$paper},
	    negated                    => $negated ? JSON::true : JSON::false,
	    internal                   => JSON::false,
	    obsolete                   => JSON::false
	};

	if ($obj->name =~ /WBStrain/) {
	    $annot->{agm_identifier} = "WB:$obj";
	} else {
	    $annot->{allele_identifier} => "WB:$obj";
	}
	
	if (@inferred_genes) {
	    $annot->{inferred_gene_identifier} = $inferred_genes[0];
	}
	
	my @condition_relations = @{get_condition_relations($phenotype, $paper)};
	$annot->{condition_relation_dtos} = \@condition_relations if @condition_relations;

	if ($obj->name =~ /WBVar/ && $obj->Phenotype_remark) {
	    my @notes;
	    for my $remark ($obj->Phenotype_remark) {
		push @notes, {
		    note_type_name => 'remark',
		    internal       => JSON::false,
		    free_text      => $obj->Phenotype_remark->name
		};
	    }
	    $annot->{note_dtos} = \@notes;
	}
    
	if ($obj->name =~ /WBStrain/) {
	    push @agm_annotations, $annot;
	} else {
	    push @allele_annotations, $annot;
	}
    }
}    

##############################################

sub print_json {
    my ($outfile, $data) = @_;

    open(my $out_fh, ">$outfile") or die("cannot open $outfile : $!\n");  
    my $json_obj = JSON->new;
    my $string = $json_obj->allow_nonref->canonical->pretty->encode($data);
    print $out_fh $string;
    close $out_fh;

    return;
}

sub get_paper_id {
    my ($wb_paper) = @_;
 
    my $pmid;
    foreach my $db ($wb_paper->Database) {
    	if ($db->name eq 'MEDLINE') {
	        $pmid = $db->right->right->name;
	        last;
	    }
    }
    my $paper_id  = $pmid ? "PMID:$pmid" : "WB:$wb_paper";
    
    return $paper_id;
}

sub get_chemical_ontology_id {
    my $obj = shift;

    for my $db ($obj->Database) {
        next unless $db->name eq 'ChEBI';
        return 'CHEBI:' . $db->right->right->name;
    }

    return 'WB:' . $obj->name;
}


sub get_condition_relations {
    my ($obj, $paper_id) = @_;
    
    my $condition_relation_type = 'has_condition';
    
    
    my (@conditions, @assays_and_molecules);
    my $pa = $obj->at('Phenotype_assay');
    if (defined $pa){
        for my $temp ($pa->at('Temperature')) {
            my $paper = $temp->at('Paper_evidence')->at() if $temp->at('Paper_evidence');
            next unless defined $paper;
            push @assays_and_molecules, {
                condition_free_text => $temp->name,
                condition_class_curie => $zeco{'temperature exposure'},
            } if $paper->name eq $paper_id;
        }
        for my $treatment ($pa->at('Treatment')) {
            my $paper = $treatment->at('Paper_evidence')->at() if $treatment->at('Paper_evidence');
            next unless defined $paper;
            push @assays_and_molecules, {
                condition_free_text => $treatment->name,
                condition_class_curie => $zeco{'experimental conditions'},
            } if $paper->name eq $paper_id;
        }
    }

    my $ab = $obj->at('Affected_by');
    if (defined $ab) {
	for my $mol ($ab->at('Molecule')) {
	    my $paper = $mol->at('Paper_evidence')->at() if $mol->at('Paper_evidence');
	    next unless defined $paper;
	    push @assays_and_molecules, {
		condition_chemical_curie => get_chemical_ontology_id($mol),
		condition_class_curie => $zeco{'chemical treatment'},
	    } if $paper eq $paper_id;
	}
    }
   
    push @conditions, {
        condition_relation_type_name => $condition_relation_type,
        conditions                   => \@assays_and_molecules
    } if @assays_and_molecules;
   
    return \@conditions;
}


sub get_chemical_ontology_id {
    my $obj = shift;
    
    for my $db ($obj->Database) {
        next unless $db->name eq 'ChEBI';
        return 'CHEBI:' . $db->right->right->name;
    }
    
    return 'WB:' . $obj->name;
}

1;
