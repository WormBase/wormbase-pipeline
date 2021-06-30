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
my ($outfile, $acedbpath, $ws_version, $out_fh, $bgi_json);

GetOptions (
    'debug=s'     => \$debug,
    'test'        => \$test,
    'verbose'     => \$verbose,
    'store:s'     => \$store,
    'database:s'  => \$acedbpath,
    'outfile:s'   => \$outfile,
    'wsversion=s' => \$ws_version,
    'bgijson=s'   => \$bgi_json,
    )||die(@!);

if ($store) {
    $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
    $wormbase = Wormbase->new(-debug => $debug,-test => $test);
}

my $tace      = $wormbase->tace;
my $date      = AGR::get_rfc_date();
my $alt_date  = join("/", $date =~ /^(\d{4})(\d{2})(\d{2})/);
my $taxid     = $wormbase->ncbi_tax_id;
my $full_name = $wormbase->full_name;

$acedbpath  = $wormbase->autoace unless $acedbpath;
$ws_version = $wormbase->get_wormbase_version_name unless $ws_version;

if (not defined $outfile) {
    $outfile = "./wormbase.agr_phenotype.${ws_version}.json";
}

my %zeco = (
    'experimental conditions' => 'ZECO:0000104',
    'chemical treatment'      => 'ZECO:0000111',
    'temperature exposure'    => 'ZECO:0000160',
    );


# restrict to genes in the BGI, if provded
my ($bgi_genes, @pheno_annots);

$bgi_genes = AGR::get_bgi_genes( $bgi_json ) if defined $bgi_json;

my $db = Ace->connect(-path => $acedbpath, -program => $tace) or die('Connection failure: '. Ace->error);

my $it = $db->fetch_many(-query => 'find Variation WHERE Live AND COUNT(Gene) == 1 AND Phenotype AND NOT Natural_variant');
process($it);

$it = $db->fetch_many(-query => 'find Transgene WHERE Phenotype');
process($it);

# create implicit mappings for BGI genes from linked RNAi and Variations
my @genes = keys %$bgi_genes;
map {s/WB://} @genes;

foreach my $g (@genes){
    my $gene = $db->fetch(Gene => $g);
    next unless $gene;
    my @objects;
#	push @objects, $gene->RNAi_result if $gene->RNAi_result;
    push @objects, $gene->Allele if $gene->Allele;
    process_genes_phenotype($g,@objects) if $objects[0];
}

my $data = {
    metaData => AGR::get_file_metadata_json( (defined $ws_version) ? $ws_version : $wormbase->get_wormbase_version_name(), $date ),
    data     => \@pheno_annots,
};

if (defined $outfile) {
    open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";
} else {
    $out_fh = \*STDOUT;
}

my $json_obj = JSON->new;
my $string   = $json_obj->allow_nonref->canonical->pretty->encode($data);
print $out_fh $string;

$db->close;

exit(0);

##############################################

# bit to process linked phenotypes
sub process_genes_phenotype{
    my $gen = shift @_;
    my @secondaries = @_;
    
    foreach my $obj(@secondaries){
	next unless $obj->Phenotype;
	my @linked_genes = $obj->Gene;
	next if scalar(@linked_genes)>1; # to remove variants with more than one gene connection
	my @phenotypes = $obj->Phenotype;
	
	foreach my $pt (@phenotypes){
	    my $phen_id   = 'WB:'.$pt->name;
	    my $phen_desc = $pt->Primary_name->name;
	    my %papers;
	    
	    foreach my $evi ($pt->col()) {
		if ($evi->name eq 'Paper_evidence') {
		    foreach my $wb_paper ($evi->col ) {
			$papers{$wb_paper} = &get_paper_json($wb_paper);
		    }
		}
	    }
	    
	    # if ($obj->class eq 'RNAi' && $obj->Reference){
	    #  foreach my $wb_paper ($obj->Reference) {
	    #     push @paper, &get_paper_json($wb_paper);
	    #  }
	    
	    foreach my $paper (keys %papers) {
		my $json_obj = {
        	    objectId                 => "WB:$gen",
	            primaryGeneticEntityIDs  => ["WB:$obj"],
        	    phenotypeTermIdentifiers => [ { termId => $phen_id, termOrder => 1 } ],
	            phenotypeStatement       => $phen_desc,
        	    dateAssigned             => $date,
	            evidence                 => $papers{$paper},
		}; 
		#my @condition_relations = @{get_condition_relations($pt, $paper)};
		#$json_obj->{conditionRelations} = \@condition_relations if @condition_relations;
		push @pheno_annots, $json_obj;
	    }
	}
    }
}




# bit to process transgenes and variations
sub process {
    my ($it)=@_;
    while (my $obj = $it->next) {
	next unless $obj->isObject();
	
        if ($obj->name=~/WBVar/){
	    my $gene = $obj->Gene->name;
	    next unless $bgi_genes->{"WB:$gene"};
	}
	
	my @phenotypes = $obj->Phenotype;
	
        foreach my $pt (@phenotypes){
	    my $phen_id   = 'WB:'.$pt->name;
	    my $phen_desc = $pt->Primary_name->name;
	    my %papers;
	    my @caused_by_genes;
	    
	    foreach my $evi ($pt->col()) {
		if ($evi->name eq 'Paper_evidence') {
		    foreach my $wb_paper ($evi->col ) {
			$papers{$wb_paper} = &get_paper_json($wb_paper);
		    }
		} elsif($evi->name eq 'Caused_by_gene' && $obj->name =~ /WBTransgene/){
		    foreach my $g ($evi->col){
			push @caused_by_genes, "$g";
		    }
		}
	    }
	    
	    foreach my $paper (keys %papers) {
		my $json_obj = {
		    objectId                 => "WB:$obj",
		    phenotypeTermIdentifiers => [ { termId => $phen_id, termOrder => 1 } ],
		    phenotypeStatement       => $phen_desc,
		    dateAssigned             => $date,
		    evidence                 => $papers{$paper},
		}; 
		my @condition_relations = @{get_condition_relations($pt, $paper)};
		$json_obj->{conditionRelations} = \@condition_relations if @condition_relations;
		push @pheno_annots, $json_obj;
		
		foreach my $g (@caused_by_genes) {
		    my $json_obj = {
			objectId                 => "WB:$g",
			primaryGeneticEntityIDs  => ["WB:$obj"],
			phenotypeTermIdentifiers => [ { termId => $phen_id, termOrder => 1 } ],
			phenotypeStatement       => $phen_desc,
			dateAssigned             => $date,
			evidence                 => $papers{$paper},
		    };
		 #   my @condition_relations = @{get_condition_relations($pt, $paper)};
		 #   $json_obj->{conditionRelations} = \@condition_relations if @condition_relations;
		    push @pheno_annots, $json_obj;
		}
		
	    }
	}
    }
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
    
    
    my @conditions;
    my $pa = $obj->at('Phenotype_assay');
    if (defined $pa){
	my @assays;
	for my $temp ($pa->at('Temperature')) {
	    my $paper = $temp->at('Paper_evidence')->at() if $temp->at('Paper_evidence');
	    next unless defined $paper;
	    push @assays, {
		conditionStatement => 'temperature exposure:' . $temp->name,
		conditionClassId => $zeco{'temperature exposure'},
	    } if $paper->name eq $paper_id;
	}
	for my $treatment ($pa->at('Treatment')) {
	    my $paper = $treatment->at('Paper_evidence')->at() if $treatment->at('Paper_evidence');
	    next unless defined $paper;
	    push @assays, {
		conditionStatement => 'experimental conditions:' . $treatment->name,
		conditionClassId => $zeco{'experimental conditions'},
            } if $paper->name eq $paper_id;
	}
        push @conditions, {conditionRelationType => $condition_relation_type,
                           conditions            => \@assays} if @assays;
    }
    
    my $ab = $obj->at('Affected_by');
    if (defined $ab) {
	my @molecules;
	for my $mol ($ab->at('Molecule')) {
	    my $paper = $mol->at('Paper_evidence')->at() if $mol->at('Paper_evidence');
	    next unless defined $paper;
	    push @molecules, {
		conditionStatement => 'chemical treatment:' . $mol->Public_name->name,
		chemicalOntologyId => get_chemical_ontology_id($mol),
		conditionClassId => $zeco{'chemical treatment'},
		conditionId => $zeco{'chemical treatment'},
	    } if $paper eq $paper_id;
	}
        push @conditions, {conditionRelationType => $condition_relation_type,
                           conditions            => \@molecules} if @molecules;
    }
    
    return \@conditions;
}


sub get_paper_json {
    my ($wb_paper) = @_;
    my $json_paper = {};
    
    my $pmid;
    foreach my $db ($wb_paper->Database) {
	if ($db->name eq 'MEDLINE') {
	    $pmid = $db->right->right->name;
	    last;
	}
    }
    $json_paper->{publicationId}  = $pmid ? "PMID:$pmid" : "WB:$wb_paper";
    $json_paper->{crossReference} = {id =>"WB:$wb_paper",pages => ['reference']};
    
    return $json_paper;
}
