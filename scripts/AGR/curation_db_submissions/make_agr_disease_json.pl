#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;
use Storable qw(dclone);
use DateTime;
use lib $ENV{CVS_DIR};
use Wormbase;
use Modules::AGR;
use lib "$ENV{CVS_DIR}/ONTOLOGY";
use Path::Class;
use Const::Fast;
use XML::LibXML;

const my $LINKML_SCHEMA = 'v1.5.0';
const my $CHEBI_PURL => 'http://purl.obolibrary.org/obo/chebi.owl';

my ($debug, $test, $verbose, $store, $wormbase, $acedbpath, $ws_version, $outfile, $schema);

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
my $date = AGR::get_rfc_date();
my $alt_date = join("/", $date =~ /^(\d{4})(\d{2})(\d{2})/);

my $chebi_name_map = get_chebi_name_map();

$acedbpath = $wormbase->autoace unless $acedbpath;
$ws_version = $wormbase->get_wormbase_version_name unless $ws_version;

$outfile = "./wormbase.disease_annotations.${ws_version}.json" unless defined $outfile;


my %go2eco = (
    IMP => 'ECO:0000315',
    IEA => 'ECO:0000265',
    ISS => 'ECO:0000250',
    ND  => 'ECO:0000307',
    IPI => 'ECO:0000021',
    EXP => 'ECO:0000269',
    IDA => 'ECO:0000314',
    IGI => 'ECO:0000316',
    IEP => 'ECO:0000270',
    TAS => 'ECO:0000304',
    NAS => 'ECO:0000303',
    IC  => 'ECO:0000305',
    ISO => 'ECO:0000266',
    ISA => 'ECO:0000247',
    ISM => 'ECO:0000255',
    IGC => 'ECO:0000317',
    IBA => 'ECO:0000318',
    IBD => 'ECO:0000319',
    IKR => 'ECO:0000320',
    IRD => 'ECO:0000321',
    RCA => 'ECO:0000245',
    IMR => 'ECO:0000320',
    );

my %zeco = (
    'experimental conditions' => 'ZECO:0000104',
    'chemical treatment'      => 'ZECO:0000111',
    'temperature exposure'    => 'ZECO:0000160',
    );

my $db = Ace->connect(-path => $acedbpath, -program => $tace) or die('Connection failure: '. Ace->error);

my (@agm_annots, @gene_annots, @allele_annots);

my $it = $db->fetch_many(-query => 'Find Disease_model_annotation');

while( my $obj = $it->next) {
    if (not $obj->Date_last_updated or
	not $obj->Paper_evidence or
	not $obj->Evidence_code or
	not $obj->Disease_term) {
	warn("Bad object $obj - missing key fields - skipping\n" . $obj->asTable . "\n\n");
	next;
    }
    
    my  $evi_date = $obj->Date_last_updated;
    $evi_date->date_style('ace');
    my ($y, $m, $d) = split(/\-/, $evi_date);
    $evi_date = sprintf('%4d-%02d-%02dT00:00:00+00:00', $y, $m, $d);
    
    my ($paper) = &get_paper( $obj->Paper_evidence );
    my @evi_codes;
    for my $ec ($obj->Evidence_code) {
	if ($ec->right =~ /^ECO:/) {
	    push @evi_codes, $ec->right->name;
	}
	elsif (exists $go2eco{$ec->right}) {
	    push @evi_codes, $go2eco{$ec->right};
	}
	else {
	    warn "No ECO conversion available for evidence code $ec for $obj\n";
	}
    }
    unless (@evi_codes) {
	warn "No evidence codes could be converted to ECO for $obj - skipping\n";
	next;
    }
    
    # modifiers
    my $modifier_type;
    my (@genetic, @exp_cond);
    if ($obj->Modifier_association_type) {
	$modifier_type = $obj->Modifier_association_type->name;
	$modifier_type =~ s/^condition_//;
	$modifier_type = 'not_' . $modifier_type if $obj->Modifier_qualifier_not;
	
	my @mod_strain    = map { 'WB:' . $_->name } $obj->Modifier_strain;
	my @mod_transgene = map { 'WB:' . $_->name } $obj->Modifier_transgene;
	my @mod_var       = map { 'WB:' . $_->name } $obj->Modifier_variation;
	my @mod_gene      = map { 'WB:' . $_->name } $obj->Modifier_gene;
	my @mod_molecule  = map { $_->name } $obj->Modifier_molecule;
	my @mod_other     = map { $_->name } $obj->Other_modifier;
	
	@genetic  = (@mod_strain, @mod_transgene, @mod_var, @mod_gene);
	@exp_cond = (@mod_molecule, @mod_other);

	if (not @genetic and not @exp_cond) {
	    warn "$obj: Genetic or Experimental modifier info must be supplied when Modifier_association_type is supplied - skipping\n";
	    next;
	}
    }
    push @genetic, 'no_modifier' unless @genetic;
    
    for my $modifier (@genetic) {
	# [200507 mh6]
	# the crossReference should be annotation specific, but as the id changes every release the 
	# linking is not possible due to the lag

	unless (@evi_codes) {
	    push @evi_codes, $go2eco{'IMP'};
	}
	
	my $annot = {
	    mod_entity_id        => $obj->name,
	    internal             => JSON::false,
	    do_term_curie        => $obj->Disease_term->name,
	    data_provider        => 'WB',
	    date_updated         => $evi_date,
	    created_by_curie     => 'WB:curator',
	    annotation_type_name => 'manually_curated',
	    evidence_code_curies => \@evi_codes,
	    reference_curie      => $paper,
	    internal             => JSON::false,
	    obsolete             => JSON::false
	};
	$annot->{updated_by_curie} = $obj->Curator_confirmed ? 'WB:' . $obj->Curator_confirmed->name : "WB:curator";
	$annot->{genetic_sex_name} = $obj->Genetic_sex->name if $obj->Genetic_sex;
	$annot->{note_dtos} = [{
	    note_type_name => 'disease_summary',
	    internal       => JSON::false,
	    free_text      => $obj->Disease_model_description->name
	}] if $obj->Disease_model_description;
	
	unless ($modifier eq 'no_modifier') {
	    $annot->{disease_genetic_modifier_curie} = $modifier;
	    $annot->{disease_genetic_modifier_relation_name} = $modifier_type; # ameliorated_by / not_ameliorated_by / exacerbated_by / not_exacerbated_by
	}

	
	
	#  $annot->{disease_qualifiers} = ; # susceptibility / disease_progression / severity / onset / sexual_dimorphism / resistance / penetrance
	
	
	# [200310 mh6]
	# based on a list of annotations to skip for AGR from 
	# https://wiki.wormbase.org/index.php/Specifications_for_data_submission_to_the_Alliance
	# * had Qualifier_not as not to be submitted, but that changed with 3.1.0
	# Inducing_chemical, Modifier_molecule, and Other_molecule were also not submitted prior to 4.0
	
	next if ($obj->Interacting_variation
		 ||$obj->Interacting_transgene
		 ||$obj->Interacting_gene
		 ||$obj->RNAi_experiment
	    );
	
	my ($strain) = $obj->Strain;
	my ($allele) = $obj->Variation;
	my ($transgene) = $obj->Transgene;
	my ($gene) = $obj->Disease_relevant_gene;
	my ($genotype) = $obj->Genotype;
	my (@asserted_genes) = map { 'WB:'.$_->name } $obj->Asserted_gene;
	if (@asserted_genes) {
	    $annot->{asserted_gene_curies} = \@asserted_genes;
	}
	my ($obj_id, $obj_name, $obj_type);

	my $subject_field;
	if (defined $strain) {
	    $subject_field = 'agm_curie';
	    $obj_type = 'strain';
	    $obj_name = $strain->Public_name ? $strain->Public_name->name : $strain->name;
	    $obj_id = 'WB:' . $strain->name;
	} elsif (defined $allele) {
	    $subject_field = 'allele_curie';
	    if ($allele->Public_name){
		$obj_type = "allele";
		$obj_name = $allele->Public_name->name;
		$obj_id = 'WB:' . $allele->name;
	    }else{
		warn "$allele is missing a public name. Skipping ${\$obj->name}\n";
		next;
	    }
	} elsif (defined $transgene) {
	    $subject_field = 'allele_curie';
	    $obj_type = 'allele';# as quick fix for 3.0
	    $obj_name = $transgene->Public_name->name;
	    $obj_id = 'WB:' . $transgene->name;
	} elsif (defined $gene) {
	    $subject_field = 'gene_curie';
	    $obj_type = 'gene';
	    $obj_name = $gene->Public_name->name;
	    $obj_id = 'WB:' . $gene->name;
	    @inferred_genes = ();
	} elsif (defined $genotype){
	    $subject_field = 'agm_curie';
	    $obj_type = 'genotype';
	    $obj_name = "${\$genotype->Genotype_name}";
	    $obj_id = 'WB:' . $genotype->name;
	} else {
	    warn "Could not identify a central object for the annotation from Disease_model_annotation ${\$obj->name}\n";
	    next;
	}
	
	my $assoc_type = $obj->Association_type->name;
	$assoc_type = 'is_model_of' if $assoc_type !~ /model_of/ && ($obj_type eq 'strain' || $obj_type eq 'genotype');
	$assoc_type = 'is_implicated_in' if $obj_type eq 'allele';
	
	$annot->{disease_relation_name} = $assoc_type;
	$annot->{$subject_field} = $obj_id;
	$annot->{negated} = JSON::true if $obj->at('Qualifier_not');
	

	my $condition_relations = get_condition_relations($obj, $chebi_name_map);
	$annot->{condition_relation_dtos} = $condition_relations if @$condition_relations;
	
	if ($obj_type eq 'gene') {
	    push @gene_annots, $annot;
	}
	elsif ($obj_type eq 'allele') {
	    push @allele_annots, $annot;
	}
	else {
	    push @agm_annots, $annot;
	}
    }
    
}

$db->close;

my $all_annots = {linkml_version => $LINKML_SCHEMA};
$all_annots->{disease_allele_ingest_set} = \@allele_annots;
$all_annots->{disease_gene_ingest_set} = \@gene_annots;
$all_annots->{disease_agm_ingest_set} = \@agm_annots;

print_json($outfile, $all_annots);
    
exit(0);

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
    

sub get_paper {
  my ($wb_paper) = @_;

  my $pmid;
  foreach my $db ($wb_paper->Database) {
    if ($db->name eq 'MEDLINE') {
      $pmid = $db->right->right->name;
      last;
    }
  }
  my $publication_id = $pmid ? "PMID:$pmid" : "WB:$wb_paper";
  

  return $publication_id;
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
    my ($obj, $chebi_name_map) = @_;

    my $condition_relation_type;
    my (@modifiers, @inducers);
    my $conditions = [];
    my $condition_relations = [];
    if ($obj->Experimental_condition){
        $condition_relation_type = 'induced_by';
        my @inducing_chemicals = map {{
            condition_chemical_curie => get_chemical_ontology_id($_),
            condition_class_curie => $zeco{'chemical treatment'},
	    internal => JSON::false
            }} $obj->Inducing_chemical;
        my @inducing_agents     = map {{
            condition_class_curie => $zeco{'experimental conditions'},
	    condition_free_text => $_->name,
	    internal => JSON::false
            }} $obj->Inducing_agent;
        @inducers = (@inducing_chemicals, @inducing_agents);
        push @$conditions, @inducers;
    }
    if ($obj->Modifier_info and ($obj->Modifier_molecule or $obj->Other_modifier)) {
	$condition_relation_type = 'has_condition';
	if ($obj->Modifier_association_type) {
	    if ($obj->Modifier_association_type->name eq 'condition_ameliorated_by') {
		$condition_relation_type = 'ameliorated_by';
	    }
	    elsif ($obj->Modifier_association_type->name eq 'condition_exacerbated_by') {
		$condition_relation_type = 'exacerbated_by';
	    }
	}
        my @modifying_molecules = map {{
            condition_chemical_curie => get_chemical_ontology_id($_),
            condition_class_curie => $zeco{'chemical treatment'},
	    internal => JSON::false
            }} $obj->Modifier_molecule;
        my @other_modifiers = map {{
	    condition_class_curie => $zeco{'experimental conditions'},
	    condition_free_text_curie => $_->name,
	    internal => JSON::false
	    }} $obj->Other_modifier;
        @modifiers = (@modifying_molecules, @other_modifiers);
        push @$conditions, @modifiers;
    }
    
    unless ($condition_relation_type eq 'has_condition'){
	$condition_relation_type = 'not_' . $condition_relation_type if $obj->at('Modifier_qualifier_not');
    }

    if ($condition_relation_type) {
	my $cr_json = {
	    condition_relation_type_name => $condition_relation_type,
	    condition_dtos => $conditions,
	    internal => JSON::false
	};
        push @$condition_relations, $cr_json;
    }

    return ($condition_relations);
}


sub get_chebi_name_map {
    my %chebi_name_map;
    
    my $chebi_dom = XML::LibXML->load_xml(location => $CHEBI_PURL);
    for my $class ($chebi_dom->findnodes('//owl:Class')) {
	my $id = $class->findvalue('./oboInOwl:id');
	my $name = $class->findvalue('./rdfs:label');
	my $deprecated = $class->findvalue('./owl:deprecated');
	if ($deprecated ne 'true') {
	    if (defined $id and defined $name) {
		$chebi_name_map{$id} = $name;
	    }
	    else {
		warn ("Could not parse ChEBI name corresponding to $id\n");
	    }
	}
    }

    return \%chebi_name_map;
}


sub get_chemical_name {
    my ($obj, $chebi_name_map) = @_;

    my $id = get_chemical_ontology_id($_);
    if (exists $chebi_name_map->{$id}) {
	return $chebi_name_map->{$id};
    }
    return $obj->Public_name->name;
}

1;
