#!/usr/bin/env perl
# Note! This script is used both by the WormBase build (autoace.pl -ontologies), and the AGR data dump, so if you make changes test on both workflows.

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
use GAF;

my ($debug, $test, $verbose, $store, $wormbase, $acedbpath, $ws_version, $agm_out, $allele_out, $gene_out, $cr_out, $exp_out);

GetOptions (
    'debug=s'      => \$debug,
    'test'         => \$test,
    'verbose'      => \$verbose,
    'store:s'      => \$store,
    'database:s'   => \$acedbpath,
    'agm_out:s'    => \$agm_out,
    'gene_out:s'   => \$gene_out,
    'allele_out:s' => \$allele_out,
    'cr_out:s'     => \$cr_out,
    'exp_out:s'    => \$exp_out,
    'wsversion=s'  => \$ws_version,
)||die("unknown command line option: $@\n");

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

$acedbpath = $wormbase->autoace unless $acedbpath;
$ws_version = $wormbase->get_wormbase_version_name unless $ws_version;

$agm_out = "./wormbase.agm_disease_association.${ws_version}.json" unless defined $agm_out;
$gene_out = "./wormbase.gene_disease_association.${ws_version}.json" unless defined $gene_out;
$allele_out = "./wormbase.allele_disease_association.${ws_version}.json" unless defined $allele_out;
$cr_out = "./wormbase.condition_relations.${ws_version}.json" unless defined $cr_out;
$exp_out = "./wormbase.experimental_conditions.${ws_version}.json" unless defined $exp_out;


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
my $condition_relations = {};
my $experimental_conditions = {};

my $it = $db->fetch_many(-query => 'Find Disease_model_annotation');

while( my $obj = $it->next) {
    if (not $obj->Date_last_updated or
	not $obj->Paper_evidence or
	not $obj->Evidence_code or
	not $obj->Disease_term) {
	warn("Bad object $obj - missing key fields - skipping\n");
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
	
	my $annot = {
	    object         => $obj->Disease_term->name,
	    data_provider => ['WB'],
	    date_last_modified => $evi_date,
	    evidence_codes => (@evi_codes) ? \@evi_codes : [$go2eco{'IMP'}],
	    annotation_reference => $paper,
	};
	$annot->{modified_by} = 'WB:' . $obj->Curator_confirmed->name if $obj->Curator_confirmed;
	$annot->{genetic_sex} = $obj->Genetic_sex->name if $obj->Genetic_sex;
	$annot->{disease_annotation_summary} = $obj->Disease_model_description->name if $obj->Disease_model_description; # or could be disease_annotation_note
	
	unless ($modifier eq 'no_modifier') {
	    $annot->{disease_genetic_modifier} = $modifier;
	    $annot->{disease_genetic_modifier_relation} = $modifier_type; # ameliorated_by / not_ameliorated_by / exacerbated_by / not_exacerbated_by
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
	my (@inferred_genes) = map { 'WB:'.$_->name } $obj->Inferred_gene;
	my ($obj_id, $obj_name, $obj_type);
	my $assoc_type = 'is_implicated_in';
	
	if (defined $strain) {
	    $obj_type = 'strain';
	    $obj_name = $strain->Public_name ? $strain->Public_name->name : $strain->name;
	    $assoc_type = 'is_model_of';
	    $obj_id = 'WB:' . $strain->name;
	} elsif (defined $allele) {
	    if ($allele->Public_name){
		$obj_type = "allele";
		$obj_name = $allele->Public_name->name;
		$obj_id = 'WB:' . $allele->name;
	    }else{
		warn "$allele is missing a public name. Skipping ${\$obj->name}\n";
		next;
	    }
	} elsif (defined $transgene) {
	    $obj_type = 'allele';# as quick fix for 3.0
	    $obj_name = $transgene->Public_name->name;
	    $obj_id = 'WB:' . $transgene->name;
	} elsif (defined $gene) {
	    $obj_type = 'gene';
	    $obj_name = $gene->Public_name->name;
	    $obj_id = 'WB:' . $gene->name;
	    @inferred_genes = ();
	} elsif (defined $genotype){
	    $obj_type = 'genotype';
	    $obj_name = "${\$genotype->Genotype_name}";
	    $obj_id = 'WB:' . $genotype->name;
	    $assoc_type = 'is_model_of';
	} else {
	    warn "Could not identify a central object for the annotation from Disease_model_annotation ${\$obj->name}\n";
	    next;
	}
	
	$assoc_type = $obj->Association_type->name if $obj->Association_type and !defined $strain and !defined $allele and !defined $genotype;
	
	$annot->{predicate} = $assoc_type;
	$annot->{subject} = $obj_id;
	$annot->{negated} = 'true' if $obj->at('Qualifier_not');
	

	my $cr_curie;
	($cr_curie, $condition_relations, $experimental_conditions) = get_condition_relations($obj, $condition_relations, $experimental_conditions);
	$annot->{conditionRelations} = [$cr_curie] if $cr_curie;

	
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

print_json($agm_out, \@agm_annots);
print_json($gene_out, \@gene_annots);
print_json($allele_out, \@allele_annots);
my @cr_annots = values %{$condition_relations};
my @exp_annots = values %{$experimental_conditions};
print_json($cr_out, \@cr_annots);
print_json($exp_out, \@exp_annots);
    
exit(0);

##############################################

sub print_json {
    my ($outfile, $data) = @_;

    open(my $out_fh, ">$outfile") or die("cannot open $outfile : $!\n");  
    my $json_obj = JSON->new;
    my $string = $json_obj->allow_nonref->canonical->pretty->encode({data => $data});
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
    my ($obj, $condition_relations, $experimental_conditions) = @_;

    my $condition_relation_type; 
    
    my (@modifiers, @inducers);
    my $conditions = [];
    if ($obj->Experimental_condition){
        $condition_relation_type = 'induced_by';
        my @inducing_chemicals = map {{
            condition_statement => 'chemical treatment:' . $_->Public_name->name,
            condition_chemical => get_chemical_ontology_id($_),
            condition_class => $zeco{'chemical treatment'},
            }} $obj->Inducing_chemical;
        my @inducing_agents     = map {{
            condition_statement => 'experimental conditions:' . $_->name,
            condition_class => $zeco{'experimental conditions'}
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
            condition_statement => 'chemical treatment:' . $_->Public_name->name,
            condition_chemical => get_chemical_ontology_id($_),
            condition_class => $zeco{'chemical treatment'}
            }} $obj->Modifier_molecule;
        my @other_modifiers = map {{
	    condition_statement => 'experimental conditions:' . $_->name,
	    condition_class => $zeco{'experimental conditions'}
	    }} $obj->Other_modifier;
        @modifiers = (@modifying_molecules, @other_modifiers);
        push @$conditions, @modifiers;
    }
    
    my ($conditions, $exp_cond_curies) = generate_experimental_condition_curies($conditions);

    unless ($condition_relation_type eq 'has_condition'){
	$condition_relation_type = 'not_' . $condition_relation_type if $obj->at('Modifier_qualifier_not');
    }

    my $cr_curie;
    if ($condition_relation_type) {
	$cr_curie = join('|', $condition_relation_type, @{$exp_cond_curies});
	my $cr_json = {
	    condition_relation_type => $condition_relation_type,
	    conditions => $exp_cond_curies,
	    curie => $cr_curie
	};
	$condition_relations->{$cr_curie} = $cr_json;

    }

    for my $condition (@$conditions) {
	$experimental_conditions->{$condition->{curie}} = $condition;
    }
    
    return ($cr_curie, $condition_relations, $experimental_conditions);
}


sub generate_experimental_condition_curies {
    my ($conditions) = @_;

    my (@conditions_with_curies, @exp_cond_curies);
    for my $condition (@$conditions) {
	my @curie_components;
	push @curie_components, $condition->{condition_chemical} if exists $condition->{condition_chemical};
	push @curie_components, $condition->{condition_class} if exists $condition->{condition_class};
	push @curie_components, $condition->{condition_statement} if exists $condition->{condition_statement};
	my $curie = join('|', @curie_components);
	$condition->{curie} = $curie;
	push @conditions_with_curies, $condition;
	push @exp_cond_curies, $curie;
    }
    @exp_cond_curies = sort @exp_cond_curies;
    
    return (\@conditions_with_curies, \@exp_cond_curies);
}
