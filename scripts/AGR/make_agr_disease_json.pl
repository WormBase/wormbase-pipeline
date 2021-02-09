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
use GAF;

my ($debug, $test, $verbose, $store, $wormbase, $build);
my ($outfile, $acedbpath, $ws_version, $outfh, $daf,$white);

GetOptions (
  'debug=s'     => \$debug,
  'test'        => \$test,
  'verbose'     => \$verbose,
  'store:s'     => \$store,
  'database:s'  => \$acedbpath,
  'outfile:s'   => \$outfile,
  'wsversion=s' => \$ws_version,
  'writedaf'    => \$daf,
  'build'       => \$build, # to enable WB/CalTech build specific changes
  'AGRwhitelist'=> \$white,
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

if (not defined $outfile) {
  if ($daf) {
    $outfile = $wormbase->ontology.'/disease_association.'.$wormbase->get_wormbase_version_name.'.daf.txt';
  } else {
    $outfile = "./wormbase.disease_association.${ws_version}.json";
  }
}

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

my ( $it, @annots);
my $wl = get_whitelist() if $white;

# New model
#
$it = $db->fetch_many(-query => 'Find Disease_model_annotation');

my %taxon_ids;
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
  my @evi_codes = map { $go2eco{$_->right} } $obj->Evidence_code;

  # [200507 mh6]
  # the crossReference should be annotation specific, but as the id changes every release the 
  # linking is not possible due to the lag

  my $annot = {
    DOid         => $obj->Disease_term->name,
    dataProvider => [
      { 
        crossReference => {
          id => $obj->Disease_term->name,
          pages => ['disease/wb'],
        },
        type => 'curated',
      },
    ],
    dateAssigned => $evi_date,
    evidence     => {
      evidenceCodes => (@evi_codes) ? \@evi_codes : [$go2eco{'IMP'}],
      publication => $paper,
    },
  };
  

  # [200310 mh6]
  # based on a list of annotations to skip for AGR from 
  # https://wiki.wormbase.org/index.php/Specifications_for_data_submission_to_the_Alliance
  # * had Qualifier_not as not to be submitted, but that changed with 3.1.0
  # Inducing_chemical, Modifier_molecule, and Other_molecule were also not submitted prior to 4.0
  unless ($build) {
	  next if ($obj->Interacting_variation
	         ||$obj->Interacting_transgene
	         ||$obj->Interacting_gene
	         ||$obj->RNAi_experiment
	         ||$obj->Modifier_transgene
	         ||$obj->Modifier_variation
	         ||$obj->Modifier_gene
	  );
  }

  my ($strain) = $obj->Strain;
  my ($allele) = $obj->Variation;
  my ($transgene) = $obj->Transgene;
  my ($gene) = $obj->Disease_relevant_gene;
  my ($genotype) = $obj->Genotype;
  my (@inferred_genes) = map { 'WB:'.$_->name } $obj->Inferred_gene;
  my ($obj_id, $obj_name, $obj_type);
  my (@with_list) = map {'WB:'.$_->name} ($obj->Interacting_variation,$obj->Interacting_gene,$obj->Interacting_transgene);
  my $assoc_type = 'is_implicated_in';
  
  if ($white && $allele){ # remove annotations unless they are in the AGR allele variations
	  next unless $wl->{"$allele"}
  }

  if (defined $strain) {
    $obj_type = 'strain';
    $obj_name = $strain->Public_name ? $strain->Public_name->name : $strain->name;
    $assoc_type = 'is_model_of';
    $obj_id = 'WB:' . $strain->name;
    if ($strain->Species) {
	$taxon_ids{$obj_id} = $strain->Species->NCBITaxonomyID;
    }
    else {
	warn "Could not identify species for $obj_id, defaulting to c. elegans\n";
	$taxon_ids{$obj_id} = 6239;
    }
    # WB/CalTech specific changes to the format
    push @with_list, "WB:" . $gene->name if (defined $gene && $build);
    push @with_list, "WBTransgene:" . $transgene->name if (defined $transgene && $build);
    push @with_list, "WBVar:" . $allele->name if (defined $allele && $build);
    
  } elsif (defined $allele) {
    if ($allele->Public_name){
	$obj_type = "allele";
	$obj_name = $allele->Public_name->name;
	$obj_id = 'WB:' . $allele->name;
	if ($allele->Species) {
	    $taxon_ids{$obj_id} = $allele->Species->NCBITaxonomyID;
	}
	else {
	    warn "Could not identify species for $obj_id, defaulting to c. elegans\n";
	    $taxon_ids{$obj_id} = 6239;
	}
    }else{
	    warn "$allele is missing a public name. Skipping ${\$obj->name}\n";
	    next;
    }

    # WB/CalTech specific changes to the format
    push @with_list, "WB:" . $gene->name if (defined $gene && $build);
    push @with_list, "WBTransgene:" . $transgene->name if (defined $transgene && $build);

  } elsif (defined $transgene) {
    $obj_type = $build ? 'transgene' : 'allele';# as quick fix for 3.0
    $obj_name = $transgene->Public_name->name;
    $obj_id = 'WB:' . $transgene->name;
    if ($transgene->Species) {
	$taxon_ids{$obj_id} = $transgene->Species->NCBITaxonomyID;
    }
    else {
        warn "Could not identify species for $obj_id, defaulting to c. elegans\n";
	$taxon_ids{$obj_id} = 6239;
    }   
    # WB/CalTech specific changes to the format
    push @with_list, "WB:" . $gene->name if (defined $gene && $build);

  } elsif (defined $gene) {
    $obj_type = 'gene';
    $obj_name = $gene->Public_name->name;
    $obj_id = 'WB:' . $gene->name;
    if ($gene->Species) {
	$taxon_ids{$obj_id} = $gene->Species->NCBITaxonomyID;
    }
    else {
        warn "Could not identify species for $obj_id, defaulting to c. elegans\n";
	$taxon_ids{$obj_id} = 6239;
    }

    @inferred_genes = ();
  } elsif (defined $genotype){
      $obj_type = 'genotype';
      $obj_name = "${\$genotype->Genotype_name}";
      $obj_id = 'WB:' . $genotype->name;
      if ($genotype->Species) {
	  $taxon_ids{$obj_id} = $genotype->Species->NCBITaxonomyID;
      }
      else {
	  warn "Could not identify species for $obj_id, defaulting to c. elegans\n";
	  $taxon_ids{$obj_id} = 6239;
      }
  } else {
    warn "Could not identify a central object for the annotation from Disease_model_annotation ${\$obj->name}\n";
    next;
  }

  $assoc_type = $obj->Association_type->name if $obj->Association_type and !defined $strain and !defined $allele;
  
  my $assoc_rel = {
    associationType => $assoc_type,
    objectType      => $obj_type,
  };
  $assoc_rel->{inferredGeneAssociation} = \@inferred_genes if @inferred_genes && ($daf || $build); # hack for Ranjana

  $annot->{objectRelation} = $assoc_rel;
  $annot->{objectId} = $obj_id;
  $annot->{objectName} = $obj_name;
  $annot->{with} = \@with_list if (@with_list && $build);
  $annot->{negation} = 'not' if $obj->at('Qualifier_not');

  # modifiers
  
  if ($obj->Modifier_association_type) {
    my ($mod_assoc_type) = $obj->Modifier_association_type->name;

    my @mod_strain    = map { 'WB:' . $_->name } $obj->Modifier_strain;
    my @mod_transgene = map { 'WB:' . $_->name } $obj->Modifier_transgene;
    my @mod_var       = map { 'WB:' . $_->name } $obj->Modifier_variation;
    my @mod_gene      = map { 'WB:' . $_->name } $obj->Modifier_gene;
    my @mod_molecule  = map { $_->name } $obj->Modifier_molecule;
    my @mod_other     = map { $_->name } $obj->Other_modifier;
 
    my $mod_annot = {
      associationType => $mod_assoc_type,
    };
    my @genetic  = (@mod_strain, @mod_transgene, @mod_var, @mod_gene);
    my @exp_cond = (@mod_molecule, @mod_other);

    die "$obj: Genetic or Experimental modifier info must be supplied when Modifier_association_type is supplied\n"
        if not @genetic and not @exp_cond;

    $mod_annot->{genetic} = \@genetic if @genetic;
    $mod_annot->{experimentalConditionsText} = \@exp_cond if @exp_cond;

    # WB/CalTech specific changes
    $annot->{modifier} = $mod_annot if $build;
    # $annot->{negation} = 'not' if $obj->at('Modifier_qualifier_not'); # according to Ranjana, it should not be included
  }

  if ($build && $obj->Experimental_condition){
    my @inducing_c     = map { $_->name } $obj->Inducing_chemical;
    my @inducing_a     = map { "$_" } $obj->Inducing_agent;
    my @exp_conditions = map {{textCondition => $_}} (@inducing_c,@inducing_a);

    # WB/CalTech specific changes
    $annot->{experimentalConditions} = \@exp_conditions if @exp_conditions;
  }

  my $conditions = get_condition_relations($obj);
  $annot->{conditionRelations} = $conditions if @$conditions;
  
  push @annots, $annot;# unless ($obj_type eq 'transgene' && ! $build);

  # here needs to be a bit of logic that adds the secondary annotations
  foreach my $secondary ($strain,$allele,$transgene,$gene){
	  next unless $secondary;

	  my $class = lc $secondary->class;
	  $class = 'allele' if $class eq 'variation'; # AGR terminology hack

	  my $secondaryPrefix = ($class eq 'transgene' && $build) ? 'WBTransgene:' : 'WB:';
	  my $sname = $secondaryPrefix.$secondary->name;

	  next if $obj_id eq $sname; # skip the primary annotation
	  $taxon_ids{$sname} = $taxon_ids{$obj_id};

	  my $secondaryAnnotation = dclone($annot);
 
          $secondaryAnnotation->{objectRelation}->{objectType} = $class;
          $secondaryAnnotation->{objectId} = $sname;

          my $public_name = $secondary->Public_name ? $secondary->Public_name->name : $secondary->name;
	  $secondaryAnnotation->{objectName}=$public_name;
          
          my $secondaryAssociationType = 'is_implicated_in';
	  $secondaryAssociationType = 'is_model_of' if $class eq 'strain';
	  $secondaryAnnotation->{objectRelation}->{associationType}=$secondaryAssociationType;

	  # add gene / variation / allele as primaryGeneicEntityIDs
          my $primaryEntity;
          if (defined $strain){$primaryEntity = 'WB:'.$strain->name}
          elsif (defined $allele){$primaryEntity = 'WB:'.$allele->name}
          elsif (defined $transgene){
		  $primaryEntity = ($build ? 'WBTransgene:' : 'WB:').$transgene->name
	  }
          if (defined $primaryEntity && $primaryEntity ne $sname){
	   $secondaryAnnotation->{primaryGeneticEntityIDs} ||=[];
	   push @{$secondaryAnnotation->{primaryGeneticEntityIDs}},$primaryEntity;
          }

          push @annots, $secondaryAnnotation; # unless ($class eq 'transgene' && ! $build);
  }

}

$db->close;

if ($outfile) {
  open($outfh, ">$outfile") or die("cannot open $outfile : $!\n");  
} else {
  $outfh = \*STDOUT;
}

if ($daf) {
  &write_DAF_header( $outfh, $wormbase->get_wormbase_version_name );
  &write_DAF_column_headers( $outfh );

  foreach my $annot (@annots) {
    # if an annotation has multiple lines of evidence, we need  a separate line for each one
    &write_DAF_line( $outfh, $annot, $taxon_ids{$annot->{objectId}} );
  }

  &make_species_files($wormbase, $outfile, 1);
  
} else {

  my $data = {
    metaData => AGR::get_file_metadata_json( (defined $ws_version) ? $ws_version : $wormbase->get_wormbase_version_name() ), 
    data => \@annots,
  };
  
  my $json_obj = JSON->new;
  my $string = $json_obj->allow_nonref->canonical->pretty->encode($data);
  print $outfh $string;
}

exit(0);

##############################################
sub get_paper {
  my ($wb_paper) = @_;

  my $json_paper = {};

  my $pmid;
  foreach my $db ($wb_paper->Database) {
    if ($db->name eq 'MEDLINE') {
      $pmid = $db->right->right->name;
      last;
    }
  }
  $json_paper->{publicationId} = $pmid ? "PMID:$pmid" : "WB:$wb_paper";
  $json_paper->{crossReference} = {id =>"WB:$wb_paper",pages => ['reference']};

  return $json_paper;
}

##############################################
sub write_DAF_header {

  my ($fh, $release) = @_;

  my $header_date = DateTime->now;
  print $fh "\!daf-version 1.0\n";
  print $fh "\!generated-by: WormBase\n";
  print $fh "\!date-generated: " . $header_date->ymd . "\n";
  print $fh "\!project-URL: https://wormbase.org\n";
  print $fh "\!specification-URL: https://wiki.wormbase.org/index.php/WormBase_gene_association_file\n";
  print $fh "\!project-release: $release\n";
  print $fh "\!Contact Email: help\@wormbase.org\n";
  print $fh "\!Funding: NHGRI at US NIH, grant number U41 HG002223\n";

}

###########################################################3
sub write_DAF_column_headers {
  my ($fh) = @_;

  my @col_headers  = (
    'Taxon',
    'DB Object Type',
    'DB',
    'DB Object ID',
    'DB Object Symbol',
    'Inferred gene association',
    'Gene Product Form ID',
    # 'Additional genetic components', ???
    'Experimental conditions',
    'Association type',
    'Qualifier',
    'DO ID',
    'With',
    'Modifier - assocation type',
    'Modifier - Qualifier',
    'Modifier - genetic',
    'Modifier - experimental conditions',
    'Evidence Code',
    'genetic sex',
    'DB:Reference',
    'Date',
    'Assigned By');

  print $fh join("\t", @col_headers), "\n";
}


###########################################################3
sub write_DAF_line {
  my ($fh, $obj, $taxid) = @_;

  my $date = $obj->{dateAssigned};
  $date =~ s/(\d{4})\-(\d{2})\-(\d{2}).+/${1}${2}${3}/; 

  my %inferred_genes;
  if ($obj->{objectRelation}->{objectType} eq 'gene') {
    $inferred_genes{$obj->{objectId}}=1;
  }elsif (exists $obj->{objectRelation}->{inferredGeneAssociation}) {
    map {$inferred_genes{$_} =1} @{$obj->{objectRelation}->{inferredGeneAssociation}};
  }

  # another crude -build block
  if ($build && exists $obj->{objectRelation}->{inferredGeneAssociation}){
    map {$inferred_genes{$_} =1} @{$obj->{objectRelation}->{inferredGeneAssociation}};
  }

  my $inferred_gene =  join(",", keys %inferred_genes);

  printf($fh "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
         $taxid,
         $obj->{objectRelation}->{objectType}, 
         split(/:/, $obj->{objectId}), 
         $obj->{objectName},
         $inferred_gene,
         '',
         (exists $obj->{experimentalConditions}) ? join(',',map{'"'.$_->{textCondition}.'"'} @{$obj->{experimentalConditions}}):'',
         $obj->{objectRelation}->{associationType},
         '',
         $obj->{DOid},
         (exists $obj->{with}) ? join(',',@{$obj->{with}}) : '',
         (exists $obj->{modifier}) ? $obj->{modifier}->{associationType} : '',
         ($obj->{negation}||''),
         (exists $obj->{modifier} and exists $obj->{modifier}->{genetic}) ? join(',', @{$obj->{modifier}->{genetic}}) : '',
         (exists $obj->{modifier} and exists $obj->{modifier}->{experimentalConditionsText}) ? join(',', @{$obj->{modifier}->{experimentalConditionsText}}) : '',
         join(",", @{$obj->{evidence}->{evidenceCodes}}),
         ($obj->{geneticSex}||''),
         $obj->{evidence}->{publication}->{publicationId},
         $date,
         'WB');
}

sub get_whitelist{
	my $it = $db->fetch_many(-query =>'find Variation WHERE Live AND COUNT(Gene) == 1 AND (Phenotype OR Disease_info OR Interactor) AND NOT Natural_variant');
	my %whitelist;
	while (my $v = $it->next){
		$whitelist{"$v"}=1;
	}
	return \%whitelist;
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
    my $obj = shift;

    my $condition_relation_type = 'has_condition';
    if ($obj->Modifier_association_type) {
        if ($obj->Modifier_association_type->name eq 'condition_ameliorated_by') {
            $condition_relation_type = 'ameliorates';
        }
        elsif ($obj->Modifier_association_type->name eq 'condition_exacerbated_by') {
            $condition_relation_type = 'exacerbates';
        }
    }

    my (@modifiers, @inducers, @conditions);
    if ($obj->Experimental_condition){
        $condition_relation_type = 'induces';
        my @inducing_chemicals = map {{
            conditionStatement => 'chemical treatment:' . $_->Public_name->name,
            chemicalOntologyId => get_chemical_ontology_id($_),
            conditionClassId => $zeco{'chemical treatment'}
            }} $obj->Inducing_chemical;
        my @inducing_agents     = map {{
            conditionStatement => 'experimental conditions:' . $_->name,
            conditionClassId => $zeco{'experimental conditions'}
            }} $obj->Inducing_agent;
        @inducers = (@inducing_chemicals, @inducing_agents);
        push @conditions, {conditionRelationType => $condition_relation_type,
                           conditions            => \@inducers};
    }
    if ($obj->Modifier_info) {
        if ($obj->Modifier_association_type) {
            if ($obj->Modifier_association_type->name eq 'condition_ameliorated_by') {
                $condition_relation_type = 'ameliorates';
            }
            elsif ($obj->Modifier_association_type->name eq 'condition_exacerbated_by') {
                $condition_relation_type = 'exacerbates';
            }
            else{
                $condition_relation_type = 'has_condition';
            }
        }
        my @modifying_molecules = map {{
            conditionStatement => 'chemical treatment:' . $_->Public_name->name,
            chemicalOntologyId => get_chemical_ontology_id($_),
            conditionClassId => $zeco{'chemical treatment'}
            }} $obj->Modifier_molecule;
        my @other_modifiers = map {{
          conditionStatement => 'experimental conditions:' . $_->name,
          conditionClassId => $zeco{'experimental conditions'}
        }} $obj->Other_modifier;
        @modifiers = (@modifying_molecules, @other_modifiers);
        push @conditions, {conditionRelationType => $condition_relation_type,
                           conditions            => \@modifiers};
    }

    return \@conditions;
}
