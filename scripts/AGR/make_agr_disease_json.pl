#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;

use lib $ENV{CVS_DIR};
use Wormbase;
use Modules::AGR;


my ($debug, $test, $verbose, $store, $wormbase, $build);
my ($outfile, $acedbpath, $ws_version, $outfh, $daf);

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
my $taxid = $wormbase->ncbi_tax_id;
my $full_name = $wormbase->full_name;

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

my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or die('Connection failure: '. Ace->error);

my ( $it, @annots);

#
# Old model
#
$it = $db->fetch_many(-query => 'find Gene Disease_info');

while (my $obj = $it->next) {
  next unless $obj->isObject();
  #next unless $obj->Species;
  #next unless $obj->Species->name eq $full_name;

  my $g = $obj->name;

  foreach my $doterm ($obj->Experimental_model) {
    my (@json_papers, $evi_date);
    foreach my $evi ($doterm->right->col) {
      if ($evi->name eq 'Paper_evidence') {
        foreach my $wb_paper ($evi->col ) {
          push @json_papers, &get_paper( $wb_paper );
        }
      } elsif ($evi->name eq 'Date_last_updated') {
        $evi_date = $evi->right;
        $evi_date->date_style('ace');
        my ($y, $m, $d) = split(/\-/, $evi_date);
        $evi_date = sprintf('%4d-%02d-%02dT00:00:00+00:00', $y, $m, $d);
      }
    }
    
    foreach my $pap (@json_papers) {
      push @annots, {
        objectId => "WB:$g",
        objectName => $obj->Public_name->name,
        DOid     => $doterm->name,
        dataProvider => [
          { 
            crossReference => {
              id => 'WB',
              pages => ['homepage'],
            },
            type => 'curated',
          },
        ],
        dateAssigned => defined $evi_date ? $evi_date : $date,
#       geneticSex   => 'hermaphrodite',
        evidence     => {
          evidenceCodes => [$go2eco{'IMP'}], # inferred from mutant phenotype; hard-coded for now
          publication => $pap,
        },
        objectRelation => {
          associationType => 'is_implicated_in',
          objectType      => 'gene',
        },
      };
    };
  }
}

#
# New model
#
$it = $db->fetch_many(-query => 'Find Disease_model_annotation');

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
  my @evi_codes = map { $go2eco{$_->name} } $obj->Evidence_code;

  my $annot = {
    DOid         => $obj->Disease_term->name,
    dataProvider => [
      { 
        crossReference => {
          id => 'WB',
          pages => ['homepage'],
        },
        type => 'curated',
      },
    ],
    dateAssigned => $evi_date,
#   geneticSex  => ($obj->Genetic_sex) ? $obj->Genetic_sex->name : 'hermaphrodite',
    evidence     => {
      evidenceCodes => (@evi_codes) ? \@evi_codes : [$go2eco{'IMP'}],
      publication => $paper,
    },
  };
  
  my ($strain) = $obj->Strain;
  my ($allele) = $obj->Variation;
  my ($transgene) = $obj->Transgene;
  my ($gene) = $obj->Disease_relevant_gene;
  my (@inferred_genes) = map { 'WB:' . $_->name } $obj->Inferred_gene;

  my ($obj_id, $obj_name, $obj_type, $assoc_type);

  my (@with_list) = map {'WB:'.$_->name} ($obj->Interacting_variation,$obj->Interacting_gene);

  if (defined $strain) {
    $obj_type = 'strain';
    $obj_name = ($strain->Genotype) ? $strain->Genotype->name : '';
    $assoc_type = 'is_model_of';
    $obj_id = 'WB:' . $strain->name;

    # WB/CalTech specific changes to the format
    push @with_list, "WB:" . $gene->name if (defined $gene && $build);
    push @with_list, "WBTransgene:" . $transgene->name if (defined $transgene && $build);
    push @with_list, "WBVar:" . $allele->name if (defined $allele && $build);
    
  } elsif (defined $allele) {
    $obj_type = "allele";
    $obj_name = $allele->Public_name->name;
    $assoc_type = 'is_implicated_in';
    $obj_id = 'WB:' . $allele->name;

    # WB/CalTech specific changes to the format
    push @with_list, "WB:" . $gene->name if (defined $gene && $build);
    push @with_list, "WB_Transgene:" . $transgene->name if (defined $transgene && $build);

  } elsif (defined $transgene) {
    $obj_type = 'transgene';
    $obj_name = $transgene->Public_name->name;
    $assoc_type = 'is_implicated_in';
    $obj_id = 'WB:' . $transgene->name;
    
    # WB/CalTech specific changes to the format
    push @with_list, "WB:" . $gene->name if (defined $gene && $build);

  } elsif (defined $gene) {
    $obj_type = 'gene';
    $obj_name = $gene->Public_name->name;
    $assoc_type = 'is_implicated_in';
    $obj_id = 'WB:' . $gene->name;

    @inferred_genes = ();
  } else {
    die "Could not identify a central object for the annotation from Disease_model_annotation $obj->name\n";
  }

  # 1. If an annotation has a gene and an allele— the primary annotation object is the allele
  # 2. If an annotation has a gene and a strain— the primary annotation object is the strain
  # 3. If an annotation has a gene and a transgene— the primary annotation object is the transgene

  # add gene / variation / allele as primaryGeneicEntityIDs
  my $primaryEntity;
  if (defined $allele){$primaryEntity = $allele->name}
  elsif (defined $strain){$primaryEntity = $strain->name}
  elsif (defined $transgene){$primaryEntity = $transgene->name}
  if (defined $primaryEntity && 'WB:'.$primaryEntity ne $obj_id){
	   $annot->{primaryGeneticEntityIDs} ||=[];
	   push @{$annot->{primaryGeneticEntityIDs}},'WB:'.$primaryEntity;
  }

  my $assoc_rel = {
    associationType => $assoc_type,
    objectType      => $obj_type,
  };
  $assoc_rel->{inferredGeneAssociation} = \@inferred_genes if @inferred_genes;

  $annot->{objectRelation} = $assoc_rel;
  $annot->{objectId} = $obj_id;
  $annot->{objectName} = $obj_name;
  $annot->{with} = \@with_list if @with_list;
  
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
    $annot->{qualifier} = 'not' if ($obj->at('Modifier_qualifier_not') && $build);
  }

  if ($obj->Experimental_condition){
    my @inducing_c     = map { $_->name } $obj->Inducing_chemical;
    my @inducing_a     = map { "$_" } $obj->Inducing_agent;
    my @exp_conditions = map {{textCondition => $_}} (@inducing_c,@inducing_a);

    # WB/CalTech specific changes
    $annot->{experimentalConditions} = \@exp_conditions if (@exp_conditions && $build);
  }

  push @annots, $annot;
}


$db->close;

if ($outfile) {
  open($outfh, ">$outfile") or die("cannot open $outfile : $!\n");  
} else {
  $outfh = \*STDOUT;
}

if ($daf) {
  &write_DAF_header( $outfh );
  &write_DAF_column_headers( $outfh );

  foreach my $annot (@annots) {
    # if an annotation has multiple lines of evidence, we need  a separate line for each one
    &write_DAF_line( $outfh, $annot );
  }
  
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

  my ($fh) = @_;

  print $fh "!daf-version 1.0\n";
  print $fh "!Date: $date\n";
  print $fh "!Project_name: WormBase (WB) Version $ws_version\n";
  print $fh "!URL: http://www.wormbase.org/\n";
  print $fh "!Contact Email: wormbase-help\@wormbase.org\n";
  print $fh "!Funding: NHGRI at US NIH, grant number U41 HG002223\n";

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
  my ($fh, $obj) = @_;

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
         ($obj->{qualifier}||''),
         (exists $obj->{modifier} and exists $obj->{modifier}->{genetic}) ? join(',', @{$obj->{modifier}->{genetic}}) : '',
         (exists $obj->{modifier} and exists $obj->{modifier}->{experimentalConditionsText}) ? join(',', @{$obj->{modifier}->{experimentalConditionsText}}) : '',
         join(",", @{$obj->{evidence}->{evidenceCodes}}),
         ($obj->{geneticSex}||''),
         $obj->{evidence}->{publication}->{publicationId},
         $date,
         'WB');
}
