#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;

use lib $ENV{CVS_DIR};
use Wormbase;


my ($debug, $test, $verbose, $store, $wormbase);
my ($outfile, $acedbpath, $ws_version, $outfh, $daf);

GetOptions (
  "debug=s"     => \$debug,
  "test"        => \$test,
  "verbose"     => \$verbose,
  "store:s"     => \$store,
  "database:s"  => \$acedbpath,
  "outfile:s"   => \$outfile,
  "wsversion=s" => \$ws_version,
  "writedaf"    => \$daf,
    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
      );
}

my $tace = $wormbase->tace;
my $date = &get_rfc_date();
my $alt_date = join("/", $date =~ /^(\d{4})(\d{2})(\d{2})/);
my $taxid = $wormbase->ncbi_tax_id;
my $full_name = $wormbase->full_name;

$acedbpath = $wormbase->autoace unless $acedbpath;
$ws_version = $wormbase->get_wormbase_version_name unless $ws_version;

if (not defined $outfile) {
  if ($daf) {
    $outfile = $wormbase->ontology."/disease_association.".$wormbase->get_wormbase_version_name.".daf.txt";
  } else {
    $outfile = "./wormbase.disease_association.${ws_version}.json";
  }
}

my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or die("Connection failure: ". Ace->error);

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
    my (@papers, $evi_date);
    foreach my $evi ($doterm->right->col) {
      if ($evi->name eq 'Paper_evidence') {
        @papers = &get_papers( $evi->col );
      } elsif ($evi->name eq 'Date_last_updated') {
        $evi_date = $evi->right;
        $evi_date->date_style('ace');
        my ($y, $m, $d) = split(/\-/, $evi_date);
        $evi_date = sprintf("%4d-%02d-%02dT00:00:00+00:00", $y, $m, $d);
      }
    }
    
    push @annots, {
      objectId => "WB:$g",
      DOid     => $doterm->name,
      taxonId  => "NCBITaxon:$taxid",
      dataProvider => "WB", 
      dateAssigned => defined $evi_date ? $evi_date : $date,
      geneticSex   => "hermaphrodite",
      evidence     => [
        {
          evidenceCode => "IMP",  # inferred from mutant phenotype; hard-coded for now
          publications => \@papers,
        },
          ],
      objectRelation => {
        associationType => "causes_condition",
        objectType      => "gene",
      },
    };
  }
}

#
# New model
#
$it = $db->fetch_many(-query => 'Find Disease_model_annotation');

while( my $obj = $it->next) {

  my  $evi_date = $obj->Date_last_updated;
  $evi_date->date_style('ace');
  my ($y, $m, $d) = split(/\-/, $evi_date);
  $evi_date = sprintf("%4d-%02d-%02dT00:00:00+00:00", $y, $m, $d);

  my @papers = &get_papers( $obj->Paper_evidence );

  my $annot = {
    DOid         => $obj->Disease_term->name,
    taxonId      => "NCBITaxon:$taxid",
    dataProvider => "WB",
    dateAssigned => $evi_date,
    geneticSex  => ($obj->Genetic_sex) ? $obj->Genetic_sex->name : "hermaphrodite",
    evidence     => [
      {
        evidenceCode => ($obj->Evidence_code) ? $obj->Evidence_code->name : "XXX",
        publications => \@papers,
      },
      ],
  };
  
  # no Experimental_condition data in WB; will need to address this one as data arrives
  
  my ($strain) = $obj->Strain;
  my ($allele) = $obj->Variation;
  my ($transgene) = $obj->Transgene;
  my ($gene) = $obj->Disease_relevant_gene;
  my ($inferred_gene) = $obj->Inferred_gene;

  my ($obj_id, $obj_type, $assoc_type, @with_list);
  if (defined $strain) {
    $obj_type = "strain";
    $assoc_type = "is_model_of";
    $obj_id = "WB_Strain:" . $strain->name;

    #push @with_list, "WB:" . $gene->name if defined $gene;
    #push @with_list, "WB_Transgene:" . $transgene->name if defined $transgene;
    #push @with_list, "WB_Var:" . $allele->name if defined $name;
  } elsif (defined $allele) {
    $obj_type = "allele";
    $assoc_type = "causes_condition";
    $obj_id = "WB_Var:" . $allele->name;

    #push @with_list, "WB:" . $gene->name if defined $gene;
    #push @with_list, "WB_Transgene:" . $transgene->name if defined $transgene;
  } elsif (defined $transgene) {
    $obj_type = "transgene";
    $assoc_type = "causes_condition";
    $obj_id = "WB_Transgene:" . $transgene->name;
    
    #push @with_list, "WB:" . $gene->name if defined $gene;
  } elsif (defined $gene) {
    $obj_type = "gene";
    $assoc_type = "causes_condition";
    $obj_id = "WB:" . $gene->name;
  } else {
    die "Could not identify a central object for the annotation from Disease_model_annotation $obj->name\n";
  }

  my $assoc_rel = {
    associationType => $assoc_type,
    objectType      => $obj_type,
  };
  $assoc_rel->{inferredGeneAssociation} = $inferred_gene->name if defined $inferred_gene;

  $annot->{objectRelation} = $assoc_rel;
  $annot->{objectId} = $obj_id;
  $annot->{with} = \@with_list if @with_list;
  
  #
  # lastly, modifiers
  #
  if ($obj->Modifier_association_type) {
    my ($mod_assoc_type) = $obj->Modifier_association_type->name;

    my @mod_strain    = map { "WB_Strain:" . $_->name } $obj->Modifier_strain;
    my @mod_transgene = map { "WB_Transgene:" . $_->name } $obj->Modifier_transgene;
    my @mod_gene      = map { "WB:" . $_->name } $obj->Modifier_gene;
    my @mod_molecule  = map { $_->name } $obj->Modifier_molecule;
    my @mod_other     = map { $_->name } $obj->Other_modifier;

    my $mod_annot = {
      associationType => $mod_assoc_type,
    };
    my @genetic  = (@mod_strain, @mod_transgene, @mod_gene);
    my @exp_cond = (@mod_molecule, @mod_other);

    die "$obj: Genetic or Experimental modifier info must be supplied when Modifier_association_type is supplied\n"
        if not @genetic and not @exp_cond;

    $mod_annot->{genetic} = \@genetic if @genetic;
    $mod_annot->{experimentalConditionsText} = \@exp_cond if @exp_cond;

    $annot->{modifier} = $mod_annot;
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
    foreach my $evidence (@{$annot->{evidence}}) {
      my $evi_code = $evidence->{evidenceCode};
      foreach my $paper (@{$evidence->{publications}}) {
        my $pmid = $paper->{pubMedId};

        &write_DAF_line( $outfh, $evi_code, $pmid, $annot );
      }
    }
  }
  
} else {
  my $meta_data = {
    dateProduced => $date,
    dataProvider => "WB", 
    release      => (defined $ws_version) ? $ws_version : $wormbase->get_wormbase_version_name(),
  };
  
  my $data = {
    metaData => $meta_data,
    data => \@annots,
  };
  
  my $json_obj = JSON->new;
  my $string = $json_obj->allow_nonref->canonical->pretty->encode($data);
  print $outfh $string;
}

exit(0);

##############################################
sub get_papers {
  my @wb_paper = @_;

  my @papers;

  foreach my $paper (@wb_paper) {
    my $pmid;
    foreach my $db ($paper->Database) {
      if ($db->name eq 'MEDLINE') {
        $pmid = $db->right->right->name;
        last;
      }
    }
    push @papers, {
      modPublicationId => "WB_REF:$paper",
    };
    if ($pmid) {
      $papers[-1]->{pubMedId} = "PMID:$pmid";
    }
  }

  return @papers;
}

##############################################
sub get_rfc_date {

  my $date;
  
  open(my $date_fh, "date --rfc-3339=seconds |");
  while(<$date_fh>) {
    if (/^(\S+)\s+(\S+)/) {
      $date = "${1}T${2}";
    }
  }
  
  return $date;
}

##############################################
sub write_DAF_header {

  my ($fh) = @_;

  print $fh "!daf-version 0.1\n";
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
  my ($fh, $evi_code, $pmid,  $obj) = @_;

  my $obj_id = $obj->{objectId};
  my ($id_source, $id) = split(/:/, $obj_id);

  my $date = $obj->{dateAssigned};
  $date =~ s/(\d{4})\-(\d{2})\-(\d{2}).+/${1}${2}${3}/; 
  
  printf($fh "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n", 
         $obj->{taxonId},
         $obj->{objectRelation}->{objectType}, 
         $id_source,
         $id,
         "",
         (exists $obj->{inferredGeneAssociation}) ? $obj->{inferredGeneAssociation} : "",
         "",
         "",
         $obj->{objectRelation}->{associationType},
         "",
         $obj->{DOid},
         (exists $obj->{with}) ? $obj->{with} : "",
         (exists $obj->{modifier}) ? $obj->{modifier}->{associationType} : "",
         "",
         (exists $obj->{modifier} and exists $obj->{modifier}->{genetic}) ? join(",", @{$obj->{modifier}->{genetic}}) : "",
         (exists $obj->{modifier} and exists $obj->{modifier}->{experimentalConditionsText}) ? join(",", @{$obj->{modifier}->{experimentalConditionsText}}) : "",
         $evi_code,
         $obj->{geneticSex},
         $pmid,
         $date,
         $obj->{dataProvider});

}
