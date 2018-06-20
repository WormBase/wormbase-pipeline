#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;

use lib $ENV{CVS_DIR};
use Wormbase;


my $MMO_MAP = {
  In_situ	          => 'MMO:0000658',
  Antibody	          => 'MMO:0000498',
  Western 	          => 'MMO:0000669',
  Reporter_gene	          => 'MMO:0000670',
  Cis_regulatory_element  => 'MMO:0000670',
  RT_PCR	          => 'MMO:0000655',
  EPIC	                  => 'MMO:0000670',
  Genome_Editing	  => 'MMO:0000670',
  Northern	          => 'MMO:0000647',
  UNDEFINED               => 'MMO:0000640',
};

my $TL_ANATOMY_TERM = 'WBbt:0000100';
my $TL_LIFESTAGE_TERM = 'WBls:0000075';

my $life_stage_term = {
  $TL_LIFESTAGE_TERM => 'Nematoda Life Stage'
};

my $anatomy_term = {
  $TL_ANATOMY_TERM =>  'C. elegans Cell and Anatomy'
};


my ($debug, $test, $verbose, $store, $wormbase);
my ($outfile, $acedbpath, $ws_version, $out_fh, $bgi_json);

GetOptions (
  "debug=s"     => \$debug,
  "test"        => \$test,
  "verbose"     => \$verbose,
  "store:s"     => \$store,
  "database:s"  => \$acedbpath,
  "outfile:s"   => \$outfile,
  "wsversion=s" => \$ws_version,
  "bgijson=s"   => \$bgi_json,
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
  $outfile = "./wormbase.agr_expression.${ws_version}.json";
}

#
# restrict to genes in the BGI, if provded
#
my (%only_these_genes, @expression_annots);

if (defined $bgi_json) {
  # read it so that we can restrict to genes that have basic info
  my %genes;

  my $json_string = "";
  open(my $json_fh, $bgi_json) or die "Could not open $bgi_json for reading\n";
  while(<$json_fh>) {
    $json_string .= $_;
  }

  my $json_reader = JSON->new;
  my $json = $json_reader->decode($json_string);
  
  foreach my $entry (@{$json->{data}}) {
    my $id = $entry->{primaryId};
    $id =~ s/WB://; 
    $only_these_genes{$id} = 1;
  }

}


my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or die("Connection failure: ". Ace->error);

my $it = $db->fetch_many(-query => 'find Expr_pattern WHERE Expression_data AND Gene AND Reference');

while (my $obj = $it->next) {
  next unless $obj->isObject();


  my ($paper) = $obj->Reference;
  my ($assay) = $obj->Type;

  if (not defined $assay or not exists $MMO_MAP->{$assay->name}) {
    $assay = $MMO_MAP->{UNDEFINED};
  } else {
    $assay = $MMO_MAP->{$assay->name};
  }
                                                   

  foreach my $gene_obj ($obj->Gene) {
    my $gene = $gene_obj->name;

    next if defined $bgi_json and not exists $only_these_genes{$gene};
    
    my %annots;

    foreach my $ls ($obj->Life_stage) {

      next if not &record_lifestage_term_name( $ls );

      my @ats;

      foreach my $thing ($ls->col()) {
        if ($thing->name eq 'Anatomy_term') {
          foreach my $at ($thing->col()) {
            next if not &record_anatomy_term_name($at);
            push @ats, $at;
          }
        }
      }
      if (@ats) {
        foreach my $at (@ats) {
          $annots{$ls->name}->{$at->name} = {};
        } 
      } else {
        $annots{$ls->name}->{$TL_ANATOMY_TERM} = {};
      }
    }

    foreach my $at ($obj->Anatomy_term) {

      next if not &record_anatomy_term_name( $at );

      my @lss;

      foreach my $thing ($at->col()) {
        if ($thing->name eq 'Life_stage') {
          foreach my $ls ($thing->col()) {
            next if not &record_lifestage_term_name($ls);
            push @lss, $ls;
          }
        }
      }
      if (@lss) {
        foreach my $ls (@lss) {
          $annots{$ls->name}->{$at->name} = {};
        } 
      } else {
        $annots{$TL_LIFESTAGE_TERM}->{$at->name} = {};
      }
    }

    foreach my $go ($obj->GO_term) {
      $annots{$TL_LIFESTAGE_TERM}->{$TL_ANATOMY_TERM}->{$go->name} = 1;
    }


    foreach my $ls (sort keys %annots) {
      foreach my $at (sort keys %{$annots{$ls}}) {

        my @where_expressed;

        if (keys %{$annots{$ls}->{$at}}) {
          foreach my $go_term (sort keys %{$annots{$ls}->{$at}}) {
            push @where_expressed, {
              anatomicalStructureTermId => $at,
              whereExpressedStatement => $anatomy_term->{$at},
              cellularComponentTermId => $go_term,
            };            
          }
        } else {
          push @where_expressed, {
            anatomicalStructureTermId => $at,
            whereExpressedStatement => $anatomy_term->{$at},
          };            
        }

        foreach my $whereEx (@where_expressed) {
          my $json_obj = {
            geneId => "WB:$gene",
            evidence => &get_paper_json( $paper ),
            assay => $assay,
            dateAssigned => $date,
            crossReference => {
              id => "WB:$obj",
              pages => ["gene/expression/annotation/detail"],
            },
            whenExpressedStage => $ls,
            wildtypeExpressionTermIdentifiers => $whereEx,
          };

          push @expression_annots, $json_obj;
        }
      }
    }
  }
}


my $meta_data = {
  dateProduced => $date,
  dataProvider => [
    { 
      crossReference => {
        id => "WB",
        pages => ["homepage"],
      },
      type => "curated",
    },
      ],
  release      => (defined $ws_version) ? $ws_version : $wormbase->get_wormbase_version_name(),
};

my $data = {
  metaData => $meta_data,
  data     => \@expression_annots,
};

if (defined $outfile) {
  open $out_fh, ">$outfile" or die "Could not open $outfile for writing\n";
} else {
  $out_fh = \*STDOUT;
}


my $json_obj = JSON->new;
my $string = $json_obj->allow_nonref->canonical->pretty->encode($data);
print $out_fh $string;

$db->close;

exit(0);


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
  $json_paper->{modPublicationId} = "WB:$wb_paper";
  if ($pmid) {
    $json_paper->{pubMedId} = "PMID:$pmid";
  }
  
  return $json_paper;
}

###############################################
sub record_anatomy_term_name {
  my $at = shift;

  if (not exists $anatomy_term->{$at->name}) {
    if (not $at->Term) {
      warn "Term field not populated for Anatomy_term $at; skipping\n";
    } else {
      $anatomy_term->{$at->name} = $at->Term->name;
    }
  }

  return exists($anatomy_term->{$at->name}) ? 1 : 0;
}

###############################################
sub record_lifestage_term_name {
  my $ls = shift;

  if (not exists $life_stage_term->{$ls->name}) {
    $life_stage_term->{$ls->name} = $ls->Public_name->name;
  }

  return (exists $life_stage_term->{$ls->name}) ? 1 : 0;
}
