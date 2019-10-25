#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;


use lib $ENV{CVS_DIR};
use Wormbase;
use Log_files;

use lib "$ENV{CVS_DIR}/ONTOLOGY";
use GAF;

my ($help, $debug, $test, $verbose, $store, $wormbase);
my ($outfile, $acedbpath);

GetOptions ("help"       => \$help,
            "debug=s"    => \$debug,
	    "test"       => \$test,
	    "verbose"    => \$verbose,
	    "store:s"    => \$store,
	    "database:s" => \$acedbpath,
	    "output:s"   => \$outfile,
	    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
      );
}

my $tace = $wormbase->tace;
my $log  = Log_files->make_build_log($wormbase);
my $date = &get_GAF_date();
my $taxid = $wormbase->ncbi_tax_id;
my $full_name = $wormbase->full_name;

$acedbpath = $wormbase->autoace unless $acedbpath;
$outfile = $wormbase->ontology."/phenotype_association.".$wormbase->get_wormbase_version_name.".wb" unless $outfile;

$log->write_to("Connecting to database $acedbpath\n");
my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or $log->log_and_die("Connection failure: ". Ace->error);

my ($count, $it);
my $gene_info = {};

foreach my $species('Caenorhabidtis elegans','Caenorhabditis briggsae','Pristionchus pacificus'){
  $gene_info =  {%$gene_info , %{&get_gene_info( $acedbpath, $wormbase, $species)}};
}

$log->write_to( "Got name information for " . scalar(keys %$gene_info) . " genes\n");

open(my $out, ">$outfile") or $log->log_and_die("cannot open $outfile : $!\n");

&print_wormbase_GAF_header($out);

$it = $db->fetch_many(-query=>'find Variation (Phenotype OR Phenotype_not_observed)');


my %g2v_via_var;

while (my $obj=$it->next) {
  next unless $obj->isObject();
  next unless $obj->Species;
# next unless $obj->Species->name eq $full_name;

  my (@affected_genes, %pheno);
  
  my $var = $obj->name;

  @affected_genes = map { $_->name } $obj->Gene;
  next if not @affected_genes;
  
  foreach my $key ('Phenotype', 'Phenotype_not_observed') {
    foreach my $pobj ($obj->$key) {
      my ($pheno) = $pobj->name =~ /(WBPhenotype:\d+)/;
      next unless $pheno;

      foreach my $thing ($pobj->col) {
        if ($thing->name eq 'Paper_evidence') {
          my $paper = $thing->right->name;
          foreach my $g (@affected_genes) {
            $g2v_via_var{$g}->{$key}->{$pheno}->{papers}->{$paper}->{$var} = 1;
          }
        } elsif ($thing->name eq 'Person_evidence') {
          my $person = $thing->right->name;
          foreach my $g (@affected_genes) {
            $g2v_via_var{$g}->{$key}->{$pheno}->{persons}->{$var}->{$person} = 1;
          }
        } elsif ($thing->name eq 'Curator_confirmed') {
          my $curator = $thing->right->name;
          foreach my $g (@affected_genes) {
            $g2v_via_var{$g}->{$key}->{$pheno}->{curators}->{$var}->{$curator} = 1;
          }
        }
      }
    }
  }
}

foreach my $g (sort keys %g2v_via_var) {
  next if not exists $gene_info->{$g};
  next if $gene_info->{$g}->{status} eq 'Dead';
  
  foreach my $key (sort keys %{$g2v_via_var{$g}}) {
    my $qual = ($key eq 'Phenotype_not_observed') ? "NOT" : "";

    foreach my $phen (sort keys %{$g2v_via_var{$g}->{$key}}) {
      my $href = $g2v_via_var{$g}->{$key}->{$phen};

      if (exists $href->{papers}) {
        foreach my $paper (sort keys %{$href->{papers}}) {
          my @vars = sort keys %{$href->{papers}->{$paper}};
          my $with_from = join("|", map { "WB:$_" } @vars);

          &print_wormbase_GAF_line($out,  
                                   $g, 
                                   $gene_info->{$g}->{public_name}, 
                                   $qual,  
                                   $phen, 
                                   "WB_REF:" . $paper,
                                   "IMP", 
                                   $with_from, 
                                   "P",
                                   $gene_info->{$g}->{sequence_name},
                                   $taxid, 
                                   $date);
          $count++;
        }
      } elsif (exists $href->{persons}) {
        # only person evidence for this association; we are not allowed to have people as evdience in GAF.
        # We therefore denote the variation as evidence, and add the person as WITH/FROM (fudge 1)
        foreach my $var (sort keys %{$href->{persons}}) {
          my @persons = sort keys %{$href->{persons}->{$var}};
          my $with_from = join("|", map { "WB:$_" } @persons);
                               
          &print_wormbase_GAF_line($out, 
                                   $g, 
                                   $gene_info->{$g}->{public_name}, 
                                   $qual, 
                                   $phen, 
                                   "WB:" . $var,
                                   "IMP", 
                                   $with_from, 
                                   "P",
                                   $gene_info->{$g}->{sequence_name},
                                   $taxid, 
                                   $date);
          $count++;
        }  
      } elsif (exists $href->{curators}) {
        # only Curator_confirmed evidence for this association. These usually come from large-scale projects
        # like NBP. What we really need here is a project identifier in the evidence field, and to add the list
        # of variations in the WITH/FROM. However, since we do not have project identifiers, we therefore
        # again denote the variation as evidence and leave WITH/FROM empty (fudge 2)
        foreach my $var (sort keys %{$href->{curators}}) {
          &print_wormbase_GAF_line($out, 
                                   $g, 
                                   $gene_info->{$g}->{public_name}, 
                                   $qual, 
                                   $phen, 
                                   "WB:" . $var,
                                   "IMP", 
                                   "",
                                   "P",
                                   $gene_info->{$g}->{sequence_name},
                                   $taxid, 
                                   $date);
          $count++;
        }
      }
    }      
  }
}

$log->write_to("Printed $count phenotype lines for Variation objects\n");

$it = $db->fetch_many(-query=>'find RNAi (Phenotype OR Phenotype_not_observed)');
$count = 0;
while (my $obj = $it->next) {
  next unless $obj->isObject();
  next unless $obj->Species;
# next unless $obj->Species->name eq $full_name;
  
  my (@affected_genes, %pheno);

  foreach my $g ($obj->Gene) {
    next if not exists $gene_info->{$g->name};
    next if $gene_info->{$g->name}->{status} eq 'Dead';

    if ($g->right(2) eq 'RNAi_primary') {
      push @affected_genes, $g;
    }
  }

  #
  # For RNAi, there should always be at most one paper as the source.
  # If more than one paper is attached, use the first one.
  # If no papers are attached, use the WB RNAi id as the primary source,
  # with the WBPerson id as the WITH/FROM
  # 
  my ($ref, $with_from);

  my (@ref) = $obj->Reference;
  if (@ref) {
    $ref = "WB_REF:$ref[0]";
    $with_from = "WB:" . $obj->name;
  } else {
    my $evi = $obj->Evidence->right;
    if ($evi->name eq 'Person_evidence') {
      $ref = "WB:" . $obj->name;
      $with_from = "WB:" . $evi->right->name;
    }
  }

  next if not defined $ref;

  foreach my $p ($obj->Phenotype) {
    foreach my $g (@affected_genes) {
      &print_wormbase_GAF_line($out, 
                               $g, 
                               $gene_info->{$g}->{public_name}, 
                               "", 
                               $p, 
                               $ref,
                               "IMP", 
                               $with_from,
                               "P",  
                               $gene_info->{$g}->{sequence_name},
                               $taxid, 
                               $date );
      $count++;
    }
  }
  foreach my $p ($obj->Phenotype_not_observed) {
    foreach my $g (@affected_genes) {
      &print_wormbase_GAF_line($out, 
                               $g, 
                               $gene_info->{$g}->{public_name}, 
                               "NOT", 
                               $p, 
                               $ref,
                               "IMP", 
                               $with_from,
                               "P",  
                               $gene_info->{$g}->{sequence_name},
                               $taxid, 
                               $date );
      $count++;
    }
  }
}

$log->write_to("Printed $count phenotype lines RNAi objects\n");

$db->close;
$log->mail;

exit(0);

	

