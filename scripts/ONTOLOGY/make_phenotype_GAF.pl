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
my $full_name = $wormbase->full_name;

$acedbpath = $wormbase->autoace unless $acedbpath;
$outfile = $wormbase->ontology."/phenotype_association.".$wormbase->get_wormbase_version_name.".wb" unless $outfile;

$log->write_to("Connecting to database $acedbpath\n");
my $db = Ace->connect(-path => $acedbpath,  -program => $tace) or $log->log_and_die("Connection failure: ". Ace->error);

my ($count, $it);
my $gene_info = {};

$gene_info = &get_gene_info( $acedbpath, $wormbase );

$log->write_to( "Got name information for " . scalar(keys %$gene_info) . " genes\n");

open(my $out, ">$outfile") or $log->log_and_die("cannot open $outfile : $!\n");

&print_wormbase_GAF_header($out, $wormbase->get_wormbase_version_name);

$it = $db->fetch_many(-query=>'find Variation WHERE (Phenotype OR Phenotype_not_observed)');

my (%g2v_via_var, %gpvs_with_pub, %gpvs_with_person, %taxon_ids); 

while (my $obj=$it->next) {
  next unless $obj->isObject();
  next unless $obj->Species;

  my (@affected_genes, %pheno);
  
  my $var = $obj->name;

  @affected_genes = map { $_->name } $obj->Gene;
  next if not @affected_genes;

  for my $g ($obj->Gene) {
      $taxon_ids{$g->name} = $obj->Species->NCBITaxonomyID;
  }
  
  foreach my $key ('Phenotype', 'Phenotype_not_observed') {
    foreach my $pobj ($obj->$key) {
      my ($pheno) = $pobj->name =~ /(WBPhenotype:\d+)/;
      next unless $pheno;

      foreach my $thing ($pobj->col) {
        if ($thing->name eq 'Paper_evidence') {
          my $paper = $thing->right->name;
          foreach my $g (@affected_genes) {
            $g2v_via_var{$g}->{$key}->{$pheno}->{papers}->{$paper}->{$var} = 1;
            $gpvs_with_pub{$key}->{$g}->{$pheno}->{$var} = 1;
          }
        } elsif ($thing->name eq 'Person_evidence') {
          my $person = $thing->right->name;
          foreach my $g (@affected_genes) {
            $g2v_via_var{$g}->{$key}->{$pheno}->{persons}->{$var}->{$person} = 1;
            $gpvs_with_person{$key}->{$g}->{$pheno}->{$var} = 1;
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
				   $taxon_ids{$g},
                                   $date);
          $count++;
        }
      } 
      if (exists $href->{persons}) {
        foreach my $var (sort keys %{$href->{persons}}) {
          #
          # only include person-based evidence for this gene/phennotype/var if there was no paper-based one
          #
          if (not exists $gpvs_with_pub{$key}->{$g}->{$phen}->{$var}) {
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
				     $taxon_ids{$g},
                                     $date);
            $count++;
          }  
        }
      } 
      if (exists $href->{curators}) {
        #
        # only include curator-based evidence for this gene/pheno/var if not supported by a paper or person
        #
        foreach my $var (sort keys %{$href->{curators}}) {
          if (not exists $gpvs_with_pub{$key}->{$g}->{$phen}->{$var} and 
              not exists $gpvs_with_person{$key}->{$g}->{$phen}->{$var}) {
            
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
				     $taxon_ids{$g},
                                     $date);
            $count++;
          }
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
  
  my (@affected_genes, %pheno);

  foreach my $g ($obj->Gene) {
    next if not exists $gene_info->{$g->name};
    next if $gene_info->{$g->name}->{status} eq 'Dead';

    if ($g->right(2) eq 'RNAi_primary') {
      push @affected_genes, $g;
      $taxon_ids{$g->name} = $obj->Species->NCBITaxonomyID;
    }
  }

  # For RNAi, there should always be at most one paper as the source.
  # If more than one paper is attached, use the first one.
  # If no papers are attached, use the WB RNAi id as the primary source,
  # with the WBPerson id as the WITH/FROM
   
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
			       $taxon_ids{$g},
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
			       $taxon_ids{$g},
                               $date );
      $count++;
    }
  }
}

close($out);

&make_species_files($wormbase, $outfile);

$log->write_to("Printed $count phenotype lines RNAi objects\n");

$db->close;
$log->mail;

exit(0);

	

