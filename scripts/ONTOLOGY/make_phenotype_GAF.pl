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

my ($gene_info, $count, $it);

$gene_info = &get_gene_info( $acedbpath, $wormbase, $full_name );

$log->write_to( "Got name information for " . scalar(keys %$gene_info) . " genes\n");

open(my $out, ">$outfile") or $log->log_and_die("cannot open $outfile : $!\n");

$it = $db->fetch_many(-query=>'find Variation (Phenotype OR Phenotype_not_observed)');

while (my $obj=$it->next) {
  next unless $obj->isObject();
  next unless $obj->Species;
  next unless $obj->Species->name eq $full_name;

  $count++;
  if ($count % 1000 == 0) {
    warn "$count Variation objects processed\n";
  }
  
  my (@affected_genes, %pheno);
  
  my $var = $obj->name;
  
  @affected_genes = map { $_->name } $obj->Gene;
  next if not @affected_genes;
  
  foreach my $key ('Phenotype', 'Phenotype_not_observed') {
    foreach my $pobj ($obj->$key) {
      my ($pheno) = $pobj->name =~ /(WBPhenotype:\d+)/;
      next unless $pheno;

      $pheno{$pheno}->{$key} = {};

      foreach my $thing ($pobj->right) {
        if ($thing->name eq 'Paper_evidence') {
          my $paper = $thing->right->name;
          
          $pheno{$pheno}->{$key}->{$paper} = 1;
        }
      }
    }
  }
  
  foreach my $g (@affected_genes) {
    next if not exists $gene_info->{$g};
    next if $gene_info->{$g}->{status} eq 'Dead';

    foreach my $pheno (keys %pheno) {
      if (exists $pheno{$pheno}->{Phenotype}) {
        my @papers = keys %{$pheno{$pheno}->{Phenotype}};
        
        &print_wormbase_GAF_line($out,  
                                 $g, 
                                 $gene_info->{$g}->{public_name}, 
                                 "",  
                                 $pheno, 
                                 join("|", map { "WB_REF:$_" } @papers), 
                                 "Variation", 
                                 $var, 
                                 "P",
                                 $gene_info->{$g}->{sequence_name},
                                 $taxid, 
                                 $date);
      }
      
      if (exists $pheno{$pheno}->{Phenotype_not_observed}) {
        my @papers = keys %{$pheno{$pheno}->{Phenotype_not_observed}};
        
        &print_wormbase_GAF_line($out, 
                                 $g, 
                                 $gene_info->{$g}->{public_name}, 
                                 "NOT", 
                                 $pheno, 
                                 join("|", map { "WB_REF:$_" } @papers), 
                                 "Variation", 
                                 $var, 
                                 "P",
                                 $gene_info->{$g}->{sequence_name},
                                 $taxid, 
                                 $date);
      }
    }      
  }
}

$log->write_to("Printed phenotype lines for $count Variation objects\n");

$it = $db->fetch_many(-query=>'find RNAi (Phenotype OR Phenotype_not_observed)');
$count = 0;
while (my $obj = $it->next) {
  next unless $obj->isObject();
  next unless $obj->Species;
  next unless $obj->Species->name eq $full_name;

  $count++;
  if ($count % 10000 == 0) {
    warn "$count RNAi objects processed\n";
  }
  
  my (@affected_genes, %pheno);

  foreach my $g ($obj->Gene) {
    next if not exists $gene_info->{$g->name};
    next if $gene_info->{$g->name}->{status} eq 'Dead';

    if ($g->right(2) eq 'RNAi_primary') {
      push @affected_genes, $g;
    }
  }

  my @ref = $obj->Reference;
  
  foreach my $p ($obj->Phenotype) {
    foreach my $g (@affected_genes) {
      &print_wormbase_GAF_line($out, 
                               $g, 
                               $gene_info->{$g}->{public_name}, 
                               "", 
                               $p, 
                               join("|", map { "WB_REF:$_" } @ref),
                               "RNAi", 
                               $obj->name, 
                               "P",  
                               $gene_info->{$g}->{sequence_name},
                               $taxid, 
                               $date );
    }
  }
  foreach my $p ($obj->Phenotype_not_observed) {
    foreach my $g (@affected_genes) {
      &print_wormbase_GAF_line($out, 
                               $g, 
                               $gene_info->{$g}->{public_name}, 
                               "NOT", 
                               $p, 
                               join("|", map { "WB_REF:$_" } @ref), 
                               "RNAi", 
                               $obj->name, 
                               "P",  
                               $gene_info->{$g}->{sequence_name},
                               $taxid, 
                               $date );
    }
  }
}

$log->write_to("Printed phenotype line sfor $count RNAi objects\n");

$db->close;
$log->mail;

exit(0);

	

