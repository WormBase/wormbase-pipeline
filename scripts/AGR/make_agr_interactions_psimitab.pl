#!/usr/bin/env perl

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;

use lib $ENV{CVS_DIR};
use Wormbase;


my $DETECTION_METHOD_MAPPING =  {
  Affinity_capture_luminescence	         => 'psi-mi:"MI:0004"(affinity chromatography technology)',
  Affinity_capture_MS                    => 'psi-mi:"MI:0004"(affinity chromatography technology)',
  Affinity_capture_Western	         => 'psi-mi:"MI:0004"(affinity chromatography technology)',
  Biochemical_activity	                 => 'psi-mi:"MI:0415"(enzymatic study)',
  Chromatin_immunoprecipitation          => 'psi-mi:"MI:0402"(chromatin immunoprecipitation assay(MI:0402)',
  Cocrystal_structure	                 => 'psi-mi:"MI:0114"(x-ray crystallography)',
  Cofractionation	                 => 'psi-mi:"MI:0401"(biochemical)',
  Colocalization	                 => 'psi-mi:"MI:0428"(imaging techniques)',
  Copurification	                 => 'psi-mi:"MI:0401"(biochemical)',
  Directed_yeast_one_hybrid	         => 'psi-mi:"MI:0432"(one hybrid)',
  DNase_I_footprinting	                 => 'psi-mi:"MI:0606"(DNase I footprinting)',
  Electrophoretic_mobility_shift_assay   => 'psi-mi:"MI:0413"(electrophoretic mobility shift assay)',
  Far_western	                         => 'psi-mi:"MI:0047"(far western blotting)',
  Fluorescence_resonance_energy_transfer => 'psi-mi:"MI:0055"(fluorescent resonance energy transfer)',
  Protein_fragment_complementation_assay => 'psi-mi:"MI:0090"(protein complementation assay)',
  Protein_peptide	                 => 'psi-mi:"MI:0686"(unspecified method)',
  Protein_RNA	                         => 'psi-mi:"MI:0686"(unspecified method)',
  Reconstituted_complex	                 => 'psi-mi:"MI:0045"(experimental interaction detection)',
  Yeast_two_hybrid	                 => 'psi-mi:"MI:0018"(two hybrid)',
  Yeast_one_hybrid	                 => 'psi-mi:"MI:0432"(one hybrid)',
  #Antibody                               => '',
};

my $psi_mi_prot = 'psi-mi:"MI:0326"(protein)';
my $psi_mi_dna = 'psi-mi:"MI:0319"(deoxyribonucleic acid)';
my $psi_mi_rna = 'psi-mi:"MI:0320"(ribonucleic acid)';


my ($debug, $test, $verbose, $store, $wormbase);
my ($outfile, $acedbpath, $ws_version, $outfh);

GetOptions (
  "debug=s"     => \$debug,
  "test"        => \$test,
  "verbose"     => \$verbose,
  "store:s"     => \$store,
  "database:s"  => \$acedbpath,
  "outfile:s"   => \$outfile,
    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug   => $debug,
                             -test    => $test,
      );
}


$acedbpath = $wormbase->autoace unless $acedbpath;
if (not defined $outfile) {
  $outfile = "./physical_interactions.".$wormbase->get_wormbase_version_name.".psimitab.txt";
}

my $db = Ace->connect(-path => $acedbpath, -program => $wormbase->tace) or die("Connection failure: ". Ace->error);

open(my $out_fh, ">$outfile") or die "Could not open $outfile for writing\n";

my $it = $db->fetch_many(-query => 'find Interaction Physical');
#$it = $db->fetch_many(-query => 'find Interaction WBInteraction000505215');

INT: while (my $obj = $it->next) {
  next unless $obj->isObject();

  my ($int_type) = $obj->at('Interaction_type.Physical');
  #
  # ignore interactions with unspecified type
  #
  if (not $int_type) {
    warn("Skipping $obj - unspecified DNA/Protein\n");
    next;
  }

  #
  # ignore ProteinRNA interactions (for now)
  #
  if ($int_type eq 'ProteinRNA') {
    warn("Skipping $obj - ProteinRNA\n");
    next;
  }


  my (%genes, %methods);
  foreach my $meth ($obj->Detection_method) {
    $methods{$meth->name} = 1;
  }

  #
  # Ignore interactions with only Antibody detection method
  #
  delete $methods{Antibody} if exists $methods{Antibody};
  if (not keys %methods) {
    warn("Skipping $obj - no non-Antibody detection method\n");
    next;
  }


  foreach my $g ($obj->Interactor_overlapping_gene) {
    my $gn = $g->name;
    my $gnp = $g->Public_name->name;
    my $gsp = sprintf("taxid:%d(%s)", $g->Species->NCBITaxonomyID->name,, $g->Species->name);

    $genes{$gn}->{id} = $gn;
    $genes{$gn}->{public_name} = $gnp;
    $genes{$gn}->{species} = $gsp;

    foreach my $tp ($g->right->col()) {
      $genes{$gn}->{roles}->{$tp} = 1;
    }
  }

  foreach my $f ($obj->Feature_interactor) {
    my @fg = $f->Associated_with_gene;
    #
    # Skip interactions for which the Feature_interactor is not associated with exactly one gene
    #
    if (scalar(@fg) != 1) {
      warn("Skipping $obj - Feature_interactor is not associated with exactly one gene\n");
      next INT;
    }

    my $fg = $fg[0]->name;
    my $fgp = $fg[0]->Public_name->name;
    my $gsp = sprintf("taxid:%d(%s)", $fg[0]->Species->NCBITaxonomyID->name, $fg[0]->Species->name);

    $genes{$fg}->{id} = $fg;
    $genes{$fg}->{public_name} = $fgp;
    $genes{$fg}->{species} = $gsp;

    foreach my $tp ($f->right->col()) {
      $genes{$fg}->{roles}->{$tp} = 1;
    
    }
  }

  my ($id_a, $name_a, $sp_a, $type_a, $id_b, $name_b, $type_b, $sp_b);

  if (scalar(keys %genes) == 1) {
    # special case of self-interactions
    my ($o) = values %genes;

    $id_a = $id_b = $o->{id};
    $name_a = $name_b = $o->{public_name};

    if (exists $o->{roles}->{Bait} and exists $o->{roles}->{Target}) {
      if ($int_type eq 'ProteinDNA') {
        $type_a = $psi_mi_prot;
        $type_b = $psi_mi_dna;
      } elsif ($int_type eq 'ProteinProtein') {
        $type_a = $psi_mi_prot;
        $type_b = $psi_mi_prot;
      } else {
        warn("Skipping  $obj - has unexpected content (bad interaction type: $int_type)\n");
      }
    } else {
      warn("Skipping $obj - has unexpected content (Bait / Target of single gene)\n");
      next;
    }
  } elsif (scalar(keys %genes) == 2) {
    my ($obj_a, $obj_b) = values %genes;

    $id_a   = $obj_a->{id};
    $name_a = $obj_a->{public_name};
    $sp_a   = $obj_a->{species};

    $id_b   = $obj_b->{id};
    $name_b = $obj_b->{public_name};
    $sp_b   = $obj_b->{species};

    if ($int_type eq 'ProteinProtein') {
      $type_a = $psi_mi_prot;
      $type_b = $psi_mi_prot;
    } elsif ($int_type eq 'ProteinDNA') {
      if (scalar(keys %{$obj_a->{roles}}) != 1 or scalar(keys %{$obj_b->{roles}}) != 1) {
        warn( "Skipping $obj - has unexpected content (multiple of missing roles for genes)\n");
        next;
      }
      if (exists $obj_a->{roles}->{Bait} and exists $obj_b->{roles}->{Target}) {
        $type_a = $psi_mi_dna;
        $type_b = $psi_mi_prot;
      } elsif (exists $obj_a->{roles}->{Target} and exists $obj_b->{roles}->{Bait}) {
        $type_b = $psi_mi_dna;
        $type_a = $psi_mi_prot;
      } else {
        warn("Skipping $obj - has unexpected content (odd roles of gene pair)\n");
        next;
      }

    }
  } else {
    warn("Skipping $obj - has unexpected content (wrong gene count)\n");
    next;
  }

  my ($pmid, $author_year) = &get_paper_stuff( $obj->Paper );

  my @l = ('-') x 42;

  $l[0]  = "wormbase:$id_a";
  $l[1]  = "wormbase:$id_b";
  $l[4]  = "wormbase:$name_a(public_name)";
  $l[5]  = "wormbase:$name_b(public_name)";
  $l[6]  = join("|", map { $DETECTION_METHOD_MAPPING->{$_} } sort keys %methods);
  $l[7]  = $author_year;
  $l[8]  = "pubmed:$pmid";
  $l[9]  = $sp_a;
  $l[10] = $sp_b;
  $l[11] = 'psi-mi:"MI:0914"(association)';
  $l[12] = 'psi-mi:"MI:0487"(wormbase)';
  $l[13] = "wormbase:$obj";

  $l[20] = $type_a;
  $l[21] = $type_b;

             
  print $out_fh join("\t", @l), "\n";
}

$db->close();
close($out_fh);

exit(0);

#####################################################################

sub get_paper_stuff {
  my ($pap) = @_;

  my ($author_year, $pmid);

  my $brief_cit = $pap->Brief_citation->name;
  ($author_year) = $brief_cit =~ /^([^\(]+\(\d+\))/; 

  foreach my $db ($pap->Database) {
    if ($db->name eq 'MEDLINE') {
      $pmid = $db->right->right->name;
      last;
    }
  }

  return ($pmid, $author_year);
}
