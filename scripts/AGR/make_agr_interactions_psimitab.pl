#!/usr/bin/env perl

# probably from: https://github.com/HUPO-PSI/miTab/blob/master/PSI-MITAB27Format.md

use strict;
use Storable;	
use Getopt::Long;
use Ace;
use JSON;

use lib $ENV{CVS_DIR};
use Wormbase;
use Modules::AGR;


my $DETECTION_METHOD_MAPPING =  {
  Affinity_capture_luminescence	         => 'psi-mi:"MI:0004"(affinity chromatography technology)',
  Affinity_capture_MS                    => 'psi-mi:"MI:0004"(affinity chromatography technology)',
  Affinity_capture_Western	         => 'psi-mi:"MI:0004"(affinity chromatography technology)',
  Biochemical_activity	                 => 'psi-mi:"MI:0415"(enzymatic study)',
  Chromatin_immunoprecipitation          => 'psi-mi:"MI:0402"(chromatin immunoprecipitation assay)',
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

my $SRC_DB_MAPPING = {
  wormbase => 'psi-mi:"MI:0487"(wormbase)',
 # BioGRID commented out for now - until AGR loader ready for multi-source interactions
 # BioGRID  => 'psi-mi:"MI:0463"(biogrid)',
};

my $MOL_TYPE_MAPPING = {
  Protein => 'psi-mi:"MI:0326"(protein)',
  DNA     => 'psi-mi:"MI:0319"(deoxyribonucleic acid)',
  RNA     => 'psi-mi:"MI:0320"(ribonucleic acid)',
};

my $ROLE_MAPPING = {
  Bait            => 'psi-mi:"MI:0496"(bait)',
  Target          => 'psi-mi:"MI:0498"(prey)',
  Non_directional => 'psi-mi:"MI:0497"(neutral component)'
}; 

my ($outfile, $acedbpath, $bgi_json,$build);

GetOptions(
  "database:s" => \$acedbpath,
  "outfile:s"  => \$outfile,
  "bgijson=s"  => \$bgi_json,
  "build"      => \$build,
  )or die(@!);

die "You must supply both -database and -outfile\n" if not defined $acedbpath or not defined $outfile;

my $bgi_genes = &get_bgi_genes($bgi_json) if (defined $bgi_json);

open(my $out_fh, ">$outfile") or die "Could not open $outfile for writing\n";

my $db = Ace->connect(-path => $acedbpath) or die("Connection failure: ". Ace->error);

my $it = $db->fetch_many(-query => 'find Interaction Physical');

INT: while (my $obj = $it->next) {
  next unless $obj->isObject();

  my ($int_type) = $obj->at('Interaction_type.Physical');
  
  # ignore interactions with unspecified type
  if (not $int_type) {
    warn("Skipping $obj - unspecified DNA/Protein\n");
    next;
  }

  my (%genes, %methods);
  foreach my $meth ($obj->Detection_method) {
    $methods{$meth->name} = 1;
  }

  # Ignore interactions with only Antibody detection method
  delete $methods{Antibody} if exists $methods{Antibody};
  if (not keys %methods) {
    warn("Skipping $obj - no non-Antibody detection method\n");
    next;
  }


  foreach my $g ($obj->Interactor_overlapping_gene) {
    my $gn      = $g->name;
    my $gnp     = $g->Public_name->name;
    my $g_tax   = $g->Species->NCBITaxonomyID->name;
    my $g_sp_nm = $g->Species->name;

    if ($g_sp_nm ne 'Caenorhabditis elegans') {
      warn "Skipping $obj - not dealing with non-C.elegans genes for now\n";
      next INT;
    }

    # taxid:6239(caeel)|taxid:6239(Caenorhabditis elegans)
    my $gsp = sprintf("taxid:%d(%s)|taxid:%d(%s)", $g_tax, "caeel", $g_tax, $g_sp_nm);

    $genes{$gn}->{id} = $gn;
    $genes{$gn}->{public_name} = $gnp;
    $genes{$gn}->{species} = $gsp;
    $genes{$gn}->{source}->{interactor_overlapping_gene} = 1;

    foreach my $tp ($g->right->col()) {
      $genes{$gn}->{roles}->{$tp} = 1;
    }
  }

  foreach my $f ($obj->Feature_interactor) {
    my @fg = $f->Associated_with_gene;
    
    # Skip interactions for which the Feature_interactor is not associated with exactly one gene
    if (scalar(@fg) != 1) {
      warn("Skipping $obj - Feature_interactor is not associated with exactly one gene\n");
      next INT;
    }

    my $fg = $fg[0]->name;
    my $fgp = $fg[0]->Public_name->name;
    my $fg_tax = $fg[0]->Species->NCBITaxonomyID->name;
    my $fg_sp_nm = $fg[0]->Species->name;

    if ($fg_sp_nm ne 'Caenorhabditis elegans') {
      warn "Skipping $obj - not dealing with non-C.elegans genes for now\n";
      next INT;
    }

    # taxid:6239(caeel)|taxid:6239(Caenorhabditis elegans)
    my $gsp = sprintf("taxid:%d(%s)|taxid:%d(%s)", $fg_tax, "caeel", $fg_tax, $fg_sp_nm);

    $genes{$fg}->{id} = $fg;
    $genes{$fg}->{public_name} = $fgp;
    $genes{$fg}->{species} = $gsp;
    $genes{$fg}->{source}->{feature_interactor} = 1;

    foreach my $tp ($f->right->col()) {
      $genes{$fg}->{roles}->{$tp} = 1;
    }
  }

  my ($id_a, $name_a, $sp_a, $type_a, $role_a, $id_b, $name_b, $type_b, $sp_b, $role_b);

  if (scalar(keys %genes) == 1) {
    # special case of self-interactions
    my ($o) = values %genes;

    if (defined $bgi_genes and not exists $bgi_genes->{"WB:$o->{id}"}) {
      warn "Skipping $obj - gene $o->{id} is not part of the live gene set for this release\n";
      next;
    }

    $id_a = $id_b = $o->{id};
    $name_a = $name_b = $o->{public_name};
    $sp_a = $sp_b = $o->{species};
    $role_a = $ROLE_MAPPING->{Bait};
    $role_b = $ROLE_MAPPING->{Target};

    if (exists $o->{roles}->{Bait} and exists $o->{roles}->{Target}) {
      if ($int_type eq 'ProteinDNA') {
        $type_a = $MOL_TYPE_MAPPING->{Protein};
        $type_b = $MOL_TYPE_MAPPING->{DNA};
      } elsif ($int_type eq 'ProteinProtein') {
        $type_a = $type_b = $MOL_TYPE_MAPPING->{Protein};
      } else {
        warn("Skipping $obj - has unexpected content (bad interaction type: $int_type)\n");
      }
    } else {
      warn("Skipping $obj - has unexpected content (Bait / Target of single gene)\n");
      next;
    }
  } elsif (scalar(keys %genes) == 2) {
    my ($obj_a, $obj_b) = values %genes;

    if (defined $bgi_genes and not exists $bgi_genes->{"WB:$obj_a->{id}"}) {
      warn "Skipping $obj - gene $obj_a->{id} is not part of the live gene set for this release\n";
      next;
    }
    if (defined $bgi_genes and not exists $bgi_genes->{"WB:$obj_b->{id}"}) {
      warn "Skipping $obj - gene $obj_b->{id} is not part of the live gene set for this release\n";
      next;
    }

    $id_a   = $obj_a->{id};
    $name_a = $obj_a->{public_name};
    $sp_a   = $obj_a->{species};

    $id_b   = $obj_b->{id};
    $name_b = $obj_b->{public_name};
    $sp_b   = $obj_b->{species};

    if (scalar(keys %{$obj_a->{roles}}) != 1 or scalar(keys %{$obj_b->{roles}}) != 1) {
      warn( "Skipping $obj - has unexpected content (multiple or missing roles for genes)\n");
      next;
    }

    # role
    if (exists $obj_a->{roles}->{Bait} and exists $obj_b->{roles}->{Target}) {
      $role_a = $ROLE_MAPPING->{Bait};
      $role_b = $ROLE_MAPPING->{Target};
    } elsif (exists $obj_b->{roles}->{Bait} and exists $obj_a->{roles}->{Target}) {
      $role_a = $ROLE_MAPPING->{Target};
      $role_b = $ROLE_MAPPING->{Bait};
    } elsif (exists $obj_a->{roles}->{Non_directional} and exists $obj_b->{roles}->{Non_directional}) {
      $role_a = $role_b = $ROLE_MAPPING->{Non_directional};
    } else {
      warn("Skipping $obj - could not unambigously determine roles for interactors\n");
      next;
    }

    # type
    if ($int_type eq 'ProteinProtein') {
      $type_a = $type_b = $MOL_TYPE_MAPPING->{Protein};
    } elsif ($int_type eq 'ProteinRNA') {
      if (exists $obj_a->{roles}->{Bait} and exists $obj_b->{roles}->{Target}) {
        $type_a = $MOL_TYPE_MAPPING->{RNA};
        $type_b = $MOL_TYPE_MAPPING->{Protein};
      } elsif (exists $obj_a->{roles}->{Target} and exists $obj_b->{roles}->{Bait}) {
        $type_b = $MOL_TYPE_MAPPING->{RNA};
        $type_a = $MOL_TYPE_MAPPING->{Protein};
      } else {
        warn("Skipping $obj - Could not unambiguously determine type (odd roles of gene pair)\n");
        next;
      }
    } elsif ($int_type eq 'ProteinDNA') {
      if (exists $obj_a->{roles}->{Bait} and exists $obj_b->{roles}->{Target}) {
        $type_a = $MOL_TYPE_MAPPING->{DNA};
        $type_b = $MOL_TYPE_MAPPING->{Protein};
        if (exists $obj_b->{source}->{feature_interactor}) {
          ($type_a, $type_b) = ($type_b, $type_a); 
        }
      } elsif (exists $obj_a->{roles}->{Target} and exists $obj_b->{roles}->{Bait}) {
        $type_b = $MOL_TYPE_MAPPING->{DNA};
        $type_a = $MOL_TYPE_MAPPING->{Protein};
        if (exists $obj_a->{source}->{feature_interactor}) {
          ($type_a, $type_b) = ($type_b, $type_a); 
        }
      } else {
        warn("Skipping $obj - Could not unambiguously determine type (odd roles of gene pair)\n");
        next;
      }
    }
  } else {
    warn("Skipping $obj - has unexpected content (wrong gene count)\n");
    next;
  }

  my ($pmid, $author_year) = &get_paper_stuff( $obj->Paper );

  my @src_db = ($SRC_DB_MAPPING->{wormbase});
  my @src_db_acc = ("wormbase:$obj");

  if ($obj->Database) {
    foreach my $db ($obj->Database) {
      if (exists $SRC_DB_MAPPING->{$db->name}) {
        my $psi_mi_name = $SRC_DB_MAPPING->{$db->name};
        my ($short_name) = $psi_mi_name =~ /\((\S+)\)/;

        my $db_acc = $db->right->right->name;
        unshift @src_db, $psi_mi_name;
        unshift @src_db_acc, "$short_name:$db_acc";
      }
    }
  }

  my @l = ('-') x 42;

  $l[0]  = "wormbase:$id_a";
  $l[1]  = "wormbase:$id_b";

  $l[4]  = "wormbase:$name_a(public_name)";
  $l[5]  = "wormbase:$name_b(public_name)";
  $l[6]  = join("|", map { $DETECTION_METHOD_MAPPING->{$_} } sort keys %methods);
  $l[7]  = $author_year;
  $l[8]  = "pubmed:$pmid";
  $l[8] .= '|wormbase:'.$obj->Paper if $build;
  $l[9]  = $sp_a;
  $l[10] = $sp_b;
  $l[11] = 'psi-mi:"MI:0914"(association)';
  $l[12] = join("|", @src_db);
  $l[13] = join("|", @src_db_acc);

  $l[18] = $role_a;
  $l[19] = $role_b;
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
  ($author_year) = $brief_cit =~ /^([^\(]+\(\d+\))/; # matches: Simske JS et al. (1996)

  foreach my $db ($pap->Database) {
    if ($db->name eq 'MEDLINE') {
      $pmid = $db->right->right->name;
      last;
    }
  }
  return ($pmid, $author_year);
}
