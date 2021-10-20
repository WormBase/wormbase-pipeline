#!/usr/bin/env perl
#
#
# Dumps protein motifs from ensembl mysql (protein) database in GFF3 format
#
# Last updated by: $Author: klh $
#

use lib $ENV{'CVS_DIR'};

use strict;
use DBI;
use Ace;
use Getopt::Long;
use Wormbase;
use Storable;
use Log_files;

use Bio::EnsEMBL::Registry;

my ($WPver, @methods);
my ($store, $test, $debug, $species, $ens_regconf, $outfile, $out_fh, $acedb);
my (%cds2wp, %cds2gene, %cds2cgc, $testcds,  $motifs, $homologies);

GetOptions(
  "methods=s"   => \@methods,
  "store:s"     => \$store,
  "test"        => \$test,
  "debug:s"     => \$debug,
  "species:s"   => \$species,
  "outfile=s"   => \$outfile,
  'ensreg=s'    => \$ens_regconf,
  'regconf=s'   => \$ens_regconf,
  "testcds=s"   => \$testcds,
  "acedb=s"     => \$acedb, 
	  );


# define the names of the methods to be dumped
my %features = (
  tmhmm       => ['TMHMM',       'transmembrane_helix'],
  signalp     => ['SignalP',     'signal_peptide'],
  seg         => ['Seg',         'compositionally_biased_region_of_peptide'],
  ncoils      => ['ncoils',      'coiled_coil'],
  mobidblite  => ['MobiDB-lite', 'intrinsically_unstructured_polypeptide_region'],
    );


my %motifs = (
  hmmpanther  => ['PANTHER',     'HMMPANTHER'],
  pirsf       => ['PIRSF',       'PIRSF'],
  prints      => ['PRINTS',      'PRINTS'],
  scanprosite => ['PROSITE',     'SCANPROSITE'],
  smart       => ['SMART',       'SMART'],
  tigrfam     => ['TIGRFAMs',    'TIGRFAM'],
  pfam        => ['Pfam',        'PFAM'],
  superfamily => ['SUPERFAMILY', 'SUPERFAMILY']
    );
 
my %pfam_sites = (
  Active_site => "catalytic_residue",
  Metal_ion_binding_site => "metal_binding_site",
    );


my @post_trans_mods = (
  "acetylation",
  "Acetylation site",
  "dimethylation",
  "methylation",
  "N-Glycosylation site",
  "Oxidation site",
  "Phosphorylation site",
  "trimethylation"
    );


my $wormbase;
if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
  $wormbase = Wormbase->new( -debug    => $debug,
                             -test     => $test,
                             -organism => $species,
			     );
}
$species = $wormbase->species;
my $log = Log_files->make_build_log($wormbase);

#
# Get Ensembl data...
#
$log->write_to("Loading registry and fetching translations...\n");

my $reg = "Bio::EnsEMBL::Registry";
$reg->load_all($ens_regconf);
my $ens_tadap = $reg->get_adaptor($species, 'Core', 'Translation'); 

my @transls;
if (defined $testcds) {
  @transls = ($ens_tadap->fetch_by_stable_id($testcds));
} else {
  @transls = @{$ens_tadap->fetch_all};
}


#
# Get Wormbase feature data...
#
$log->write_to("Fetching WormBase feature data...\n");

$acedb = $wormbase->autoace if not defined $acedb;
$wormbase->FetchData('cds2wormpep', \%cds2wp, "$acedb/COMMON_DATA");
$wormbase->FetchData('cds2wbgene_id', \%cds2gene, "$acedb/COMMON_DATA");
$wormbase->FetchData('cds2cgc', \%cds2cgc, "$acedb/COMMON_DATA");

my $ptm_feats = &get_post_translational_feats_from_ace();
my $pfam_feats = &get_pfam_sites_from_ace();
my $mass_spec_feats = ($wormbase->species eq 'elegans') ? &get_mass_spec_homols_from_ace() : {};

if ($outfile) {
  if ($outfile eq '-') {
    $out_fh = \*STDOUT;
  } else {
    open $out_fh, ">$outfile" or $log->log_and_die("Could not open $outfile for writing\n");
  }
} else {
  $outfile = $wormbase->sequences . "/" . join(".", $wormbase->species, "protein_annotation", "gff3");
  open $out_fh, ">$outfile" or $log->log_and_die("Could not open $outfile for writing\n");
}

$log->write_to("Dumping data...\n");
&write_gff_header();

my $db = Ace->connect(-path => $acedb, -program => $wormbase->tace) or die ('Connection failure: ' . Ace->error);

while (my $transl = shift @transls) {
  my $wormpep = $cds2wp{$transl->stable_id};
  my $gene = $cds2gene{$transl->stable_id};
  my $gname = $transl->stable_id;
  my $uniprot_id = get_uniprot_id($wormpep, $db);
  $gname =~ s/[a-z]$//; 
  $gname = $cds2cgc{$transl->stable_id} if exists $cds2cgc{$transl->stable_id};

  #
  # get the exon structure and turn it into peptide coords
  #
  my %pep_attr = (ID => $transl->stable_id, wormbase_protein => $wormpep, wormbase_geneid => $gene, wormbase_genename => $gname);
  $pep_attr{uniprot_id} = $uniprot_id if $uniprot_id;
  &write_pep_gff(
      $transl->stable_id,
      "WormBase",
      "polypeptide", 
      1, 
      $transl->length,
      ".", 
      \%pep_attr
      );

  my %cds_attr = (ID => $transl->stable_id . "_exon_boundaries", wormbase_protein => $wormpep, wormbase_geneid => $gene, wormbase_genename => $gname);
  $cds_attr{uniprot_id} = $uniprot_id if $uniprot_id;
  &write_pep_gff(
      $transl->stable_id,
      "WormBase",
      "CDS", 
      1, 
      $transl->length,
      ".",
      \%cds_attr
    );

  my $tran = $transl->transcript;
  foreach my $seg ($transl->transcript->genomic2pep($tran->seq_region_start, $tran->seq_region_end, $tran->strand)) {
    next if $seg->isa("Bio::EnsEMBL::Mapper::Gap");
    &write_pep_gff(
       $transl->stable_id,
       "WormBase",
       "exon", 
       $seg->start,  
       $seg->end,
       ".", 
       { Parent => $transl->stable_id . "_exon_boundaries" });
    
  }


  #
  # PTMs 
  #
  if (exists $ptm_feats->{$wormpep}) {
    foreach my $f (sort { $a->{start} <=> $b->{start} } @{$ptm_feats->{$wormpep}}) {
      &write_pep_gff(
         $transl->stable_id,
         "Mass-spec",
         "post_translationally_modified_region",
         $f->{start}, 
         $f->{end},
         ".",
         { Name => $f->{feat} });         
    }
  }

  #
  # Pfam feats 
  #
  if (exists $pfam_feats->{$wormpep}) {
    foreach my $f (sort { $a->{start} <=> $b->{start} } @{$pfam_feats->{$wormpep}}) {
      &write_pep_gff(
         $transl->stable_id,
         $f->{source}, 
         $pfam_sites{$f->{feat}},
         $f->{start}, 
         $f->{end},
         ".",
         { Name => $pfam_sites{$f->{feat}} } );         
    }
  }

  #
  # mass-spec homols
  #
  if (exists $mass_spec_feats->{$wormpep}) {
    foreach my $f (sort { $a->{start} <=> $b->{start} } @{$mass_spec_feats->{$wormpep}}) {
      &write_pep_gff(
         $transl->stable_id,
         "Mass_spec_peptide",
         "match",
         $f->{start}, 
         $f->{end},
         ".",
         { Name => $f->{peptide}, wormbase_protein => $f->{peptide} } );
    }
  }         


  #
  # features and motifs
  #

  foreach my $logic (sort keys %features) {
    foreach my $pf (sort { $a->start <=> $b->start} @{$transl->get_all_ProteinFeatures($logic)}) {
      &write_pep_gff(
         $transl->stable_id,
         $features{$logic}->[0],
         $features{$logic}->[1],
         $pf->start,  
         $pf->end,
         $pf->score,
         { Name =>$features{$logic}->[0] });
    }
  }
    

  foreach my $logic (sort keys %motifs) {
    foreach my $pf (sort { $a->start <=> $b->start} @{$transl->get_all_ProteinFeatures($logic)}) {
      my $name = sprintf("%s %s (%s)", $motifs{$logic}->[0], $pf->hseqname, $pf->hdescription);
      my $wb_obj = sprintf("%s:%s", $motifs{$logic}->[1], $pf->hseqname);
      &write_pep_gff(
        $transl->stable_id,
        $motifs{$logic}->[0],
        "motif", 
        $pf->start,  
        $pf->end,
        $pf->score,
        { Name => $name, wormbase_motif =>$wb_obj });
    }
  }  

  print $out_fh "###\n";

}

close($out_fh) or $log->log_and_die("Could not cleanly close output file $outfile\n");
$db->close();
$log->mail();

exit(0);

#################################


sub write_gff_header {
  print $out_fh "##gff-version 3\n";
}

sub write_pep_gff {
  my ($seq, $source, $feat, $st, $en, $score, $group) = @_;

  printf($out_fh "%s\t%s\t%s\t%d\t%d\t%s\t.\t.\t", $seq, $source, $feat, $st, $en, $score ? $score : ".");
  my @attr;
  foreach my $k (sort keys %$group) {
    push @attr, sprintf("%s=%s", $k, $group->{$k});
  }
  print $out_fh join(";", @attr);
  print $out_fh "\n";
}


sub get_post_translational_feats_from_ace {
  
  my %ptm_feats;

  my $tace = $wormbase->tace;

  my $def = &ptm_table_maker_def(@post_trans_mods);
  my $tm_cmd = "Table-maker -p \"$def\"\nquit\n";
  
  open(my $tace_fh, "echo '$tm_cmd' | $tace $acedb |");
  while(<$tace_fh>) {
    /^\"(\S+)\"\s+\"(.+)\"\s+(\d+)\s+(\d+)/ and do {
      my ($protein, $feature, $start, $end) = ($1, $2, $3, $4);
      $protein =~ s/^\S+://; 
      push @{$ptm_feats{$protein}}, {
        feat => $feature, 
        start => $start, 
        end => $end, 
      };
    }      
  }
  
  unlink $def;

  #foreach my $prot (sort keys %ptm_feats) {
  #  print STDERR "Prot: $prot\n";
  #  foreach my $f ( sort { $a->{feat} cmp $b->{feat} or $a->{start} <=> $b->{start} } @{$ptm_feats{$prot}}) {
  #    print STDERR "  $f->{feat} $f->{start} $f->{end}\n";
  #  }
  #}

  return \%ptm_feats;
}

sub get_pfam_sites_from_ace {
  
  my %pfam_feats;

  my $tace = $wormbase->tace;

  my $def = &pfam_sites_table_maker_def(keys %pfam_sites);
  my $tm_cmd = "Table-maker -p \"$def\"\nquit\n";
  
  open(my $tace_fh, "echo '$tm_cmd' | $tace $acedb |");
  while(<$tace_fh>) {
    /^\"(\S+)\"\s+\"(.+)\"\s+\"(\S+)\"\s+(\d+)\s+(\d+)\s+(\d+)/ and do {
      my ($protein, $feature, $method, $score, $start, $end) = ($1, $2, $3, $4, $5, $6);
      $protein =~ s/^\S+://; 
      push @{$pfam_feats{$protein}}, {
        feat   => $feature, 
        source => $method,
        start  => $start, 
        end    => $end, 
      };
    }      
  }
  
  unlink $def;

  #foreach my $prot (sort keys %pfam_feats) {
  #  print STDERR "Prot: $prot\n";
  #  foreach my $f ( sort { $a->{feat} cmp $b->{feat} or $a->{start} <=> $b->{start} } @{$pfam_feats{$prot}}) {
  #    print STDERR "  $f->{feat} $f->{source} $f->{start} $f->{end}\n";
  #  }
  #}

  return \%pfam_feats;
}

sub get_mass_spec_homols_from_ace {
  my %mass_spec_homols;

  my $tace = $wormbase->tace;

  my $def = &mass_spec_table_maker_def();
  my $tm_cmd = "Table-maker -p \"$def\"\nquit\n";
  
  open(my $tace_fh, "echo '$tm_cmd' | $tace $acedb |");
  while(<$tace_fh>) {
    /^\"(\S+)\"\s+\"(.+)\"\s+\"(\S+)\"\s+(\d+)\s+(\d+)\s+(\d+)/ and do {
      my ($protein, $target, $method, $score, $start, $end) = ($1, $2, $3, $4, $5, $6);
      $protein =~ s/^\S+://; 
      push @{$mass_spec_homols{$protein}}, {
        peptide   => $target, 
        start  => $start, 
        end    => $end, 
      };
    }      
  }
  
  unlink $def;

  return \%mass_spec_homols;
}

#############

sub ptm_table_maker_def {
  my (@mods) = @_;

  my $tm_def_file = "/tmp/wb_protein_features.ptm.$$.def";
  open(my $tm_def_fh, ">$tm_def_file") or $log->log_and_die("Could not open $tm_def_file\n");

  my $cond1 = join(" OR ", map { "Feature = \"$_\"" } @mods); 
  my $cond2 = join(" OR ", map { "\"$_\"" } @mods); 

  my $query = <<"END";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Protein 
From 1 

Condition $cond1
 
Colonne 2 
Width 12 
Optional 
Visible 
Class 
Class Method 
From 1 
Tag Feature 
Condition $cond2
 
Colonne 3 
Width 12 
Optional 
Visible 
Integer 
Right_of 2 
Tag  HERE  
 
Colonne 4 
Width 12 
Optional 
Visible 
Integer 
Right_of 3 
Tag  HERE  

END

  print $tm_def_fh "$query\n";
  close($tm_def_fh);

  return $tm_def_file;

}


sub pfam_sites_table_maker_def {
  my (@pfam_sites) = @_;

  my $tm_def_file = "/tmp/wb_protein_features.pfsite.$$.def";
  open(my $tm_def_fh, ">$tm_def_file") or $log->log_and_die("Could not open $tm_def_file\n");

  my $cond1 = join(" OR ", map { "Motif_homol = \"$_\"" } @pfam_sites); 
  my $cond2 = join(" OR ", map { "\"$_\"" } @pfam_sites); 

  my $query = <<"END";

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Protein 
From 1 
Condition $cond1
 
Colonne 2 
Width 12 
Optional 
Visible 
Class 
Class Motif 
From 1 
Tag Motif_homol 
Condition $cond2
 
Colonne 3 
Width 12 
Optional 
Visible 
Class 
Class Method 
Right_of 2 
Tag  HERE  
 
Colonne 4 
Width 12 
Optional 
Visible 
Float 
Right_of 3 
Tag  HERE  
 
Colonne 5 
Width 12 
Optional 
Visible 
Integer 
Right_of 4 
Tag  HERE  
 
Colonne 6 
Width 12 
Optional 
Visible 
Integer 
Right_of 5 
Tag  HERE  
 
 
END

  print $tm_def_fh "$query\n";
  close($tm_def_fh);

  return $tm_def_file;

}


sub mass_spec_table_maker_def {
  my (@pfam_sites) = @_;

  my $tm_def_file = "/tmp/wb_protein_features.mass_spec.$$.def";
  open(my $tm_def_fh, ">$tm_def_file") or $log->log_and_die("Could not open $tm_def_file\n");

  my $query = <<"END";

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Protein 
From 1 
Condition CE*
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class Protein 
From 1 
Tag Pep_homol 
Condition MSP*
 
Colonne 3 
Width 12 
Mandatory 
Visible 
Class 
Class Method 
Right_of 2 
Tag  HERE  
Condition mass-spec
 
Colonne 4 
Width 12 
Mandatory 
Visible 
Float 
Right_of 3 
Tag  HERE  
 
Colonne 5 
Width 12 
Mandatory 
Visible 
Integer 
Right_of 4 
Tag  HERE  
 
Colonne 6 
Width 12 
Mandatory 
Visible 
Integer 
Right_of 5 
Tag  HERE  
 
END

  print $tm_def_fh "$query\n";
  close($tm_def_fh);

  return $tm_def_file;

}


sub get_uniprot_id {
    my ($wormpep, $db) = @_;

    my $obj = $db->fetch(Protein => $wormpep);
    return unless $obj->Database;
    my @databases = $obj->at('DB_info.Database');
    for my $database (@databases) {
	next unless $database->name eq 'UniProt' and $database->right->name eq 'UniProtAcc';
	return $database->right->right->name;
    }

    return;
}
