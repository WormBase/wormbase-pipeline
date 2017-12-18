#!/usr/bin/env perl
#
# promote_xrefs_to_gene
# 
# Inspects xrefs for CDSs and proteins, and for each type, if there
# is only one for that type, promote it up to the gene
#
# Version: $Version: $
# Last updated by: $Author: klh $
# Last updated on: $Date: 2013-10-03 14:51:23 $

use strict;
use warnings;
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Ace;

###############
# variables   #
###############

my (
  $store,
  $test,
  $debug,
  $species,
  $noload,
  $acefile,
  $ace_fh,
  $wormbase,
    );

GetOptions(
  "debug=s"       => \$debug,
  'store=s'       => \$store,
  'species=s'     => \$species,
  "test"          => \$test,
  "noload"        => \$noload,
  "acefile=s"     => \$acefile,
    );

#################################################
# config 
#################################################
my $to_promote = {
  #UniProt => {
  #  UniProtAcc => [['UniProt', 'UniProtAcc']],
  #},
  SwissProt => {
    UniProtAcc => [['SwissProt', 'UniProtAcc']],
  },
  TrEMBL => {
    UniProtAcc => [['TrEMBL', 'UniProtAcc']],
  },
  UniProt_GCRP => {
    UniProtAcc => [['UniProt_GCRP', 'UniProtAcc']],
  },
  KEGG => {
    KEGG_id => [['KEGG', 'KEGG_id'],
                ['NemaPath', 'KEGG_id']],
  },
  TREEFAM => {
    TREEFAM_ID => [['TREEFAM', 'TREEFAM_ID']],
  },
};
  

############################
# recreate configuration   #
############################
if ($store) { 
  $wormbase = Storable::retrieve($store) or croak("cant restore wormbase from $store\n"); 
}
else { 
  $wormbase = Wormbase->new(-debug => $debug, 
                            -test => $test,
                            -organism => $species,
      );
}

my $log = Log_files->make_build_log($wormbase);

##############
# Paths etc. #
##############

$acefile = $wormbase->acefiles . "/protein_xref_promotion.ace" if not defined $acefile;

my %cds2gene = $wormbase->FetchData('cds2wbgene_id');
my %gene_info;

my $tmpdef = &write_TM_def;
my $tm_query = $wormbase->table_maker_query($wormbase->autoace,$tmpdef);
while(<$tm_query>) {
  s/\"//g;
  next if (/acedb/ or /\/\// or /^\s*$/);
  chomp;
  my ($wormpep_id, $cds_id, $db_name, $db_field, $db_val) = split(/\s+/, $_);

  next if not exists $cds2gene{$cds_id} or not defined $db_name;
  my $gene_id = $cds2gene{$cds_id};

  next if not exists $to_promote->{$db_name};
  next if not exists $to_promote->{$db_name}->{$db_field};

  $gene_info{$gene_id}->{$db_name}->{$db_field}->{$db_val} = 1;
}
unlink $tmpdef;

open($ace_fh, ">$acefile") or $log->log_and_die("Could not open $acefile for writing\n");
foreach my $g (sort keys %gene_info) {
  print $ace_fh "\nGene : $g\n";
  foreach my $db (sort keys %{$gene_info{$g}}) {
    foreach my $dbf (sort keys %{$gene_info{$g}->{$db}}) {
      my @vals = keys %{$gene_info{$g}->{$db}->{$dbf}};
      foreach my $val (@vals) {
        foreach my $pair (@{$to_promote->{$db}->{$dbf}}) {
          my ($to_db, $to_dbf) = @$pair;
          print $ace_fh "Database $to_db $to_dbf $val\n";
        }
      }
    }
  }
}
close($ace_fh) or $log->log_and_die("Could not close output file $ace_fh\n");

unless ($noload) {
  $wormbase->load_to_database($wormbase->autoace, $acefile, "promote_prot_xrefs_to_gene", $log);
}

$log->mail();
exit(0);

###############################################
sub write_TM_def {
  my $def = '/tmp/promote_xrefs_to_gene';
  open my $tmpdefh,">$def" or $log->log_and_die("cant write $def: $!\n");
  my $txt = <<END;

Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Protein 
From 1 
 
Colonne 2 
Width 12 
Mandatory 
Visible 
Class 
Class CDS 
From 1 
Tag Corresponding_CDS 
 
Colonne 3 
Width 12 
Optional 
Visible 
Class 
Class Database 
From 2 
Tag Database 
 
Colonne 4 
Width 12 
Optional 
Visible 
Class 
Class Database_field 
Right_of 3 
Tag  HERE  
 
Colonne 5 
Width 12 
Optional 
Visible 
Class 
Class Text 
Right_of 4 
Tag  HERE  

END

  print $tmpdefh $txt;
  close($tmpdefh);
  return $def;
}
############################################

__END__
