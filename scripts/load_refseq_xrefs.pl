#!/usr/bin/env perl
#
# load_refseq_xrefs
# 
# A script to load xrefs to RefSeq. The source file is expected to live here:
#
# BUILD_DATA/MISC_STATIC/XREFS/RefSeq/refseq_xrefs.txt
# 
# This file is supplied by RefSeq each time they update (approx once per year)
# and corresponds to a specific WormBase freeze. Since WormBase will have moved
# on since then, we take care not to add xrefs to objects that are no longer live
#
# Last edited by: $Author: klh $
# Last edited on: $Date: 2013-10-14 09:54:27 $
#
#==================

use strict;
use lib $ENV{CVS_DIR};
use Wormbase;
use Getopt::Long;
use Log_files;
use Storable;

##############################
# command-line options       #
##############################

my ($test,$debug,$wormbase,$store,$noload,$database,$xreffile,$acefile);

GetOptions (
  "debug:s"     => \$debug,
  "store"       => \$store,
  "test"        => \$test, 
  "database:s"  => \$database,
  "acefile=s"   => \$acefile,
  "xreffile=s"  => \$xreffile,
  "noload"      => \$noload,
	   );

if ($store) {
  $wormbase = retrieve($store) or croak ("Can't restore wormbase from $store\n");
}
else {
  $wormbase = Wormbase->new( -debug => $debug,
                             -test => $test,
                           );
}

my $log = Log_files->make_build_log($wormbase);

$database = $wormbase->autoace if not defined $database;
$xreffile = $wormbase->misc_static . "/XREFS/RefSeq/refseq_xrefs.txt" if not defined $xreffile;
$acefile = $wormbase->acefiles . "/refseq_xrefs.ace" if not defined $acefile;

my (%locus_tag, %xrefs);

$log->write_to("Reading gene information...\n");

my $def = &get_tm_def();
my $tm_query = $wormbase->table_maker_query($database,$def);
while(<$tm_query>) {
  chomp;
  next if (/>/ or /\/\//);
  s/\"//g;

  my($gene,$locus_tag) = split(/\t/, $_);
  next unless $locus_tag;

  $locus_tag{$locus_tag} = $gene;
}
unlink $def;

$log->write_to("Reading xref file...\n");

open(my $fh, $xreffile) or $log->log_and_die("Could not open $xreffile for reading\n");
while(<$fh>) {
  /^\#/ and next;
  my @l = split(/\t/, $_);
  my ($entrez_gene, $locus_tag, $refseq_trans, $refseq_prot) = @l[0,1,2,4];

  if (exists $locus_tag{$locus_tag}) {
    my $g = $locus_tag{$locus_tag};
    $xrefs{$g}->{NCBI}->{gene}->{$entrez_gene} = 1;
    $xrefs{$g}->{RefSeq}->{mRNA}->{$refseq_trans} = 1 if $refseq_trans ne "-";
    $xrefs{$g}->{RefSeq}->{protein}->{$refseq_prot} = 1 if $refseq_prot ne "-";

  } else {
    $log->write_to("Ignoring $locus_tag - does not seem to be associated with live gene\n");
  }
}

open(my $outfh, ">$acefile") or $log->log_and_die("Could not open $acefile for writing\n");
foreach my $gene (sort keys %xrefs) {
  print $outfh "\nGene : \"$gene\"\n";
  foreach my $db (sort keys %{$xrefs{$gene}}) {
    foreach my $field (sort keys %{$xrefs{$gene}->{$db}}) {
      foreach my $acc (keys %{$xrefs{$gene}->{$db}->{$field}}) {
        print $outfh "Database \"$db\" \"$field\" \"$acc\"\n";
      }
    }
  }
}
close($outfh) or $log->log_and_die("Could not close $acefile\n");

if (not $noload) {
  $log->write_to("Loading results...\n");
  $wormbase->load_to_database($database,$acefile,'load_refseq_xrefs', $log);
}

$log->mail();
exit(0);

##########################################
sub get_tm_def {

  my $tm_def = "/tmp/load_refseq_xrefs.$$.tm_def";
  open(my $outfh, ">$tm_def") 
      or $log->log_and_die("Could not open $tm_def for writing\n");

  my $query =  <<"EOF";
Sortcolumn 1

Colonne 1 
Width 12 
Optional 
Visible 
Class 
Class Gene 
From 1 
Condition Live
 
Colonne 2 
Width 12 
Mandatory 
Hidden 
Class 
Class Database 
From 1 
Tag Database 
Condition "NDB"
 
Colonne 3 
Width 12 
Optional 
Hidden 
Class 
Class Database_field 
Right_of 2 
Tag  HERE  
Condition "locus_tag"
 
Colonne 4 
Width 12 
Optional 
Visible 
Class 
Class Text 
Right_of 3 
Tag  HERE  

EOF

  print $outfh $query;
  close($outfh);

  return $tm_def;
}


__END__
