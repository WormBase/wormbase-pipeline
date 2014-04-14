#!/usr/bin/env perl

use strict;
use Getopt::Long;

use FindBin qw($Bin);
use lib "$Bin/../lib";

use WormBase2Ensembl;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


my ($species, $debug, $verbose,
    $dbhost, $dbname, $dbuser, $dbpass, $dbport,
    $gff_file, $analysis_logic, @gff_source,
    $ignore_phases,
    );

if (not &GetOptions(
      'species=s'         => \$species,
      'host=s'            => \$dbhost,
      'user=s'            => \$dbuser, 
      'pass=s'            => \$dbpass,
      'port=s'            => \$dbport,
      'dbname=s'          => \$dbname,
      'gff=s'             => \$gff_file,
      'analysis=s'        => \$analysis_logic,
      'debug'             => \$debug,
      'verbose'           => \$verbose,
      'gffsource=s@'      => \@gff_source,
      'ignorephases'      => \$ignore_phases)) {
  die "Could not process options. Quitting\n";
}
  

$analysis_logic = "wormbase" if not defined $analysis_logic;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
    -host   => $dbhost,
    -user   => $dbuser,
    -dbname => $dbname,
    -pass   => $dbpass,
    -port   => $dbport,
);


my $wb2ens = WormBase2Ensembl->new(
  -species         => $species,
  -debug           => $debug,
  -verbose         => $verbose,
  -ignoregffphases => $ignore_phases,
  -dbh             => $db);

my $analysis = $db->get_AnalysisAdaptor->fetch_by_logic_name($analysis_logic);

if (not defined $analysis) {
  $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => $analysis_logic,
                                          -gff_source => "WormBase_imported");
}


my $source_hash;
if (@gff_source) {
  $source_hash = {};
  map { $source_hash->{$_} = 1 } @gff_source;
}

my $genes = $wb2ens->parse_genes_gff3( $gff_file, $analysis, $analysis, $analysis, $source_hash );

$verbose and printf STDERR "Parsed %d genes from GFF3. Writing genes to database...\n", scalar(@$genes);      
$wb2ens->write_genes( $genes, 1 );
    
$verbose and print STDERR "Checking translations...\n";
my @genes = @{$db->get_GeneAdaptor->fetch_all_by_biotype('protein_coding')};
&check_translations(\@genes);

exit 0;


##############################################################
sub check_translations {
  my ($genes) = @_;

  foreach my $gene (@$genes) {
    my @non_translate;
    foreach my $tran (@{$gene->get_all_Transcripts}) {
      if ($tran->biotype eq 'protein_coding' and not $wb2ens->translation_check($tran)) {
        push @non_translate, $tran;
      }
    }
    
    if (@non_translate) {
      print "Transcripts from " . $gene->stable_id . " do not translate\n";
    }
  }
}
