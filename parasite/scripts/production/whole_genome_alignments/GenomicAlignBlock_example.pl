#!/usr/bin/env perl
use warnings;
use strict;

use Getopt::Long;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

$ENV{COMPARA_HAL_DIR} = '/hps/nobackup/flicek/wormbase/parasite/cactus_alignments/WBPS18/Anc028';

my $mlss_id = 4822;
my $PARASITE_STAGING_MYSQL = $ENV{PARASITE_STAGING_MYSQL};

my ( $species, $chr, $start, $end ) = ('caenorhabditis_angaria_prjna51225', 'Cang_2012_03_13_00001', 231967, 232967);

Bio::EnsEMBL::Registry->load_registry_from_url("mysql://ensro\@mysql-ps-staging-2.ebi.ac.uk:4467/108");

my $compara_url = "mysql://ensro\@mysql-ps-staging-2.ebi.ac.uk:4467/ensembl_compara_parasite_18_108";
my $compara_dba = Bio::EnsEMBL::Compara::DBSQL::DBAdaptor->go_figure_compara_dba($compara_url);
my $cactus_mlss = $compara_dba->get_MethodLinkSpeciesSetAdaptor->fetch_by_dbID($mlss_id);
print "cactus_mlss = " . $cactus_mlss->dbID . "\n";

print "Fetching $chr:$start-$end from $species (" . ($end-$start+1) . " bp)\n";

my $sliceAdaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'Core', 'Slice');
my $slice = $sliceAdaptor->fetch_by_region('supercontig', $chr, $start, $end);
my $slice_gabs = $compara_dba->get_GenomicAlignBlockAdaptor->fetch_all_by_MethodLinkSpeciesSet_Slice( $cactus_mlss, $slice );

my $c = 0;
foreach my $gab ( @$slice_gabs ) {
    print "$c : " . $gab->toString . "\n";
    $c++;
}
print "Got $c blocks! Yay!\n\n";
