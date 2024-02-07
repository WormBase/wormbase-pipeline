#!/usr/bin/env perl
use warnings;
use strict;

use Getopt::Long;

use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Compara::DBSQL::DBAdaptor;

my $mlss_id;
my $compara_hal_dir;

GetOptions(
	    "compara_hal_dir=s" => \$compara_hal_dir,
	        "mlss_id=i"         => \$mlss_id
	) or die("Error in command line arguments.\n");

# Check if required options are provided
unless (defined $compara_hal_dir && defined $mlss_id) {
	die("Usage: $0 --compara_hal_dir=<HAL_DIR> --mlss_id=<MLSS_ID>\n");
}

$ENV{COMPARA_HAL_DIR} = $compara_hal_dir;

my $PARASITE_STAGING_MYSQL = $ENV{PARASITE_STAGING_MYSQL};

# my ( $species, $chr, $start, $end ) = ('caenorhabditis_angaria_prjna51225', 'Cang_2012_03_13_00001', 231967, 232967);
my ( $species, $chr, $start, $end ) = ('caenorhabditis_remanei_prjna53967', 'Crem_Contig0', 2412208, 2413208);
# my ( $species, $chr, $start, $end ) = ('schmidtea_mediterranea_s2f19h1prjna885486', '1_h1', 2428871, 2440280); #1
# my ( $species, $chr, $start, $end ) = ('trichinella_zimbabwensis_prjna257433', 'scaffold1s', 46637, 62176); #2
# my ( $species, $chr, $start, $end ) = ('ptycholaimellus_gst110_prjna953805','QYR16.scaffold0007',89018, 137765); #3
# my ( $species, $chr, $start, $end ) = ('syphacia_muris_prjeb524','SMUV_scaffold0000001',250977,256553); #4
# my ( $species, $chr, $start, $end ) = ('Rhabditophanes_kr3021_prjeb1297','RSKR_scaffold0000001',4104838,4106098); #6
# my ( $species, $chr, $start, $end ) = ('enoplolaimus_lenunculus_prjna953805','QYR23.scaffold0007',776043,791852); #5
# my ( $species, $chr, $start, $end ) = ('Micoletzkya_japonica_prjeb27334','scaffold48',717227,722613); #9
# my ( $species, $chr, $start, $end ) = ('Mesorhabditis_belari_prjeb61636','Mbe_Ptig_118',5091943,5103000); #7
# my ( $species, $chr, $start, $end ) = ('heterorhabditis_bacteriophora_prjna13977','Scaffold1250',29645,32693); #8



Bio::EnsEMBL::Registry->load_registry_from_url('mysql://ensro@mysql-ps-staging-1.ebi.ac.uk:4451/111');

my $compara_url = 'mysql://ensro@mysql-ps-staging-1.ebi.ac.uk:4451/ensembl_compara_parasite_19_111';
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
