#!/usr/local/bin/perl5.8.0 -w
#
# autoace_builder.pl
#
# based on the original autoace_minder.pl
#
# Usage : autoace_builder.pl [-options]
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2006-01-10 14:00:02 $

my $script_dir = $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'};

use strict;
use Wormbase;
use Getopt::Long;
use File::Copy;
use Coords_converter;
use Log_files;

my ( $debug, $test, $database );
my ( $initiate, $prepare_databases, $acefile, $build, $first_dumps );
my ( $make_wormpep, $finish_wormpep );
my ( $run_blat,     $finish_blat );
my ( $gff_dump,     $processGFF, $gff_split );
my $gene_span;
my ( $load, $tsuser, $map, $transcripts );

GetOptions(
    'debug:s'        => \$debug,
    'test'           => \$test,
    'database:s'     => \$database,
    'initiate:s'     => \$initiate,
    'prepare'        => \$prepare_databases,
    'acefiles'       => \$acefile,
    'build'          => \$build,
    'first_dumps'    => \$first_dumps,
    'make_wormpep'   => \$make_wormpep,
    'finish_wormpep' => \$finish_wormpep,
    'gff_dump:s'     => \$gff_dump,
    'processGFF:s'   => \$processGFF,
    'gff_split'      => \$gff_split,
    'gene_span'      => \$gene_span,
    'load=s'         => \$load,
    'run_blat'       => \$run_blat,
    'finish_blat'    => \$finish_blat,
    'tsuser=s'       => \$tsuser,
    'map'            => \$map,
    'transcripts'    => \$transcripts
);

my $wormbase = Wormbase->new(
    -test    => $test,
    -debug   => $debug,
    -version => $initiate
);

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$wormbase->run_script( "initiate_build.pl",                 $log ) if defined($initiate);
$wormbase->run_script( 'prepare_primary_databases.pl',      $log ) if $prepare_databases;
$wormbase->run_script( 'make_acefiles.pl',                  $log ) if $acefile;
$wormbase->run_script( 'make_autoace.pl',                   $log ) if $build;
#//--------------------------- batch job submission -------------------------//
$wormbase->run_script( "build_dumpGFF.pl -stage $gff_dump", $log ) if $gff_dump;
$wormbase->run_script( "processGFF.pl -$processGFF", $log ) if $processGFF;    #clone_acc
&first_dumps if $first_dumps;


$wormbase->run_script( 'make_wormpep.pl -initial', $log ) if $make_wormpep;

#########   BLAT  ############
$wormbase->run_script( 'BLAT_controller.pl -mask -dump -run', $log ) if $run_blat;

#//--------------------------- batch job submission -------------------------//

$wormbase->run_script( 'BLAT_controller.pl -virtual -process -postprocess -ace -load', $log ) if $finish_blat;

# mapping part
&map_features if $map;

######## WBGene spans ########
$wormbase->run_script( 'WBGene_span.pl -prepare', $log ) if $gene_span;

if ($load) {
    $log->("loading $load to $database\n");
    $log->("\ttsuser = $tsuser\n\n");
    $wormbase->load_to_database( $database, $load, $tsuser ) if ( -e $load );
}

$log->mail;
exit(0);

############################
#       SUBROUTINES        #
############################

sub first_dumps {
    $wormbase->run_script( "chromosome_dump.pl --dna --composition", $log );
    $wormbase->run_script( "make_agp_file.pl",                       $log );
    $wormbase->run_script( "agp2dna.pl",                             $log ); #dependant on processGFF producing clone_acc files.

    my $agp_errors = 0;

    my @chrom = qw( I II III IV V X);
    foreach my $chrom (@chrom) {
        open( AGP, ">" . $wormbase->autoace . "/yellow_brick_road/CHROMOSOME_${chrom}.agp_seq.log" )
          or die "Couldn't open agp file : $!";
        while (<AGP>) {
            $agp_errors++ if (/ERROR/);
        }
        close(AGP);
    }
    $log->write_to("ERRORS ( $agp_errors ) in agp file\n");
}

sub map_features {

    # features
    $wormbase->run_script( 'map_features.pl -all -build', $log );

    # PCR products
    $wormbase->run_script( 'map_PCR_products.pl', $log );

    #Oligo_sets
    $wormbase->run_script( 'map_Oligo_set.pl', $log );

    # RNAi experiments
    $wormbase->run_script( 'map_RNAi.pl -load', $log );

    # alleles
    $wormbase->run_script( 'map_Alleles.pl', $log );

    # Y2H objects
    $wormbase->run_script( 'map_Y2H.pl -load', $log );

    # microarray connections
    $wormbase->run_script( 'map_microarray.pl -load', $log );

    # TSL features
    $wormbase->run_script( 'map_feature2gene.pl -load', $log );
}

#__ end map_features __#
