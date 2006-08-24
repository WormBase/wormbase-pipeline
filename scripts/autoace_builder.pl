#!/nfs/disk100/wormpub/bin/perl
#
# autoace_builder.pl
#
# based on the original autoace_minder.pl
#
# Usage : autoace_builder.pl [-options]
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2006-08-24 17:02:48 $

my $script_dir = $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'};

use strict;
use Wormbase;
use Getopt::Long;
use File::Copy;
use Coords_converter;
use Log_files;
use Storable;

my ( $debug, $test, $database );
my ( $initiate, $prepare_databases, $acefile, $build, $first_dumps );
my ( $make_wormpep, $finish_wormpep );
my ( $run_blat,     $finish_blat );
my ( $gff_dump,     $processGFF, $gff_split );
my $gene_span;
my ( $load, $tsuser, $map_features, $remap_misc_dynamic, $map, $transcripts, $intergenic, $data_sets, $nem_contigs);
my ( $GO_term, $rna , $dbcomp, $confirm, $operon ,$repeats, $remarks, $names, $treefam, $cluster);
my ( $utr, $agp, $gff_munge, $extras ,$interpolate, $check);
my ( $data_check, $buildrelease, $public,$finish_build, $release);


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
	   'remap_misc_dynamic' => \$remap_misc_dynamic,
	   'map_features'   => \$map_features,
	   'transcripts'    => \$transcripts,
	   'intergenic'     => \$intergenic,
	   'nem_contig'     => \$nem_contigs,
	   'data_sets'      => \$data_sets,
	   'go_term'        => \$GO_term,
	   'rna'            => \$rna,
	   'dbcomp'         => \$dbcomp,
	   'confirm'        => \$confirm,
	   'operon'         => \$operon,
	   'repeats'        => \$repeats,
	   'remarks'        => \$remarks,
	   'names'          => \$names,
	   'treefam'        => \$treefam,
	   'cluster'        => \$cluster,
	   'utr'            => \$utr,
	   'interpolation'  => \$interpolate,
	   'agp'            => \$agp,
	   'gff_munge'      => \$gff_munge,
	   'extras'         => \$extras,
	   'buildrelease'   => \$buildrelease,
	   'public'         => \$public,
	   'finish_build'   => \$finish_build,
	   'release'        => \$release,
	   'check'    	    => \$check,
	   'data_check'     => \$data_check
	  );

my $wormbase = Wormbase->new(
    -test    => $test,
    -debug   => $debug,
    -version => $initiate
);

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$wormbase->run_script( "initiate_build.pl -version $initiate",$log ) if defined($initiate);
$wormbase->run_script( 'prepare_primary_databases.pl',      $log ) if $prepare_databases;
$wormbase->run_script( 'make_acefiles.pl',                  $log ) if $acefile;
$wormbase->run_script( 'make_autoace.pl',                   $log ) if $build;

#//--------------------------- batch job submission -------------------------//
$wormbase->run_script( "build_dumpGFF.pl -stage $gff_dump", $log ) if $gff_dump;      #init

$wormbase->run_script( "processGFF.pl -$processGFF",        $log ) if $processGFF;    #clone_acc
&first_dumps                                                       if $first_dumps;   # dependant on clone_acc for agp
$wormbase->run_script( 'make_wormpep.pl -initial',          $log ) if $make_wormpep;
$wormbase->run_script( 'map_features.pl -all',              $log ) if $map_features;


#########   BLAT  ############
$wormbase->run_script( 'BLAT_controller.pl -mask -dump -run', $log ) if $run_blat;
#//--------------------------- batch job submission -------------------------//
$wormbase->run_script( 'BLAT_controller.pl -virtual -process -postprocess -ace -load', $log ) if $finish_blat;
# $build_dumpGFF.pl; (blat) is run chronologically here but previous call will operate

$wormbase->run_script( 'batch_transcript_build.pl', $log) if $transcripts;
#requires GFF dump of transcripts (done within script if all goes well)

$wormbase->run_script( 'WBGene_span.pl'                   , $log ) if $gene_span;
&make_UTR($log)                                                    if $utr;

$wormbase->run_script( 'find_intergenic.pl'               , $log ) if $intergenic;
$wormbase->run_script( 'inherit_GO_terms.pl -phenotype'   , $log ) if $GO_term;

##  Horrid Geneace related stuff  ##########
#make_pseudo_map_positions.pl -load
#get_interpolated_gmap.pl
#update_inferred_multi_pt.pl -load

####### mapping part ##########
&map_features                                                            if $map;

&remap_misc_dynamic                                                      if $remap_misc_dynamic;

&get_repeats                                                             if $repeats; # loaded with homols
#must have farm complete by this point.
$wormbase->run_script( 'load_data_sets.pl -homol -briggsae -misc', $log) if $data_sets;
# $build_dumpGFF.pl; (homol) is run chronologically here but previous call will operate
$wormbase->run_script( 'make_wormrna.pl'                         , $log) if $rna;
$wormbase->run_script( 'confirm_genes.pl -load'                  , $log) if $confirm;
$wormbase->run_script( 'map_operons.pl'                          , $log) if $operon;
$wormbase->run_script( 'make_wormpep.pl -final'                  , $log) if $finish_wormpep;
$wormbase->run_script( 'write_DB_remark.pl'                      , $log) if $remarks;
$wormbase->run_script( 'molecular_names_for_genes.pl'            , $log) if $names;
$wormbase->run_script( 'get_treefam.pl'                          , $log) if $treefam;
$wormbase->run_script( 'cluster_gene_connection.pl'              , $log) if $cluster;

# $build_dumpGFF.pl; (final) is run chronologically here but previous call will operate
# $wormbase->run_script( "processGFF.pl -$processGFF",        $log ) if $processGFF;    #nematode - to add species to nematode BLATs

$wormbase->run_script( "interpolation_manager.pl"                , $log) if $interpolate;
$wormbase->run_script( "make_agp_file.pl"                        , $log) if $agp;
$wormbase->run_script( "landmark_genes2gff.pl"                   , $log) if $gff_munge;
$wormbase->run_script( "GFFmunger.pl -all"                       , $log) if $gff_munge;
&make_extras                                                             if $extras;

#run some checks
$wormbase->run_script( "post_build_checks.pl -a"                 , $log) if $check;
$wormbase->run_script( "data_checks.pl -ace -gff"                , $log) if $data_check;
$wormbase->run_script( "dbcomp.pl"                               , $log) if $data_check;
$wormbase->run_script( "build_release_files.pl"                  , $log) if $buildrelease;
&public_sites                                                            if $public;
$wormbase->run_script( "distribute_letter.pl"                    , $log) if $release;

$wormbase->run_script("finish_build.pl"                          , $log) if $finish_build;
$wormbase->run_script("update_gffdb.csh"                         , $log) if $finish_build;

if ($load) {
    $log->write_to("loading $load to ".$wormbase->autoace."\n");
    $log->write_to("\ttsuser = $tsuser\n\n");
    $wormbase->load_to_database( $wormbase->autoace, $load, $tsuser ,$log) if ( -e $load );
}

$log->mail;

exit(0);


############################
#       SUBROUTINES        #
############################

sub first_dumps {
    $wormbase->run_script( "chromosome_dump.pl --dna --composition", $log );

    my $version = $wormbase->get_wormbase_version;
    $wormbase->run_script( "inspect-old-releases.pl -version $version -database1 ".$wormbase->database('current')." -database2 ".$wormbase->autoace, $log );

    $wormbase->run_script( "make_agp_file.pl",                       $log );
    $wormbase->run_script( "agp2dna.pl",                             $log ); #dependant on processGFF producing clone_acc files.

    my $agp_errors = 0;

    my @chrom = qw( I II III IV V X);
    foreach my $chrom (@chrom) {
        open( AGP, "<" . $wormbase->autoace . "/yellow_brick_road/CHROMOSOME_${chrom}.agp_seq.log" )
          or die "Couldn't open agp file : $!";
        while (<AGP>) {
            $agp_errors++ if (/ERROR/);
        }
        close(AGP);
    }
    
	$log->write_to("ERRORS ( $agp_errors ) in agp file\n");
}

sub map_features {

    # PCR products  - requires UTR GFF files
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

    # writes tables listing microarrays to genes
    $wormbase->run_script( 'make_oligo_set_mapping_table.pl -all', $log );
}

#__ end map_features __#

sub remap_misc_dynamic {

  my $release = $wormbase->get_wormbase_version;
  my $previous_release = $release - 1;

  # remap twinscan
  my $twinscan = $wormbase->misc_dynamic."/misc_twinscan.ace";
  my $backup_twinscan = "$twinscan.$previous_release";

  if (-e $backup_twinscan) {$log->log_and_die("$backup_twinscan already exists - please move it to be $twinscan before running this again\n");}
  $wormbase->run_command("cp $twinscan $backup_twinscan", $log);
  $wormbase->run_script( "remap_twinscan_between_releases.pl -release1 $previous_release -release2 $release -twinscan $backup_twinscan -out $twinscan", $log);

  # remap genefinder
  my $genefinder = $wormbase->misc_dynamic."/misc_genefinder.ace";
  my $backup_genefinder = "$genefinder.$previous_release";

  if (-e $backup_genefinder) {$log->log_and_die("$backup_genefinder already exists - please move it to be $genefinder before running this again\n");}
  $wormbase->run_command("cp $genefinder $backup_genefinder", $log);
  $wormbase->run_script( "remap_genefinder_between_releases.pl -input $backup_genefinder -out $genefinder", $log);

  # remap fosmids
  my $fosmids = $wormbase->misc_dynamic."/fosmids.ace";
  my $backup_fosmids = "$fosmids.$previous_release";

  if (-e $backup_fosmids) {$log->log_and_die("$backup_fosmids already exists - please move it to be $fosmids before running this again\n");}
  $wormbase->run_command("cp $fosmids $backup_fosmids", $log);
  $wormbase->run_script( "remap_fosmids_between_releases.pl -input $backup_fosmids -out $fosmids", $log);
   
  # remap and copy over the SUPPLEMENTARY_GFF dir from BUILD_DATA
  my $sup_dir = $wormbase->build_data."/SUPPLEMENTARY_GFF";
  my $release = $wormbase->version;
  my $old_release = $release - 1;
  opendir(DIR,$sup_dir) or $log->log_and_die("cant open $sup_dir: $!\n");
  while ( my $file = readdir( DIR ) ) {
  		next unless( $file =~ /gff$/ );
	  	$wormbase->run_script("remap_gff_between_releases.pl -gff $sup_dir/$file -output /tmp/remap_GFF -release1 $old_release -release2 $release", $log);
	  	$wormbase->run_command("mv /tmp/remap_GFF $sup_dir/$file", $log);
  	}
  	closedir DIR;
  	$wormbase->run_command("cp -R ".$wormbase->build_data."/SUPPLEMENTARY_GFF ".$wormbase->chromosomes."/");
   
  # remap waba (takes a long time ~24 hours)
#  my $waba = $wormbase->misc_dynamic."/waba.ace";
#  my $backup_waba = "$waba.$previous_release";
#
#  if (-e $backup_waba) {$log->log_and_die("$backup_waba already exists - please move it to be $waba before running this again\n");}
#  $wormbase->run_command("cp $waba $backup_waba", $log);
#  $wormbase->run_script( "remap_waba_between_releases.pl -input $backup_waba -out $waba", $log);


}

#__ end remap_misc_dynamic __#

sub make_UTR {
  my ($log)=@_;
  foreach (qw( I II III IV V X ) ) {
	  # crude ... should be beautified
	  my $store = $wormbase->autoace . "/wormbase.store";
	  $wormbase->run_command("bsub -J make_UTRs -o /dev/null make_UTR_GFF.pl -chromosome $_ -store $store",$log)
  }
}


sub get_repeats {
  #repeatmasked chromosomes
  my $wormpipe= glob("~wormpipe");
  my $release = $wormbase->get_wormbase_version;
  my $agp = $wormpipe."/Elegans/WS$release.agp";
  $wormbase->run_script("get_repeatmasked_chroms.pl -agp $agp", $log);

  #inverted
  $wormbase->run_script("run_inverted.pl -all" , $log);
}


sub make_extras {
  my $version = $wormbase->get_wormbase_version;
  $wormbase->run_script( "make_keysets.pl -all -history $version", $log);
  $wormbase->run_script( "genestats.pl" , $log);
}


sub public_sites {
  # gets everything on the to FTP and websites and prepares release letter ready for final edit and sending.
  $wormbase->run_script( "make_FTP_sites.pl -all", $log);
  $wormbase->run_script( "update_website.pl -all", $log);
  $wormbase->run_script( "release_letter.pl -l"  , $log);
}
