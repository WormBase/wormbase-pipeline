#!/usr/bin/env perl

# autoace_builder.pl
#
# based on the original autoace_minder.pl
#
# Usage : autoace_builder.pl [-options]
#
# Last edited by: $Author: klh $
# Last edited on: $Date: 2015-05-20 11:07:12 $

my $script_dir = $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'};

use strict;
use Wormbase;
use Getopt::Long;
use File::Copy;
use Coords_converter;
use Modules::Remap_Sequence_Change;
use Log_files;
use Storable;
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;

my ( $debug, $test, $database, $species);
my ( $initiate, $prepare_databases, $acefile, $build, $build_check, $assembly, $first_dumps );
my ( $make_wormpep, $finish_wormpep );
my ( $prep_blat, $run_blat,     $finish_blat );
my ( $gff_dump,     $processGFF, $gff_split );
my $gene_span;
my ( $load, $big_load, $tsuser );
my ($map_features, $remap_misc_dynamic, $map, $map_alleles, $transcripts, $cdna_files, $misc_data_sets, $homol_data_sets, $nem_contigs);
my ( $GO_term, $rna , $dbcomp, $confirm, $operon ,$repeats, $treefam, $ncbi_xrefs, $load_interpro, $RRID, $omim);
my ( $utr, $agp, $gff_munge, $gff3_munge, $extras , $ontologies, $interpolate, $check, $enaseqxrefs, $enagenexrefs, $enaprotxrefs, $xrefs);
my ( $data_check, $buildrelease, $public,$finish_build, $gffdb, $autoace, $user, $kegg, $prepare_gff_munge, $post_merge, $gtf);


GetOptions(
	   'debug:s'        => \$debug,
	   'test'           => \$test,
	   'database:s'     => \$database,
	   'initiate:s'     => \$initiate,
	   'prepare'        => \$prepare_databases,
	   'acefiles'       => \$acefile,
	   'build'          => \$build,
           'buildcheck'     => \$build_check,
           'first_dumps'    => \$first_dumps,
	   'assembly'       => \$assembly,
	   'make_wormpep'   => \$make_wormpep,
           'loadinterpro'   => \$load_interpro,
	   'finish_wormpep' => \$finish_wormpep,
	   'gff_dump:s'     => \$gff_dump,
	   'processGFF:s'   => \$processGFF,
	   'gff_split'      => \$gff_split,
	   'gene_span'      => \$gene_span,
	   'load=s'         => \$load,
           'bigload=s'      => \$big_load,
	   'prep_blat'      => \$prep_blat,
	   'run_blat'       => \$run_blat,
	   'finish_blat'    => \$finish_blat,
	   'tsuser=s'       => \$tsuser,
	   'map'            => \$map,
	   'remap_misc_dynamic' => \$remap_misc_dynamic,
	   'map_features'   => \$map_features,
           'map_alleles'    => \$map_alleles,
	   'transcripts'    => \$transcripts,
	   'cdnafiles'      => \$cdna_files,
	   'nem_contig'     => \$nem_contigs,
	   'misc_data_sets' => \$misc_data_sets,
	   'homol_data_sets'=> \$homol_data_sets,
           'sequencexrefs'  => \$enaseqxrefs,
           'genexrefs'      => \$enagenexrefs,
           'proteinxrefs'   => \$enaprotxrefs,
           'xrefs'          => \$xrefs,
	   'rrid'           => \$RRID,
	   'omimxref'       => \$omim,
	   'rna'            => \$rna,
	   'dbcomp'         => \$dbcomp,
	   'confirm'        => \$confirm,
	   'operon'         => \$operon,
	   'repeats'        => \$repeats,
	   'treefam'        => \$treefam,
           'ncbi_xrefs'     => \$ncbi_xrefs,
	   'utr'            => \$utr,
	   'interpolation'  => \$interpolate,
	   'agp'            => \$agp,
           'prepmunge'      => \$prepare_gff_munge,
	   'gff_munge'      => \$gff_munge,
	   'gff3_munge'     => \$gff3_munge,
           'gtf'            => \$gtf,
	   'extras'         => \$extras,
	   'ontologies'     => \$ontologies,
	   'buildrelease'   => \$buildrelease,
	   'public'         => \$public,
	   'finish_build'   => \$finish_build,
	   'gffdb'          => \$gffdb,
	   'autoace'        => \$autoace,
	   'check'    	    => \$check,
	   'data_check'     => \$data_check,
	   'species:s'      => \$species,
	   'user:s'         => \$user,
	   'kegg'           => \$kegg,
           'postmerge'      => \$post_merge,
	  )||die(@!);


my $wormbase = Wormbase->new(
    -test    => $test,
    -debug   => $debug,
    -version => $initiate,
    -organism=> $species
);

# establish log file.
my $log = Log_files->make_build_log($wormbase);

$wormbase->run_script( "initiate_build.pl -version $initiate",$log ) if defined($initiate);

$wormbase->run_script( 'prepare_primary_databases.pl',      $log ) if $prepare_databases;
$wormbase->run_script( "check_primary_database.pl -organism ${\$wormbase->species}", $log ) if $prepare_databases;

$wormbase->run_script( 'make_acefiles.pl',                  $log ) if $acefile;
$wormbase->run_script( 'make_autoace.pl',                   $log ) if $build;

if ($build_check) {
  # Check for missing curation by checking for Live genes that have a Sequence name but aren't connected to a current gene model.
  if ($wormbase->species eq 'elegans')  {
    $wormbase->run_script("check_predicted_genes.pl -database ".$wormbase->autoace." -build", $log);
    $wormbase->run_script("check_class.pl -camace -genace -caltech -misc_static -briggsae -stage init", $log);
    $wormbase->run_script("check_class.pl -incomplete -stage incomplete", $log);
  }
}


$wormbase->run_script( "chromosome_dump.pl --dna --composition", $log ) if $first_dumps;
$wormbase->run_script("update_Common_data.pl -clone2centre -clone2acc -clone2size -clone2dbid -genes2lab -worm_gene2cgc -worm_gene2geneID -worm_gene2class -est -est2feature -gene_id -clone2type -cds2cgc -rna2cgc -pseudo2cgc ", $log ) if $first_dumps;

$wormbase->run_script( "build_dumpGFF.pl -stage $gff_dump",                $log ) if $gff_dump;

$wormbase->run_script( "get_ena_submission_xrefs.pl -sequencexrefs", $log)  if $enaseqxrefs;
$wormbase->run_script( "get_ena_submission_xrefs.pl -genexrefs", $log)  if $enagenexrefs;

$wormbase->run_script( "processGFF.pl -$processGFF",                       $log ) if $processGFF;

&make_remap_data()    if $first_dumps;

&do_assembly_stuff() if $assembly;   # dependant on clone_acc for agp

$wormbase->run_script("GetPFAM_motifs.pl", $log) if $load_interpro;
$wormbase->run_script("GetInterPro_motifs.pl", $log) if $load_interpro;
$wormbase->run_script("check_class.pl -stage load_interpro -classes Motif", $log) if $load_interpro;

$wormbase->run_script( 'make_wormpep.pl -initial -all',                    $log ) if $make_wormpep;
$wormbase->run_script("check_class.pl -stage make_wormpep -classes Protein,Peptide", $log) if $make_wormpep;

&map_features_to_genome() if $map_features;

#########   BLAT  ############
$wormbase->run_script( 'BLAT_controller.pl -dump', $log ) if $prep_blat;
#//--------------------------- batch job submission -------------------------//
$wormbase->run_script( 'BLAT_controller.pl -run', $log )        if $run_blat;
#//--------------------------- batch job submission -------------------------//
$wormbase->run_script( 'BLAT_controller.pl -postprocess -process -intron -load', $log ) if $finish_blat;
$wormbase->run_script("check_class.pl -stage finish_blat -classes Homol_data", $log) if $finish_blat;
#//--------------------------- batch job submission -------------------------//
# $build_dumpGFF.pl; (blat) is run chronologically here but previous call will operate

$wormbase->run_script( 'batch_transcript_build.pl', $log) if $transcripts;
$wormbase->run_script("check_class.pl -stage transcripts -classes Transcript", $log) if $transcripts;
#requires GFF dump of transcripts (done within script if all goes well)

$wormbase->run_script( 'WBGene_span.pl'                   , $log ) if $gene_span;
&make_UTR($log)                                                    if $utr;

if ($cdna_files) {
  my $seqdir = $wormbase->sequences;
  $wormbase->run_script( 'find_intergenic.pl'               , $log );
  $wormbase->run_script( "fasta_dumper.pl -classmethod Transcript:Coding_transcript -output $seqdir/mRNA_transcripts.dna", $log);
  $wormbase->run_script( "fasta_dumper.pl -classmethod Pseudogene:Pseudogene -output $seqdir/pseudogenic_transcripts.dna", $log);
  if ($wormbase->species eq 'elegans') {
    my @options = "-classmethod CDS:Transposon_CDS:Transposon-mRNA -classmethod Pseudogene:Transposon_Pseudogene:Transposon-pseudogenic_transcript -classmethod Transcript:Transposon_ncRNA:Transposon-non-coding_transcript";
    $wormbase->run_script( "fasta_dumper.pl @options -output $seqdir/transposon_transcripts.dna", $log);
    $wormbase->run_script( "fasta_dumper.pl -classmethod Transposon:Transposon -output $seqdir/transposons.dna", $log);
  }
}

####### mapping part ##########
&map_features                                                                           if $map;

&map_alleles                                                                            if $map_alleles;

&remap_misc_dynamic                                                                     if $remap_misc_dynamic;

$wormbase->run_script("run_inverted.pl -all" , $log)                                    if $repeats;
$wormbase->run_script("check_class.pl -stage run_inverted -classes Feature_data", $log) if $repeats;


#must have farm complete by this point.
$wormbase->run_script( 'load_data_sets.pl -misc', $log) if $misc_data_sets;
$wormbase->run_script( 'load_data_sets.pl -homol', $log) if $homol_data_sets;
# $build_dumpGFF.pl; (homol) is run chronologically here but previous call will operate
$wormbase->run_script( 'make_wormrna.pl'                         , $log) if $rna;
if ($confirm) {
  $wormbase->run_script( 'confirm_genes.pl', $log);
  $wormbase->run_script( 'update_Common_data.pl -cds2status ', $log);
}    
$wormbase->run_script( 'map_operons.pl'                          , $log) if $operon;
if ($enaprotxrefs) {
  $wormbase->run_script( "get_ena_submission_xrefs.pl -proteinxrefs", $log);
  $wormbase->run_script( "propagate_cds_xrefs_to_protein.pl", $log );
  if ($wormbase->species eq 'elegans') {
    $wormbase->run_script( 'load_panther_xrefs.pl', $log);
    $wormbase->run_script( 'generate_dbxref_file.pl -nocodingtrans -ebiupload', $log);  
  }
}

$wormbase->run_script('make_wormpep.pl -all -final', $log) if $finish_wormpep;

$wormbase->run_script( 'get_treefam.pl'                          , $log) if $treefam;
$wormbase->run_script( 'KEGG.pl', $log )                                 if $kegg;

# $build_dumpGFF.pl; (final) is run chronologically here but previous call will operate
# $wormbase->run_script( "processGFF.pl -$processGFF",        $log ) if $processGFF;    #nematode - to add species to nematode BLATs

$wormbase->run_script( "interpolation_manager.pl -fix -pseudo", $log) if $interpolate;
$wormbase->run_script( "make_agp_file.pl"                        , $log) if $agp;

if ($ncbi_xrefs and $wormbase->species eq 'elegans') {
  $wormbase->run_script( 'get_GI.pl', $log );
  $wormbase->run_script( 'load_refseq_xrefs.pl', $log);
}

if ($prepare_gff_munge) {
  if ($wormbase->species eq 'elegans') {
    $wormbase->run_script( 'landmark_genes2gff.pl', $log);
    $wormbase->run_script( 'landmark_genes2gff.pl -gff3', $log);
    $wormbase->run_script( 'web_data/interpolate_gmap2pmap.pl', $log);
    $wormbase->run_script( 'web_data/interpolate_gmap2pmap.pl -gff3 -balancers', $log);
  }
  $wormbase->run_script( 'web_data/map_translated_features_to_genome.pl', $log);
  $wormbase->run_script( 'web_data/map_translated_features_to_genome.pl -gff3', $log);
}

#several GFF manipulation steps
if ($gff_munge or $gff3_munge) {
    #Used for finding the previous gff files in the staging area of the sanger file system
    #previously looked at the internal FTP dir but we no longer store the data at Sanger. 
  my $prev_gff_prefix = 
      join("/", 
           $wormbase->ftp_staging, 
           "releases", 
           "WS" . ($wormbase->version - 1),
           "species",
           $wormbase->full_name(-g_species => 1),
           $wormbase->ncbi_bioproject,
           join(".", 
                $wormbase->full_name(-g_species => 1),
                $wormbase->ncbi_bioproject, 
                "WS" . ($wormbase->version - 1),
                "annotations")
      );

  if ($gff_munge) {
    $wormbase->run_script( 'GFF_post_process/GFF_post_process.pl -all', $log); 
    my $new_gff = $wormbase->processed_GFF_file;
    if (not -e $new_gff) {
      $new_gff .= ".gz";
    }
    my $prev_gff =  $prev_gff_prefix . ".gff2.gz";
    $wormbase->run_script("generate_gff_report.pl -currentgff $new_gff -previousgff $prev_gff", $log);
  }

  if ($gff3_munge) {
    $wormbase->run_script( 'GFF_post_process/GFF_post_process.pl -all -gff3', $log); 
    my $new_gff = $wormbase->processed_GFF_file(1);
    if (not -e $new_gff) {
      $new_gff .= ".gz";
    }
    my $prev_gff =  $prev_gff_prefix . ".gff3.gz";
    $wormbase->run_script("generate_gff_report.pl -currentgff $new_gff -previousgff $prev_gff", $log);
  }
}

if ($gtf) {
  $wormbase->run_script("ENSEMBL/scripts/dump_gtf_from_gff3.pl", $log);
}

if ($xrefs) {
  $wormbase->run_script( 'generate_dbxref_file.pl', $log);
}

if ($RRID) {
  $wormbase->run_script( 'generate_RRID_data.pl', $log);
}

if ($omim) {
  $wormbase->run_script( 'generate_OMIM_xref.pl', $log);
}

&post_merge_steps                                                        if $post_merge;

&ontologies								 if $ontologies;
&make_extras                                                             if $extras;
#run some checks
$wormbase->run_script( "post_build_checks.pl -a"                 , $log) if $check;
$wormbase->run_script( "data_checks.pl -ace -gff"                , $log) if $data_check;
$wormbase->run_script( "dbcomp.pl -post_gff"                     , $log) if $data_check;
&build_release                                                           if $buildrelease;
$wormbase->run_script("finish_build.pl"                          , $log) if $finish_build;
&go_public                                                               if $public;


if ($load) {
  my $db_to_load = ($database) ? $database : $wormbase->autoace;

  $log->write_to("loading $load to $db_to_load\n");
  $log->write_to("\ttsuser = $tsuser\n\n");
  $wormbase->load_to_database( $db_to_load, $load, $tsuser ,$log ); 
  #appropriate checks are made in the Wormbase.pm
} elsif ($big_load) {
  my $db_to_load = ($database) ? $database : $wormbase->autoace;

  $log->write_to("Big-loading $big_load to $db_to_load\n");
  $log->write_to("\ttsuser = $tsuser\n\n");
  $wormbase->load_to_database( $db_to_load, $big_load, $tsuser ,$log, 1); 
}

$log->mail;

exit(0);


############################
#       SUBROUTINES        #
############################


sub make_remap_data {
  if ($wormbase->species eq 'elegans'){
    my $version = $wormbase->get_wormbase_version;
    $wormbase->run_script( "inspect-old-releases.pl -version $version -database1 ".$wormbase->database('current')." -database2 ".$wormbase->autoace, $log );
  }
}

sub do_assembly_stuff {
  if ($wormbase->species eq 'elegans'){
    $wormbase->run_script( "make_agp_file.pl",                       $log );
    $wormbase->run_script( "agp2dna.pl",                             $log ); #dependant on processGFF producing clone_acc files.
  }
}

sub map_features_to_genome {
  #
  # feature mapping, direct in perl using flanks
  #
  $wormbase->run_script( 'map_features.pl -all', $log );

  
  #
  # Oligo_sets
  #
  my $oligo_set_mapping_file = sprintf("%s/oligo_set_mappings/oligo_set_mappings.%s.ace",
                                       $wormbase->misc_dynamic,
                                       $wormbase->species);
  if (-e $oligo_set_mapping_file) {
    $log->write_to("Loading existing mappings: $oligo_set_mapping_file\n");
    $wormbase->load_to_database( $wormbase->autoace, $oligo_set_mapping_file, "OLIGO_SET_TO_GENOME", $log);
  } else{
    unless ($wormbase->species eq 'brugia' || $wormbase->species eq 'ovolvulus' || $wormbase->species eq 'sratti') {
     $log->log_and_die("Could not find Oligo_set mapping file ($oligo_set_mapping_file). Bad. Exiting\n");
    }
  }

  # check count of classes loaded
  $wormbase->run_script("check_class.pl -stage map_features -classes Sequence,Feature", $log);
 


  my $release = $wormbase->get_wormbase_version;
  my $previous_release = $release - 1;

  my $assembly_mapper = Remap_Sequence_Change->new($previous_release, $release, $wormbase->species, $wormbase->genome_diffs);

  if ($wormbase->species eq 'elegans') {    
    #
    # PCR_products - only C. elegans, but must be run before RNAiGenome
    #

    my $pcr_file = "misc_PCR_mappings_" . $wormbase->species . ".ace";
    my $pcr_mappings = $wormbase->misc_dynamic . "/" . $pcr_file;
    
    if (not $assembly_mapper->remap_test and -e $pcr_mappings) {
      # genome has not changed. Load the current mappings, and supplement with new data
      $wormbase->load_to_database( $wormbase->autoace, $pcr_mappings, "PCR_product_MISC_DYN", $log );
      $wormbase->run_script( 'PCR_product2Genome.pl -onlyunmapped', $log );
    } else {
      # genome has changed. Need to remap everything
      unlink $pcr_mappings if -e $pcr_mappings;
      $wormbase->run_script( "PCR_product2Genome.pl -acefile $pcr_mappings", $log );
    }

    #
    # RNAi
    #
    my $rnai_file = "misc_RNAi_homols_" . $wormbase->species . ".ace";
    my $rnai_mappings = $wormbase->misc_dynamic . "/" . $rnai_file;
    
    if (not $assembly_mapper->remap_test and -e $rnai_mappings) {
      # genome has not changed. Load the current mappings, and supplement with new data
      $wormbase->load_to_database( $wormbase->autoace, $rnai_mappings, "RNAi_MISC_DYN", $log );
      $wormbase->run_script( 'RNAi2Genome.pl -onlyunmapped', $log );
    } else {
      # genome has changed. Need to remap everything
      unlink $rnai_mappings if -e $rnai_mappings;
      $wormbase->run_script( "RNAi2Genome.pl -acefile $rnai_mappings", $log );
    }
    
    # check count of classes loaded
    $wormbase->run_script("check_class.pl -stage map_features_elegans -classes PCR_product,Homol_data,Sequence", $log) 
  }
 
}

sub map_features {

  # PCR products  - requires UTR GFF files
  $wormbase->run_script( 'map_PCR_products.pl', $log );
  
  #Oligo_sets
  $wormbase->run_script( 'map_Oligo_set.pl', $log );
  
  # microarray connections
  $wormbase->run_script( 'map_microarray.pl', $log );
  
  
  ## elegans only stuff
  if ($wormbase->species eq 'elegans') {
    
    # writes tables listing microarrays to genes
    $wormbase->run_script( 'make_oligo_set_mapping_table.pl -all', $log );
    
    # maps SAGE tags to the genes and to the genome
    $wormbase->run_script( 'map_tags.pl', $log );
    
    # Y2H objects
    $wormbase->run_script( 'map_Interaction.pl', $log );
    
    # RNAi experiments
    $wormbase->run_script( 'map_RNAi.pl', $log );
  }
  
  ## all species
  # TSL features
  $wormbase->run_script( 'map_feature2gene.pl', $log );
  
  # attach 'other nematode' ESTs to the genes they BLAT to best
  $wormbase->run_script( 'attach_other_nematode_ests.pl', $log );
  
}


sub map_alleles {

  ## elegans or briggsae
  if (($wormbase->species eq 'elegans') or ($wormbase->species eq 'briggsae')){
    # alleles
    $wormbase->run_script( 'split_alleles.pl', $log );
  }
  
}

#__ end map_features __#

sub remap_misc_dynamic {

  my $release = $wormbase->get_wormbase_version;
  my $previous_release = $release - 1;

  my $assembly_mapper = Remap_Sequence_Change->new($previous_release, $release, $wormbase->species, $wormbase->genome_diffs);
  if ($assembly_mapper->remap_test) {

    # remap ace files with homol_data mapped to clones
    my %clone_data = (
		      'misc_21urna_homol.ace'                 => '21_urna',
		      'misc_Expression_pattern_homol.ace'     => 'expression_pattern',
		      'misc_Tijsterman_G4.ace'                => 'Tijsterman_G4'
		      );

    foreach my $clone_data_file (keys %clone_data) {
      my $data_file = $wormbase->misc_dynamic."/$clone_data_file";
      my $backup_file = $wormbase->misc_dynamic."/BACKUP/$clone_data_file.$previous_release";
      if (-e $backup_file) {$wormbase->run_command("mv -f $backup_file $data_file", $log);}
      $wormbase->run_command("mv -f $data_file $backup_file", $log);
      $wormbase->run_script( "remap_clone_homol_data.pl -input $backup_file -out $data_file -data_type $clone_data{$clone_data_file}", $log);
    }

    # remap ace files with Feature_data mapped to clones
    my @clone_feature_data = (
		      'misc_modENCODE_Tiling_array_TARs.ace',
		      'misc_Lamm_polysomes.ace',
		      #'misc_RNASeq_hits_elegans.ace',
		      );

    foreach my $clone_feature_data_file (@clone_feature_data) {
      my $data_file = $wormbase->misc_dynamic."/$clone_feature_data_file";
      my $backup_file = $wormbase->misc_dynamic."/BACKUP/$clone_feature_data_file.$previous_release";
      if (-e $backup_file) {$wormbase->run_command("mv -f $backup_file $data_file", $log);}
      $wormbase->run_command("mv -f $data_file $backup_file", $log);
      $wormbase->run_script( "remap_clone_feature_data.pl -input $backup_file -out $data_file", $log);
    }


    # remap twinscan
    my $twinscan = $wormbase->misc_dynamic."/misc_twinscan.ace";
    my $backup_twinscan = $wormbase->misc_dynamic."/BACKUP/misc_twinscan.ace.$previous_release";
    if (-e $backup_twinscan) {$wormbase->run_command("mv -f $backup_twinscan $twinscan", $log);}
    $wormbase->run_command("mv -f $twinscan $backup_twinscan", $log);
    $wormbase->run_script( "remap_twinscan_between_releases.pl -release1 $previous_release -release2 $release -twinscan $backup_twinscan -out $twinscan", $log);

    # remap genefinder
    my $genefinder = $wormbase->misc_dynamic."/misc_genefinder.ace";
    my $backup_genefinder = $wormbase->misc_dynamic."/BACKUP/misc_genefinder.ace.$previous_release";
    if (-e $backup_genefinder) {$wormbase->run_command("mv -f $backup_genefinder $genefinder", $log);}
    $wormbase->run_command("mv -f $genefinder $backup_genefinder", $log);
    $wormbase->run_script( "remap_genefinder_between_releases.pl -input $backup_genefinder -out $genefinder", $log);

    # remap jigsaw - uses the remap genefinder script
    my $jigsaw = $wormbase->misc_dynamic."/misc_jigsaw.ace";
    my $backup_jigsaw = $wormbase->misc_dynamic."/BACKUP/misc_jigsaw.ace.$previous_release";
    if (-e $backup_jigsaw) {$wormbase->run_command("mv -f $backup_jigsaw $jigsaw", $log);}
    $wormbase->run_command("mv -f $jigsaw $backup_jigsaw", $log);
    $wormbase->run_script( "remap_genefinder_between_releases.pl -input $backup_jigsaw -out $jigsaw", $log);

    # remap mGene - uses the remap genefinder script
    my $mgene = $wormbase->misc_dynamic."/misc_mgene.ace";
    my $backup_mgene = $wormbase->misc_dynamic."/BACKUP/misc_mgene.ace.$previous_release";
    if (-e $backup_mgene) {$wormbase->run_command("mv -f $backup_mgene $mgene", $log);}
    $wormbase->run_command("mv -f $mgene $backup_mgene", $log);
    $wormbase->run_script( "remap_genefinder_between_releases.pl -input $backup_mgene -out $mgene", $log);

    # remap Hillier RNASEQ_CDS - uses the remap genefinder script
    my $rnaseq_cds = $wormbase->misc_dynamic."/misc_RNASEQ_CDS.ace";
    my $backup_rnaseq_cds = $wormbase->misc_dynamic."/BACKUP/misc_RNASEQ_CDS.ace.$previous_release";
    if (-e $backup_rnaseq_cds) {$wormbase->run_command("mv -f $backup_rnaseq_cds $rnaseq_cds", $log);}
    $wormbase->run_command("mv -f $rnaseq_cds $backup_rnaseq_cds", $log);
    $wormbase->run_script( "remap_genefinder_between_releases.pl -input $backup_rnaseq_cds -out $rnaseq_cds", $log);

    # remap modENCODE misc_modENCODE_aggregate_transcripts.ace - uses the remap genefinder script
    my $aggregate_cds = $wormbase->misc_dynamic."/misc_modENCODE_aggregate_transcripts.ace";
    my $backup_aggregate_cds = $wormbase->misc_dynamic."/BACKUP/misc_modENCODE_aggregate_transcripts.ace.$previous_release";
    if (-e $backup_aggregate_cds) {$wormbase->run_command("mv -f $backup_aggregate_cds $aggregate_cds", $log);}
    $wormbase->run_command("mv -f $aggregate_cds $backup_aggregate_cds", $log);
    $wormbase->run_script( "remap_genefinder_between_releases.pl -input $backup_aggregate_cds -out $aggregate_cds", $log);

    # remap fosmids
    my $fosmids = $wormbase->misc_dynamic."/fosmids.ace";
    my $backup_fosmids = $wormbase->misc_dynamic."/BACKUP/fosmids.ace.$previous_release";
    if (-e $backup_fosmids) {$wormbase->run_command("mv -f $backup_fosmids $fosmids", $log);}
    $wormbase->run_command("mv -f $fosmids $backup_fosmids", $log);
    $wormbase->run_script( "remap_fosmids_between_releases.pl -input $backup_fosmids -out $fosmids", $log);
   
    # the TEC-REDs are placed back on the genome by using the location of the Features they defined
    $wormbase->run_script( "map_tec-reds.pl", $log);

    # remap and copy over the SUPPLEMENTARY_GFF dir from BUILD_DATA
    my $sup_dir = $wormbase->misc_dynamic."/SUPPLEMENTARY_GFF";
    my $release = $wormbase->version;
    my $old_release = $release - 1;
    my $backup_dir = "$sup_dir/BACKUP_${old_release}";

    foreach my $file (glob("$sup_dir/elegans.*.gff2"), glob("$sup_dir/elegans.*.gff3")) {
      my ($fname) = $file =~ /$sup_dir\/(\S+)$/;
      my $backup_file = "$sup_dir/$fname";
      $wormbase->run_command("mv -f $file $backup_file", $log);
      $wormbase->run_script("remap_gff_between_releases.pl -gff $backup_file -output $file -release1 $old_release -release2 $release", $log);
    }
  } else {
    $log->write_to("Assembly has apparently not changed, so did not remap the MISC_DYNAMIC data\n");
  }   
}

#__ end remap_misc_dynamic __#

sub make_UTR {
  my ($log)=@_;
  
  $wormbase->run_script("make_UTR_GFF.pl", $log);
  $wormbase->run_script("make_UTR_GFF.pl -gff3", $log);
}

sub post_merge_steps {
  $wormbase->run_script("molecular_names_for_genes.pl", $log);
  $wormbase->run_script("ONTOLOGY/update_GO_terms.pl", $log); 
  $wormbase->run_script('transfer_interpro_GO_terms.pl', $log );
  $wormbase->run_script("cluster_gene_connection.pl", $log);
  $wormbase->run_script("write_DB_remark.pl", $log);
  $wormbase->run_script("tier3_stubs.pl", $log);
}

sub ontologies {
  $wormbase->run_script( "ONTOLOGY/make_anatomy_GAF.pl", $log);
  $wormbase->run_script( "ONTOLOGY/make_phenotype_GAF.pl", $log);
  $wormbase->run_script( "ONTOLOGY/make_disease_GAF.pl", $log);
  $wormbase->run_script( "ONTOLOGY/make_lifestage_GAF.pl", $log);
  $wormbase->run_script( "ONTOLOGY/make_GO_GAF.pl", $log);
  $wormbase->run_script( "ONTOLOGY/get_easy_phenotypes.pl", $log);
  $wormbase->run_script( "AGR/make_agr_disease_json.pl -writedaf", $log);
}

sub make_extras {
  my $version = $wormbase->get_wormbase_version;
  $wormbase->run_script( "make_keysets.pl -all", $log);
  $wormbase->run_script( "genestats.pl" , $log);
}

sub build_release {
  $wormbase->run_script( "build_release_files.pl", $log);
  $wormbase->run_script( "make_assembly_manifest.pl", $log);
  $wormbase->run_script( "release_letter.pl"  , $log);
}


sub go_public {

  my $ftp_release_dir = $wormbase->ftp_site . "/releases";
  my $ftp_staging_dir = $wormbase->ftp_site . "/staging/releases";
  my $db_dir = $wormbase->wormpub . "/DATABASES";
  my $rel   = $wormbase->get_wormbase_version_name;
  
  if (not -d "$ftp_staging_dir/$rel") {
    $log->log_and_die("Did not find $ftp_staging_dir/$rel. Something wrong. Not going public\n");
  }
  if (not -d "$db_dir/$rel") {
    $log->log_and_die("Did not find $db_dir/$rel. Something wrong. Not going public\n");
  }
  if (not -e "$ftp_staging_dir/$rel/letter.$rel") {
    $log->log_and_die("Did not find $ftp_staging_dir/$rel/letter.$rel. Something wrong. Not going public\n");
  }
  
  $log->write_to("Moving the release folder from staging to live\n");
  $wormbase->run_command("mv $ftp_staging_dir/$rel $ftp_release_dir/", $log) 
      and $log->log_and_die("Failed to mv release folder into place - aborting\n");
  
  $log->write_to("Updating the current_development symlink\n");
  eval {
    $wormbase->run_command("rm -f $ftp_release_dir/current-development-release", $log) and die;
    $wormbase->run_command("cd $ftp_release_dir && ln -s $rel current-development-release", $log) and die;
  };
  $@ and $log->write_to("WARNING: Failed to update FTP development symlink - need to fix manually\n");
  
  $log->write_to("Updating the current_DB symlink\n");
  eval {
    $wormbase->run_command("rm -f $db_dir/current_DB", $log) and die;
    $wormbase->run_command("cd $db_dir && ln -s $rel current_DB", $log) and die;
  };
  $@ and $log->write_to("WARNING: Failed to update current_DB symlink - need to fix manually\n");
  
  $log->write_to("Sending release letter to staff\n");
  my $letter = "$ftp_release_dir/$rel/letter.$rel";
  $wormbase->mail_maintainer( "WormBase $rel release",
                              'staff@wormbase.org',
                              $letter);
  $log->write_to("\n\nWormBase $rel has been (data) has been released!\n");
}
