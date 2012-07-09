#!/software/bin/perl -w
#
# autoace_builder.pl
#
# based on the original autoace_minder.pl
#
# Usage : autoace_builder.pl [-options]
#
# Last edited by: $Author: klh $
# Last edited on: $Date: 2012-07-09 08:19:44 $

my $script_dir = $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'};
use lib '/software/worm/lib/site_perl';

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
my ( $load, $tsuser, $map_features, $remap_misc_dynamic, $map, $transcripts, $intergenic, $misc_data_sets, $homol_data_sets, $nem_contigs);
my ( $GO_term, $rna , $dbcomp, $confirm, $operon ,$repeats, $remarks, $names, $treefam, $cluster);
my ( $utr, $agp, $gff_munge, $extras , $ontologies, $interpolate, $check, $enaseqxrefs, $enaprotxrefs, $xrefs);
my ( $data_check, $buildrelease, $public,$finish_build, $gffdb, $autoace, $release, $user, $kegg);


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
	   'finish_wormpep' => \$finish_wormpep,
	   'gff_dump:s'     => \$gff_dump,
	   'processGFF:s'   => \$processGFF,
	   'gff_split'      => \$gff_split,
	   'gene_span'      => \$gene_span,
	   'load=s'         => \$load,
	   'prep_blat'      => \$prep_blat,
	   'run_blat'       => \$run_blat,
	   'finish_blat'    => \$finish_blat,
	   'tsuser=s'       => \$tsuser,
	   'map'            => \$map,
	   'remap_misc_dynamic' => \$remap_misc_dynamic,
	   'map_features'   => \$map_features,
	   'transcripts'    => \$transcripts,
	   'intergenic'     => \$intergenic,
	   'nem_contig'     => \$nem_contigs,
	   'misc_data_sets' => \$misc_data_sets,
	   'homol_data_sets'=> \$homol_data_sets,
           'sequencexrefs'  => \$enaseqxrefs,
           'proteinxrefs'   => \$enaprotxrefs,
           'xrefs'          => \$xrefs,
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
	   'ontologies'     => \$ontologies,
	   'buildrelease'   => \$buildrelease,
	   'public'         => \$public,
	   'finish_build'   => \$finish_build,
	   'gffdb'          => \$gffdb,
	   'autoace'        => \$autoace,
	   'release'        => \$release,
	   'check'    	    => \$check,
	   'data_check'     => \$data_check,
	   'species:s'      => \$species,
	   'user:s'         => \$user,
	   'kegg'           => \$kegg,
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
    $wormbase->run_script("check_class.pl -stlace -camace -genace -csh -caltech -misc_static -brigace -stage init", $log);
    $wormbase->run_script("check_class.pl -incomplete -stage incomplete", $log);
  }
}


$wormbase->run_script( "chromosome_dump.pl --dna --composition", $log ) if $first_dumps;
$wormbase->run_script("update_Common_data.pl -clone2centre -clone2acc -clone2size -clone2dbid -clone2seq $species -genes2lab -worm_gene2cgc -worm_gene2geneID -worm_gene2class -est -est2feature -gene_id -clone2type -cds2cgc -rna2cgc -pseudo2cgc ", $log ) if $first_dumps;

$wormbase->run_script( "build_dumpGFF.pl -stage $gff_dump",                $log ) if $gff_dump;

$wormbase->run_script( "get_ena_submission_xrefs.pl -sequencexrefs -load", $log)  if $enaseqxrefs;
$wormbase->run_script( "processGFF.pl -$processGFF",                       $log ) if $processGFF;

&do_assembly_stuff() if $assembly;   # dependant on clone_acc for agp

$wormbase->run_script( 'make_wormpep.pl -initial -all',                    $log ) if $make_wormpep;

&map_features_to_genome() if $map_features;

#########   BLAT  ############
$wormbase->run_script( 'BLAT_controller.pl -dump', $log ) if $prep_blat;
#//--------------------------- batch job submission -------------------------//
$wormbase->run_script( 'BLAT_controller.pl -run', $log )        if $run_blat;
#//--------------------------- batch job submission -------------------------//
$wormbase->run_script( 'BLAT_controller.pl -postprocess -process -intron -load', $log ) if $finish_blat;
#//--------------------------- batch job submission -------------------------//
# $build_dumpGFF.pl; (blat) is run chronologically here but previous call will operate

$wormbase->run_script( 'batch_transcript_build.pl', $log) if $transcripts;
#requires GFF dump of transcripts (done within script if all goes well)

$wormbase->run_script( 'WBGene_span.pl'                   , $log ) if $gene_span;
&make_UTR($log)                                                    if $utr;

$wormbase->run_script( 'find_intergenic.pl'               , $log ) if $intergenic;

##  Horrid Geneace related stuff  ##########
#make_pseudo_map_positions.pl -load
#get_interpolated_gmap.pl
#update_inferred_multi_pt.pl -load

####### mapping part ##########
&map_features                                                            if $map;

&remap_misc_dynamic                                                      if $remap_misc_dynamic;

&get_repeats                                                             if $repeats; # loaded with homols
#must have farm complete by this point.
$wormbase->run_script( 'load_data_sets.pl -misc', $log) if $misc_data_sets;
$wormbase->run_script( 'load_data_sets.pl -homol', $log) if $homol_data_sets;
# $build_dumpGFF.pl; (homol) is run chronologically here but previous call will operate
$wormbase->run_script( 'make_wormrna.pl'                         , $log) if $rna;
$wormbase->run_script( 'confirm_genes.pl -load'                  , $log) if $confirm;
$wormbase->run_script( 'map_operons.pl'                          , $log) if $operon;
$wormbase->run_script( "get_ena_submission_xrefs.pl -proteinxrefs -load", $log) if $enaprotxrefs;
$wormbase->run_script( 'make_wormpep.pl -all -final'                  , $log) if $finish_wormpep;
$wormbase->run_script( 'write_DB_remark.pl'                      , $log) if $remarks;
$wormbase->run_script( 'molecular_names_for_genes.pl'            , $log) if $names;
$wormbase->run_script( 'get_treefam.pl'                          , $log) if $treefam;
$wormbase->run_script( 'cluster_gene_connection.pl'              , $log) if $cluster;
$wormbase->run_script( 'inherit_GO_terms.pl -phenotype -motif -tmhmm', $log ) if $GO_term;
$wormbase->run_script( 'KEGG.pl', $log )                                 if $kegg;

# $build_dumpGFF.pl; (final) is run chronologically here but previous call will operate
# $wormbase->run_script( "processGFF.pl -$processGFF",        $log ) if $processGFF;    #nematode - to add species to nematode BLATs

$wormbase->run_script( "interpolation_manager.pl"                , $log) if $interpolate;
$wormbase->run_script( "make_agp_file.pl"                        , $log) if $agp;

#several GFF manipulation steps
if ($gff_munge) {

  if ($wormbase->species eq 'elegans') {
    $wormbase->run_script( 'landmark_genes2gff.pl', $log); # this has to be run before GFFmunger.pl
    $wormbase->run_script( 'web_data/interpolate_gmap2pmap.pl', $log);
  }
  $wormbase->run_script( 'web_data/map_translated_features_to_genome.pl', $log);

  ###
  $wormbase->run_script( 'GFFmunger.pl -all', $log); # GFFmunger uses the files created by the previous scripts
  ###

  $wormbase->run_script( 'over_load_SNP_gff.pl' , $log);
  $wormbase->run_script( 'overload_rnai.pl'     , $log);
  $wormbase->run_script( 'overload_operon.pl' , $log);
  
  if ($wormbase->assembly_type eq 'chromosome') {
    $wormbase->run_script( "Map_pos_GFFprocess.pl", $log);
  }
  if ($wormbase->species eq 'elegans') {
    $wormbase->run_script( "chromosome_script_lsf_manager.pl -command '/software/bin/perl $ENV{'CVS_DIR'}/process_sage_gff.pl' -mito -prefix", $log);
  }
}
if ($xrefs) {
  if ($wormbase->species eq 'elegans' or $wormbase->species eq 'briggsae') {
    $wormbase->run_script( 'generate_dbxref_file.pl', $log);
  } else {
    $log->write_to("This should only be run for elegans and briggsae, so not doing anything\n");
  }
}

&ontologies								 if $ontologies;
&make_extras                                                             if $extras;
#run some checks
$wormbase->run_script( "post_build_checks.pl -a"                 , $log) if $check;
$wormbase->run_script( "data_checks.pl -ace -gff"                , $log) if $data_check;
$wormbase->run_script( "dbcomp.pl"                               , $log) if $data_check;
$wormbase->run_script( "build_release_files.pl"                  , $log) if $buildrelease;
&public_sites                                                            if $public;
&release                                                                 if $release;

$wormbase->run_script("finish_build.pl"                          , $log) if $finish_build;
# Update the gffdb to reflect the database being built.
if  ($gffdb && $autoace) {
  $wormbase->run_command("update_gffdb.csh -autoace"               , $log);
}
else {$wormbase->run_command("update_gffdb.csh"                  , $log) if $gffdb;
}

if ($load) {
    $log->write_to("loading $load to ".$wormbase->autoace."\n");
    $log->write_to("\ttsuser = $tsuser\n\n");
    $wormbase->load_to_database( $wormbase->autoace, $load, $tsuser ,$log); #appropriate checks are made in the Wormbase.pm
}

$log->mail;

exit(0);


############################
#       SUBROUTINES        #
############################

sub do_assembly_stuff {

  if ($wormbase->species eq 'elegans'){
    my $version = $wormbase->get_wormbase_version;
    $wormbase->run_script( "inspect-old-releases.pl -version $version -database1 ".$wormbase->database('current')." -database2 ".$wormbase->autoace, $log );
    
    $wormbase->run_script( "make_agp_file.pl",                       $log );
    $wormbase->run_script( "agp2dna.pl",                             $log ); #dependant on processGFF producing clone_acc files.
    
    my $agp_errors = 0;
    
    foreach my $chrom ($wormbase->get_chromosome_names) {
      open( AGP, "<" . $wormbase->autoace . "/yellow_brick_road/CHROMOSOME_${chrom}.agp_seq.log" )
          or die "Couldn't open agp file : $!";
      while (<AGP>) {
        $agp_errors++ if (/ERROR/);
      }
      close(AGP);
    }
    
    $log->write_to("ERRORS ( $agp_errors ) in agp file\n");
  }
}

sub map_features_to_genome {
  #
  # feature mapping, direct in perl using flanks
  #
  $wormbase->run_script( 'map_features.pl -all', $log );

  my $release = $wormbase->get_wormbase_version;
  my $previous_release = $release - 1;

  my $assembly_mapper = Remap_Sequence_Change->new($previous_release, $release, $wormbase->species, $wormbase->genome_diffs);

  #
  # PCR_products
  #
  my $pcr_file = "misc_PCR_mappings_" . $wormbase->species . ".ace";
  my $pcr_mappings = $wormbase->misc_dynamic . "/" . $pcr_file;

  if (not $assembly_mapper->remap_test and -e $pcr_mappings) {
    # genome has not changed. Load the current mappings, and supplement with new data
    $wormbase->load_to_database( $wormbase->autoace, $pcr_mappings, "PCR_product_MISC_DYN", $log );
    $wormbase->run_script( 'PCR_product2Genome.pl', $log );
  } else {
    # genome has changed. Need to remap everything
    unlink $pcr_mappings if -e $pcr_mappings;
    $wormbase->run_script( 'PCR_product2Genome.pl -acefile $pcr_mappings', $log );
  }
  
  #
  # RNAi
  #
  my $rnai_file = "misc_RNAi_homols_" . $wormbase->species . "ace";
  my $rnai_mappings = $wormbase->misc_dynamic . "/" . $rnai_file;

  if (not $assembly_mapper->remap_test and -e $rnai_mappings) {
    # genome has not changed. Load the current mappings, and supplement with new data
    $wormbase->load_to_database( $wormbase->autoace, $rnai_mappings, "RNAi_MISC_DYN", $log );
    $wormbase->run_script( 'RNAi2Genome.pl', $log );
  } else {
    # genome has changed. Need to remap everything
    unlink $rnai_mappings if -e $rnai_mappings;
    $wormbase->run_script( 'RNAi2Genome.pl -acefile $rnai_mappings', $log );
  }
  
}

sub map_features {

    ## elegans only stuff
    if ($wormbase->species eq 'elegans') {
	# PCR products  - requires UTR GFF files
	$wormbase->run_script( 'map_PCR_products.pl', $log );

	#Oligo_sets
	$wormbase->run_script( 'map_Oligo_set.pl', $log );


	# writes tables listing microarrays to genes
	$wormbase->run_script( 'make_oligo_set_mapping_table.pl -all', $log );

	# maps SAGE tags to the genes and to the genome
	$wormbase->run_script( 'map_tags.pl -load', $log );

	# microarray connections
	$wormbase->run_script( 'map_microarray.pl -load', $log );
	
	# Y2H objects
	$wormbase->run_script( 'map_Y2H.pl -load', $log );


    	# RNAi experiments
	$wormbase->run_script( 'map_RNAi.pl -load', $log );
    }

    ## elegans or briggsae
    if (($wormbase->species eq 'elegans') or ($wormbase->species eq 'briggsae')){
	# alleles
	$wormbase->run_script( 'split_alleles.pl', $log );
    }

    ## all species
    # TSL features
    $wormbase->run_script( 'map_feature2gene.pl -load', $log );
    
    # attach 'other nematode' ESTs to the genes they BLAT to best
    $wormbase->run_script( 'attach_other_nematode_ests.pl -load', $log );
    
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
		      'misc_mass_spec_GenniferMerrihew.ace'   => 'mass_spec',
		      );

    foreach my $clone_data_file (keys %clone_data) {
      my $data_file = $wormbase->misc_dynamic."/$clone_data_file";
      my $backup_file = $wormbase->misc_dynamic."/BACKUP/$clone_data_file.$previous_release";
      if (-e $backup_file) {$wormbase->run_command("mv -f $backup_file $data_file", $log);}
      $wormbase->run_command("mv -f $data_file $backup_file", $log);
      $wormbase->run_script( "remap_clone_homol_data.pl -input $backup_file -out $data_file -data_type $clone_data{$clone_data_file}", $log);
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

    # remap fosmids
    my $fosmids = $wormbase->misc_dynamic."/fosmids.ace";
    my $backup_fosmids = $wormbase->misc_dynamic."/BACKUP/fosmids.ace.$previous_release";
    if (-e $backup_fosmids) {$wormbase->run_command("mv -f $backup_fosmids $fosmids", $log);}
    $wormbase->run_command("mv -f $fosmids $backup_fosmids", $log);
    $wormbase->run_script( "remap_fosmids_between_releases.pl -input $backup_fosmids -out $fosmids", $log);
   
    # the TEC-REDs are placed back on the genome by using the location of the Features they defined
    $wormbase->run_script( "map_tec-reds.pl", $log);

    # remap and copy over the SUPPLEMENTARY_GFF dir from BUILD_DATA
    my $sup_dir = $wormbase->build_data."/SUPPLEMENTARY_GFF";
    my $backup_dir = "$sup_dir/BACKUP";
    my $release = $wormbase->version;
    my $old_release = $release - 1;
    opendir(DIR,$sup_dir) or $log->log_and_die("cant open $sup_dir: $!\n");
    while ( my $file = readdir( DIR ) ) {
      next unless( $file =~ /gff$/ );
      my $gff = "$sup_dir/$file";
      my $backup_gff = "$backup_dir/$file.$old_release";
      if (-e $backup_gff) {$wormbase->run_command("mv -f $backup_gff $gff", $log);}
      $wormbase->run_command("mv -f $gff $backup_gff", $log);
      $wormbase->run_script("remap_gff_between_releases.pl -gff $backup_gff -output $gff -release1 $old_release -release2 $release", $log);
    }
    closedir DIR;

  } else {
    $log->write_to("Assembly has apparently not changed, so did not remap the MISC_DYNAMIC data\n");
  }

  # the SUPPLEMENTARY_GFF directory is copied over whether or not it has been remapped
  $wormbase->run_command("cp -R ".$wormbase->build_data."/SUPPLEMENTARY_GFF ".$wormbase->sequences."/", $log);
   
}

#__ end remap_misc_dynamic __#

sub make_UTR {
  my ($log)=@_;
  
  $log->write_to("bsub commands . . . . \n\n");
  my $lsf = LSF::JobManager->new(-J => "make_UTRs", -o => "/dev/null");

  my $out_dir = $wormbase->gff_splits;
  my (@commands, @out_files);

  if ($wormbase->assembly_type eq 'contig') {
    my @chrs = $wormbase->get_chromosome_names;
    my $chunk_total = 24;
    $chunk_total = scalar(@chrs) if $chunk_total > scalar(@chrs);

    foreach my $chunk_id (1..$chunk_total) {
      my $outfile = "$out_dir/UTR.chunk_${chunk_id}.gff";
      push @commands, "make_UTR_GFF.pl -chunkid $chunk_id -chunktotal $chunk_total -outfile $outfile";
      push @out_files, $outfile;
    }
  } else {
    foreach my $chrom ($wormbase->get_chromosome_names(-mito => 1, -prefix => 1)) {
      my $outfile = "$out_dir/${chrom}_UTR.gff";
      push @commands, "make_UTR_GFF.pl -chrom $chrom -outfile $outfile";
    }
  }

  foreach my $cmd (@commands) {
    $log->write_to("$cmd\n");
    $cmd = $wormbase->build_cmd($cmd);
    $lsf->submit($cmd);
  }
  $lsf->wait_all_children( history => 1 );
  $log->write_to("All UTR jobs have completed!\n");
  for my $job ( $lsf->jobs ) {
    $log->error("Job $job (" . $job->history->command . ") exited with LSF error code: ". $job->history->exit_status ."\n") if $job->history->exit_status != 0;
  }
  $lsf->clear;   
  
  foreach my $outfile (@out_files) {
    $log->error("$outfile is missing or empty") if not -e $outfile or not -s $outfile;
  }

  #merge into single file.
  if($wormbase->assembly_type eq 'contig') {
    $wormbase->run_command("cat @out_files  > $out_dir/UTR.gff",$log);
    $wormbase->run_command("rm -f @out_files",$log);
	
# check the files      
#Crem_Contig35   Coding_transcript       three_prime_UTR 74394   74463   .       -       .       Transcript "CRE24456"
#Crem_Contig35   Coding_transcript       coding_exon     74464   75307   .       -       1       Transcript "CRE24456" ; CDS "CRE24456"
#Crem_Contig35   Coding_transcript       five_prime_UTR  75458   75463   .       -       .       Transcript "CRE24456"

    $wormbase->check_file($wormbase->gff_splits."/UTR.gff", $log,
			  lines => ['^##', 
				    "^\\S+\\s+Coding_transcript\\s+(three_prime_UTR|coding_exon|five_prime_UTR)\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+\\s+Transcript\\s+\\S+"],
			  gff => 1,
			 );   
  } else {
  
    # check the files      
    foreach my $sequence ( $wormbase->get_chromosome_names(-prefix => 1, -mito => 1) ) {
      if($wormbase->species eq 'elegans') {

	my %sizes = (
		     'CHROMOSOME_I'       => 3400000,
		     'CHROMOSOME_II'      => 3500000,
		     'CHROMOSOME_III'     => 3000000,
		     'CHROMOSOME_IV'      => 3500000,
		     'CHROMOSOME_V'       => 4200000,
		     'CHROMOSOME_X'       => 3400000,
		     'CHROMOSOME_MtDNA'   =>    2000,
		    );
	$wormbase->check_file($wormbase->gff_splits."/${sequence}_UTR.gff", $log,
			      minsize => $sizes{$sequence},
			      lines => ['^##', 
					"^\\S+\\s+Coding_transcript\\s+(three_prime_UTR|coding_exon|five_prime_UTR)\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+\\s+Transcript\\s+\\S+"],
			      gff => 1,
			     );   
      } elsif ($wormbase->species eq 'briggsae') {
       
	my %sizes = (
		     'chr_I'          => 1300000,
		     'chr_I_random'   =>  300000,
		     'chr_II'         => 1500000,
		     'chr_II_random'  =>  200000,
		     'chr_III'        => 1500000,
		     'chr_III_random' =>   70000,
		     'chr_IV'         => 1600000,
		     'chr_IV_random'  =>   50000,
		     'chr_V'          => 1900000,
		     'chr_V_random'   =>  300000,
		     'chr_X'          => 2100000,
		     'chr_Un'         =>  600000,
		    );
	$wormbase->check_file($wormbase->gff_splits."/${sequence}_UTR.gff", $log,
			      minsize => $sizes{$sequence},
			      lines => ['^##', 
					"^\\S+\\s+Coding_transcript\\s+(three_prime_UTR|coding_exon|five_prime_UTR)\\s+\\d+\\s+\\d+\\s+\\S+\\s+[-+\\.]\\s+\\S+\\s+Transcript\\s+\\S+"],
			      gff => 1,
			     );   
      }
    }

  }
}


sub get_repeats {
  #repeatmasked chromosomes
  my $species= lc (ref $wormbase);
  $wormbase->run_script("get_repeatmasked_chroms.pl -database worm_ensembl_$species", $log);
  $wormbase->run_script("get_repeatmasked_chroms.pl -database worm_ensembl_$species -softmask", $log);

  #inverted
  $wormbase->run_script("run_inverted.pl -all" , $log);
}


sub ontologies {
	$wormbase->run_script( "ONTOLOGY/parse_expr_pattern_new.pl", $log);
	$wormbase->run_script( "ONTOLOGY/parse_go_terms_new.pl -rnai -gene", $log);
	$wormbase->run_script( "ONTOLOGY/parse_phenotype_new.pl", $log);
	$wormbase->run_script( "ONTOLOGY/get_easy_phenotypes.pl", $log);
}

sub make_extras {
  my $version = $wormbase->get_wormbase_version;
  $wormbase->run_script( "make_keysets.pl -all -history $version", $log);
  $wormbase->run_script( "genestats.pl" , $log);
}


sub public_sites {
  # gets everything on the to FTP and websites and prepares release letter ready for final edit and sending.
  $wormbase->run_script( "make_FTP_sites.pl -all", $log);
  $wormbase->run_script( "release_letter.pl -l"  , $log);
}


sub release {
  # copy and send the release letter around
  $wormbase->run_script( "distribute_letter.pl", $log);
  
  # Make data on FTP site available
  $log->write_to("Updating symlink on FTP site\n");

  my $ftpdir = $wormbase->ftp_site;
  my $rel   = $wormbase->get_wormbase_version_name;
  $wormbase->run_command("rm -f $ftpdir/development_release", $log);
  $wormbase->run_command("cd $ftpdir; ln -s releases/$rel development_release", $log);  
}
