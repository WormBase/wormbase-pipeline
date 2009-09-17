#!/nfs/disk100/wormpub/bin/perl
#
# autoace_builder.pl
#
# based on the original autoace_minder.pl
#
# Usage : autoace_builder.pl [-options]
#
# Last edited by: $Author: ar2 $
# Last edited on: $Date: 2009-09-17 10:52:08 $

my $script_dir = $ENV{'CVS_DIR'};
use lib $ENV{'CVS_DIR'};

use strict;
use Wormbase;
use Getopt::Long;
use File::Copy;
use Coords_converter;
use Log_files;
use Storable;
use LSF RaiseError => 0, PrintError => 1, PrintOutput => 0;
use LSF::JobManager;

my ( $debug, $test, $database, $species);
my ( $initiate, $prepare_databases, $acefile, $build, $first_dumps );
my ( $make_wormpep, $finish_wormpep );
my ( $prep_blat, $run_blat,     $finish_blat );
my ( $gff_dump,     $processGFF, $gff_split );
my $gene_span;
my ( $load, $tsuser, $map_features, $remap_misc_dynamic, $map, $transcripts, $intergenic, $data_sets, $nem_contigs);
my ( $GO_term, $rna , $dbcomp, $confirm, $operon ,$repeats, $remarks, $names, $treefam, $cluster);
my ( $utr, $agp, $gff_munge, $extras , $ontologies, $interpolate, $check);
my ( $data_check, $buildrelease, $public,$finish_build, $release, $user, $kegg);


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
	   'ontologies'     => \$ontologies,
	   'buildrelease'   => \$buildrelease,
	   'public'         => \$public,
	   'finish_build'   => \$finish_build,
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
$wormbase->run_script( 'make_acefiles.pl',                  $log ) if $acefile;
$wormbase->run_script( 'make_autoace.pl',                   $log ) if $build;

# check the primary databases once all files have been loaded and you are ready for the 1st gff dumps.
$wormbase->run_script( "check_primary_database.pl -organism ${\$wormbase->species}", $log ) if ($gff_dump eq 'init');

#//--------------------------- batch job submission -------------------------//
$wormbase->run_script( "build_dumpGFF.pl -stage $gff_dump", $log ) if $gff_dump;      #init
$wormbase->run_script( "processGFF.pl -$processGFF",        $log ) if $processGFF;    #clone_acc
&first_dumps                                                       if $first_dumps;   # dependant on clone_acc for agp
$wormbase->run_script( 'make_wormpep.pl -initial -all',          $log ) if $make_wormpep;
$wormbase->run_script( 'map_features.pl -all',              $log ) if $map_features;


#########   BLAT  ############
$wormbase->run_script( 'BLAT_controller.pl -dump', $log ) if $prep_blat;
#//--------------------------- batch job submission -------------------------//
$wormbase->run_script( 'BLAT_controller.pl -run', $log )        if $run_blat;
#//--------------------------- batch job submission -------------------------//
$wormbase->run_script( 'BLAT_controller.pl -virtual -process -postprocess -intron -load', $log ) if $finish_blat;
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
$wormbase->run_script( 'load_data_sets.pl -homol -misc', $log) if $data_sets;
# $build_dumpGFF.pl; (homol) is run chronologically here but previous call will operate
$wormbase->run_script( 'make_wormrna.pl'                         , $log) if $rna;
$wormbase->run_script( 'confirm_genes.pl -load'                  , $log) if $confirm;
$wormbase->run_script( 'map_operons.pl'                          , $log) if $operon;
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
 $wormbase->run_script( "landmark_genes2gff.pl"                   , $log) if ($wormbase->species eq 'elegans');
 $wormbase->run_script( "GFFmunger.pl -all"                       , $log);
 $wormbase->run_script( "over_load_SNP_gff.pl"                    , $log);
 # these are currently elegans specific.
 if ($wormbase->species eq 'elegans') {
   $wormbase->run_script( "Map_pos_GFFprocess.pl"                  , $log);	# must be run before CGH_GFFprocess.pl as it looks for 'Allele' lines
   $wormbase->run_script( "CGH_GFFprocess.pl"                      , $log);
   $wormbase->run_script( "chromosome_script_lsf_manager.pl -command '/software/bin/perl $ENV{'CVS_DIR'}/process_sage_gff.pl' -mito -prefix", $log);
 }
}

&ontologies								if $ontologies;
&make_extras                                                             if $extras;
#run some checks
$wormbase->run_script( "post_build_checks.pl -a"                 , $log) if $check;
$wormbase->run_script( "data_checks.pl -ace -gff"                , $log) if $data_check;
$wormbase->run_script( "dbcomp.pl"                               , $log) if $data_check;
$wormbase->run_script( "build_release_files.pl"                  , $log) if $buildrelease;
&public_sites                                                            if $public;
$wormbase->run_script( "distribute_letter.pl"                    , $log) if $release;

$wormbase->run_script("finish_build.pl"                          , $log) if $finish_build;
$wormbase->run_command("update_gffdb.csh"                         , $log) if $finish_build;

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

sub first_dumps {

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

sub map_features {

    ## elegans only stuff
    if ($wormbase->species eq 'elegans') {
	# PCR products  - requires UTR GFF files
	$wormbase->run_script( 'map_PCR_products.pl', $log );

	#Oligo_sets
	$wormbase->run_script( 'map_Oligo_set.pl', $log );

	# RNAi experiments
	$wormbase->run_script( 'map_RNAi.pl -load', $log );

	# writes tables listing microarrays to genes
	$wormbase->run_script( 'make_oligo_set_mapping_table.pl -all', $log );

	# maps SAGE tags to the genes and to the genome
	$wormbase->run_script( 'map_tags.pl -load', $log );

	# microarray connections
	$wormbase->run_script( 'map_microarray.pl -load', $log );
	
	# Y2H objects
	$wormbase->run_script( 'map_Y2H.pl -load', $log );
    }

    ## elegans or briggsae
    if (($wormbase->species eq 'elegans') or ($wormbase->species eq 'briggsae')){
	# alleles
	$wormbase->run_script( 'map_Alleles.pl', $log );
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

  # test to see if we need to run the remapping programs
  $wormbase->run_script( "test_remap_between_releases.pl -release1 $previous_release -release2 $release", $log );
  my $flag = "/tmp/remap_elegans_data";
  open(FLAG, "< $flag") || die "Could not open the file $flag\n";
  my $answer = <FLAG>;
  close(FLAG);
  chomp $answer;

  if ($answer eq "yes") {	# we do want to remap the data

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

  }

  # the SUPPLEMENTARY_GFF directory is copied over whether or not it has been remapped
  $wormbase->run_command("cp -R ".$wormbase->build_data."/SUPPLEMENTARY_GFF ".$wormbase->chromosomes."/", $log);
   
}

#__ end remap_misc_dynamic __#

sub make_UTR {
  my ($log)=@_;
  
  $log->write_to("bsub commands . . . . \n\n");
  my $lsf = LSF::JobManager->new(-J => "make_UTRs", -o => "/dev/null");
  foreach (@{$wormbase->get_binned_chroms} ) {
    my $cmd = "make_UTR_GFF.pl -chromosome $_";
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
  
  #merge into single file.
  if($wormbase->assembly_type eq 'contig') {
    $wormbase->run_command("cat ".$wormbase->gff_splits."/UTR.gff* >> ".$wormbase->gff_splits."/allUTR.gff",$log);
    $wormbase->run_command("rm -f ".$wormbase->gff_splits."/UTR.gff*",$log);
    $wormbase->run_command("mv -f ".$wormbase->gff_splits."/allUTR.gff ".$wormbase->gff_splits."/UTR.gff",$log);
	
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
		     'CHROMOSOME_III'     => 3200000,
		     'CHROMOSOME_IV'      => 3500000,
		     'CHROMOSOME_V'       => 4200000,
		     'CHROMOSOME_X'       => 3600000,
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

  #inverted
  $wormbase->run_script("run_inverted.pl -all" , $log);
}


sub ontologies {
	$wormbase->run_script( "ONTOLOGY/parse_expr_pattern_new.pl", $log);
	$wormbase->run_script( "ONTOLOGY/parse_go_terms_new.pl -rnai -gene", $log);
	$wormbase->run_script( "ONTOLOGY/parse_phenotype_new.pl", $log);
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
  $wormbase->run_script( "update_web_gene_names.pl", $log);
}
