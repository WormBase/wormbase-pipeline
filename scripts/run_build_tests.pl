#!/usr/bin/perl

use strict;
use Wormbase;
use Wormtest;
use Log_files;
use Getopt::Long;

my ($debug, $test, $database, $species, $log_file, $prev_release);
my ($recent_citace, $primary_seq_dumps, $elegans_first, $build_contents, $dna_headers, $split_gffs,
    $tier2_contigs, $masked_files, $blat_files, $homology_loaded, $tier2_blat_contigs, $est_dat,
    $cache_size, $uniprot_ids, $vep, $dbxref, $final_gff, $species_merge);
GetOptions(
    'debug:s'            => \$debug,
    'test'               => \$test,
    'database:s'         => \$database,
    'species:s'          => \$species,
    'logfile:s'          => \$log_file,
    'prev_release:s'     => \$prev_release,
    'recent_citace'      => \$recent_citace, 
    'primary_seq_dumps'  => \$primary_seq_dumps,
    'elegans_first'      => \$elegans_first,
    'build_contents'     => \$build_contents,
    'dna_headers'        => \$dna_headers,
    'split_gffs:s'       => \$split_gffs,
    'tier2_contigs:s'    => \$tier2_contigs,
    'masked_files'       => \$masked_files,
    'blat_files'         => \$blat_files,
    'homology_loaded'    => \$homology_loaded,
    'est_dat'            => \$est_dat,
    'cache_size'         => \$cache_size,
    'uniprot_ids'        => \$uniprot_ids,
    'vep'                => \$vep,
    'dbxref'             => \$dbxref,
    'final_gff'          => \$final_gff,
    'species_merge'      => \$species_merge
    );

my $wormbase = Wormbase->new(
    -test     => $test,
    -debug    => $debug,
    -organism => $species,
    -autoace  => $database
    );

my $log = $log_file ? Log_files->make_log($log_file, $debug) : Log_files->make_build_log($wormbase);

my $wormtest = Wormtest->make_build_tester($wormbase, $log, $prev_release);

if ($recent_citace) {
    $wormtest->recent_citace_dump;
}
if ($primary_seq_dumps) {
    $wormtest->primary_seq_dumps_present;
}
if ($elegans_first) {
    $wormtest->elegans_loaded_first;
}
if ($build_contents) {
    $wormtest->build_folder_contents_present;
}
if ($dna_headers) {
    $wormtest->dna_files_have_headers;
}
if ($split_gffs) {
    $split_gffs = lc($split_gffs);
    die "Value passed to -split_gffs parameter must be one of 'init', 'blat', 'homol', or 'variation\n"
	unless $split_gffs eq 'init' or $split_gffs eq 'blat'
	or $split_gffs eq 'homol' or $split_gffs eq 'variation';
    $wormtest->split_gffs_present($split_gffs);
}
if ($tier2_contigs) {
    die "Value passed to -tier2_contigs parameter must be either 'init' or 'blat'\n"
	unless $split_gffs eq 'init' or $split_gffs eq 'blat';
    $wormtest->tier2_contigs_dumped($tier2_contigs);
}
if ($masked_files) {
    $wormtest->masked_files_present;
}
if ($blat_files) {
    $wormtest->blat_files_present;
}
if ($homology_loaded) {
    $wormtest->homology_data_loaded;
}
if ($est_dat) {
    $wormtest->create_est_dat_files_if_required;
}
if ($cache_size) {
    $wormtest->cache_size_sufficient;
}
if ($uniprot_ids) {
    $wormtest->uniprot_ids_in_wormpep;
}
if ($vep) {
    $wormtest->vep_output_present;
}
if ($dbxref) {
    $wormtest->dbxref_report_correctly_formatted;
}
if ($final_gff) {
    $wormtest->final_gff_dumps_present;
}
if ($species_merge) {
    $wormtest->species_merge_successful;
}

$log->mail;

exit(0);
