#!/usr/bin/perl

# A script to clean up empty files from FTP staging directory before release
# Unexpected empty files will be flagged and not deleted

use strict;
use warnings;
use Getopt::Long;
use File::Find;
use Const::Fast;
use Wormbase;
use Storable;

const my %EXPECTED_EMPTY => (
    'c_briggsae.PRJNA10731.WSXXX.ncRNA_transcripts.fa.gz' => 1,
    'c_brenneri.PRJNA20035.WSXXX.variations.vcf.gz' => 1,
    'c_brenneri.PRJNA20035.WSXXX.ncRNA_transcripts.fa.gz' => 1,
    'p_pacificus.PRJNA12644.WSXXX.variations.vcf.gz' => 1,
    'p_pacificus.PRJNA12644.WSXXX.ncRNA_transcripts.fa.gz' => 1,
    's_ratti.PRJEB125.WSXXX.variations.vcf.gz' => 1,
    's_ratti.PRJEB125.WSXXX.ncRNA_transcripts.fa.gz' => 1,
    's_ratti.PRJEB125.WSXXX.pseudogenic_transcripts.fa.gz' => 1,
    's_ratti.PRJEB125.WSXXX.uniprot_papers.txt.gz' => 1,
    'o_volvulus.PRJEB513.WSXXX.ncRNA_transcripts.fa.gz' => 1,
    'o_volvulus.PRJEB513.WSXXX.variations.vcf.gz' => 1,
    'o_volvulus.PRJEB513.WSXXX.uniprot_papers.txt.gz' => 1,
    'c_elegans.PRJNA275000.WSXXX.canonical_geneset.gtf.gz' => 1,
    'c_elegans.PRJNA275000.WSXXX.mRNA_transcripts.fa.gz' => 1,
    'c_elegans.PRJNA275000.WSXXX.protein.fa.gz' => 1,
    'c_japonica.PRJNA12591.WSXXX.ncRNA_transcripts.fa.gz' => 1,
    'c_japonica.PRJNA12591.WSXXX.variations.vcf.gz' => 1,
    'b_malayi.PRJNA10729.WSXXX.variations.vcf.gz' => 1,
    't_muris.PRJEB126.WSXXX.variations.vcf.gz' => 1,
    't_muris.PRJEB126.WSXXX.ncRNA_transcripts.fa.gz' => 1,
    't_muris.PRJEB126.WSXXX.pseudogenic_transcripts.fa.gz' => 1,
    't_muris.PRJEB126.WSXXX.uniprot_papers.txt.gz' => 1,
    't_muris.PRJEB126.WSXXX.RNASeq_controls_FPKM.dat' => 1,
    'c_remanei.PRJNA53967.WSXXX.ncRNA_transcripts.fa.gz' => 1,
    'c_remanei.PRJNA53967.WSXXX.variations.vcf.gz' => 1
    );

my ($debug, $store, $test, $test_dir, $ws_version, $wormbase);

GetOptions(
    "debug=s"     => \$debug,
    "store=s"     => \$store,
    "test"        => \$test,
    "testdir=s"   => \$test_dir, # directory containing test FTP data,
    "wbversion=s" => \$ws_version 
    );

if ( $store ) {
  $wormbase = retrieve( $store ) or croak("Can't restore wormbase from $store\n");
} else {
    $wormbase = Wormbase->new(
	-debug   => $debug,
	-test    => $test
	);
}

my $log = Log_files->make_build_log($wormbase);

if (not $ws_version) {
  $ws_version = $wormbase->get_wormbase_version();
}
my $ws_version_name = "WS${ws_version}";

my $target_dir = ($test_dir) 
    ? $test_dir
    : $wormbase->ftp_site . "/releases/.${ws_version_name}";

$log->write_to("Checking files in $target_dir\n");

find(\&cleanup, "$target_dir");

$log->mail;
exit(0);

sub cleanup {
    if (-f) {
	my $filepath = $File::Find::name;
	my $read_cmd = "cat $filepath |";
	if ($filepath =~ /\.gz$/) {
	    $read_cmd = "gunzip -c $filepath |";
	}
	open (IN, $read_cmd) or die $!;
	my $non_header_lines = 0;
	while(<IN>) {
	    next if $_ =~ /^\s*$/;
	    next if $_ =~ /^#/;
	    next if $_ =~ /^Data not available/;
	    $non_header_lines = 1;
	    last;
	}
	close(IN);
	if (!$non_header_lines) {
	    my ($filename) = $filepath =~ /\/([^\/]+)$/;
	    my $generic_filename = $filename;
	    $generic_filename =~ s/WS\d\d\d/WSXXX/g;
	    if (!exists $EXPECTED_EMPTY{$generic_filename}) {
		$log->error("ERROR: Unexpected empty file ($filepath) has not been removed, check and remove manually if necessary\n");
	    } else {
		$log->write_to("Removing empty file $filepath\n");
		unlink $filepath or $log->error("ERROR: Could not delete $filepath: $!\n");
	    }
	}
    }
}
