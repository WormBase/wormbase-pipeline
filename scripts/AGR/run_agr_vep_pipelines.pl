#!/bin/env perl
use strict;
use Getopt::Long;
use Compress::Zlib;
use Const::Fast;
use JSON;
use LWP::Simple;
use File::Path qw(make_path);
use Time::Piece;
use File::Slurp;
use Digest::MD5;
use Log_files;
use Wormbase;
use DateTime;
use Modules::WormSlurm;

const my $FMS_LATEST_PREFIX => 'https://fms.alliancegenome.org/api/datafile/by/';
const my $FMS_LATEST_SUFFIX => '?latest=true';
const my @FMS_DATATYPES => ('FASTA', 'GFF', 'HTVCF', 'MOD-GFF-BAM-KNOWN', 'MOD-GFF-BAM-MODEL', 'VARIATION');
const my %DATATYPE_EXTENSIONS => ('FASTA'             => 'fa',
				  'GFF'               => 'gff',
				  'HTVCF'             => 'vcf',
				  'MOD-GFF-BAM-KNOWN' => 'bam',
				  'MOD-GFF-BAM-MODEL' => 'bam',
				  'VARIATION'         => 'json',
    );
const my @CHECKSUM_SUFFIXES => ('FASTA.fa', 'GFF.gff', 'VCF.vcf', 'HTVCF.vcf', 'BAM.bam');
const my %ASSEMBLIES => ('GRCm39'    => 'MGI',
			 'R6'        => 'FB',
			 'R627'      => 'FB',
			 'mRatBN7.2' => 'RGD',
			 'SGDr64'    => 'SGD',
			 'WBcel235'  => 'WB',
			 'GRCz11'    => 'ZFIN',
			 'HUMAN'     => 'HUMAN'
    );
const my $BASE_DIR => $ENV{'AGR_VEP_BASE_DIR'} . '/' . $ENV{'AGR_RELEASE'} . '/' . $ENV{'DOWNLOAD_DATE'};
const my $HGNC_FILE_URL => 'http://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt';
const my @IDS_TO_MAP => ('symbol', 'entrez_id', 'ensembl_gene_id', 'vega_id', 'ucsc_id', 'refseq_accession', 'mgd_id', 'rgd_id', 'omim_id', 'mim_id', 'agr');
const my $CHECKSUMS_FILE => $ENV{'AGR_VEP_BASE_DIR'} . '/mod_file_checksums.txt';
const my $HUMAN_FILES_DIR => $ENV{'AGR_VEP_BASE_DIR'} . '/human_vep_input_files';
const my $MOUSE_FILES_DIR => $ENV{'AGR_VEP_BASE_DIR'} . '/mouse_vep_input_files';
const my $RESOURCES_DIR => $ENV{'AGR_BASE_DIR'} . '/resources';

const my %REFSEQ_CHROMOSOMES => (
    'FB'   => {
	'2L' => 'NT_033779.5',
	'2R' => 'NT_033778.4',
	'3L' => 'NT_037436.4',
	'3R' => 'NT_033777.3',
	'4'  => 'NC_004353.4',
	'X'  => 'NC_004354.4',
	'Y'  => 'NC_024512.1',
	'mitochondrion_genome' => 'NC_024511.2',
	'Unmapped_Scaffold_8_D1580_D1567' => 'NW_007931083.1',
	'211000022278279' => 'NW_007931104.1',
	'211000022278436' => 'NW_001845431.1',
	'211000022278449' => 'NW_001845819.1',
	'211000022278760' => 'NW_001846712.1',
	'211000022279165' => 'NW_001846812.1',
	'211000022279188' => 'NW_001845284.1',
	'211000022279264' => 'NW_001847227.1',
	'211000022279392' => 'NW_001846198.1',
	'211000022279681' => 'NW_001845031.1',
	'211000022280328' => 'NW_001844935.1',
	'211000022280341' => 'NW_001846187.1',
	'211000022280347' => 'NW_001845870.1',
	'211000022280481' => 'NW_001845220.1',
	'211000022280494' => 'NW_001845164.1',
	'211000022280703' => 'NW_001845199.1',
	'rDNA' => 'NW_007931121.1',
    },
    'MGI'  => {
	# GRCm38
	#'1'  => 'NC_000067.6',
	#'2'  => 'NC_000068.7',
	#'3'  => 'NC_000069.6',
	#'4'  => 'NC_000070.6',
	#'5'  => 'NC_000071.6',
	#'6'  => 'NC_000072.6',
	#'7'  => 'NC_000073.6',
	#'8'  => 'NC_000074.6',
	#'9'  => 'NC_000075.6',
	#'10' => 'NC_000076.6',
	#'11' => 'NC_000077.6',
	#'12' => 'NC_000078.6',
	#'13' => 'NC_000079.6',
	#'14' => 'NC_000080.6',
	#'15' => 'NC_000081.6',
	#'16' => 'NC_000082.6',
	#'17' => 'NC_000083.6',
	#'18' => 'NC_000084.6',
	#'19' => 'NC_000085.6',
	#'X'  => 'NC_000086.7',
	#'Y'  => 'NC_000087.7',
	#'MT' => 'NC_005089.1',
	# GRCm39
	'1'  => 'NC_000067.7',
	'2'  => 'NC_000068.8',
	'3'  => 'NC_000069.7',
	'4'  => 'NC_000070.7',
	'5'  => 'NC_000071.7',
	'6'  => 'NC_000072.7',
	'7'  => 'NC_000073.7',
	'8'  => 'NC_000074.7',
	'9'  => 'NC_000075.7',
	'10' => 'NC_000076.7',
	'11' => 'NC_000077.7',
	'12' => 'NC_000078.7',
	'13' => 'NC_000079.7',
	'14' => 'NC_000080.7',
	'15' => 'NC_000081.7',
	'16' => 'NC_000082.7',
	'17' => 'NC_000083.7',
	'18' => 'NC_000084.7',
	'19' => 'NC_000085.7',
	'X'  => 'NC_000086.8',
	'Y'  => 'NC_000087.8',
	'MT' => 'NC_005089.1',
    },
    'RGD'  => {
	# mRatBN7.2
	'1'  => 'NC_051336.1',
	'2'  => 'NC_051337.1',
	'3'  => 'NC_051338.1',
	'4'  => 'NC_051339.1',
	'5'  => 'NC_051340.1',
	'6'  => 'NC_051341.1',
	'7'  => 'NC_051342.1',
	'8'  => 'NC_051343.1',
	'9'  => 'NC_051344.1',
	'10' => 'NC_051345.1',
	'11' => 'NC_051346.1',
	'12' => 'NC_051347.1',
	'13' => 'NC_051348.1',
	'14' => 'NC_051349.1',
	'15' => 'NC_051350.1',
	'16' => 'NC_051351.1',
	'17' => 'NC_051352.1',
	'18' => 'NC_051353.1',
	'19' => 'NC_051354.1',
	'20' => 'NC_051355.1',
	'X'  => 'NC_051356.1',
	'Y'  => 'NC_051357.1',
	'MT' => 'NC_001665.2',
	# Rnor60	    
	#'1'  => 'NC_005100.4',
	#'2'  => 'NC_005101.4',
	#'3'  => 'NC_005102.4',
	#'4'  => 'NC_005103.4',
	#'5'  => 'NC_005104.4',
	#'6'  => 'NC_005105.4',
	#'7'  => 'NC_005106.4',
	#'8'  => 'NC_005107.4',
	#'9'  => 'NC_005108.4',
	#'10' => 'NC_005109.4',
	#'11' => 'NC_005110.4',
	#'12' => 'NC_005111.4',
	#'13' => 'NC_005112.4',
	#'14' => 'NC_005113.4',
	#'15' => 'NC_005114.4',
	#'16' => 'NC_005115.4',
	#'17' => 'NC_005116.4',
	#'18' => 'NC_005117.4',
	#'19' => 'NC_005118.4',
	#'20' => 'NC_005119.4',
	#'X'  => 'NC_005120.4',
	#'Y'  => 'NC_024475.1',
	#'MT' => 'NC_001665.2',
    },
    'SGD'  => {
	'chrI'    => 'NC_001133.9',
	'chrII'   => 'NC_001134.8',
	'chrIII'  => 'NC_001135.5',
	'chrIV'   => 'NC_001136.10',
	'chrV'    => 'NC_001137.3',
	'chrVI'   => 'NC_001138.5',
	'chrVII'  => 'NC_001139.9',
	'chrVIII' => 'NC_001140.6',
	'chrIX'   => 'NC_001141.2',
	'chrX'    => 'NC_001142.9',
	'chrXI'   => 'NC_001143.9',
	'chrXII'  => 'NC_001144.5',
	'chrXIII' => 'NC_001145.3',
	'chrXIV'  => 'NC_001146.8',
	'chrXV'   => 'NC_001147.6',
	'chrXVI'  => 'NC_001148.4',
	'chrmt'   => 'NC_001224.1',
    },
    'WB'   => {
	'I'     => 'NC_003279.8',
	'II'    => 'NC_003280.10',
	'III'   => 'NC_003281.10',
	'IV'    => 'NC_003282.8',
	'V'     => 'NC_003283.11',
	'X'     => 'NC_003284.9',
	'MtDNA' => 'NC_001328.1',
    },
    'ZFIN' => {
	'1'  => 'NC_007112.7',
	'2'  => 'NC_007113.7',
	'3'  => 'NC_007114.7',
	'4'  => 'NC_007115.7',
	'5'  => 'NC_007116.7',
	'6'  => 'NC_007117.7',
	'7'  => 'NC_007118.7',
	'8'  => 'NC_007119.7',
	'9'  => 'NC_007120.7',
	'10' => 'NC_007121.7',
	'11' => 'NC_007122.7',
	'12' => 'NC_007123.7',
	'13' => 'NC_007124.7',
	'14' => 'NC_007125.7',
	'15' => 'NC_007126.7',
	'16' => 'NC_007127.7',
	'17' => 'NC_007128.7',
	'18' => 'NC_007129.7',
	'19' => 'NC_007130.7',
	'20' => 'NC_007131.7',
	'21' => 'NC_007132.7',
	'22' => 'NC_007133.7',
	'23' => 'NC_007134.7',
	'24' => 'NC_007135.7',
	'25' => 'NC_007136.7',
	'MT' => 'NC_002333.2',
    },
    'HUMAN' => {
	'1'  => 'NC_000001.11',
	'2'  => 'NC_000002.12',
	'3'  => 'NC_000003.12',
	'4'  => 'NC_000004.12',
	'5'  => 'NC_000005.10',
	'6'  => 'NC_000006.12',
	'7'  => 'NC_000007.14',
	'8'  => 'NC_000008.11',
	'9'  => 'NC_000009.12',
	'10' => 'NC_000010.11',
	'11' => 'NC_000011.10',
	'12' => 'NC_000012.12',
	'13' => 'NC_000013.11',
	'14' => 'NC_000014.9',
	'15' => 'NC_000015.10',
	'16' => 'NC_000016.10',
	'17' => 'NC_000017.11',
	'18' => 'NC_000018.10',
	'19' => 'NC_000019.10',
	'20' => 'NC_000020.11',
	'21' => 'NC_000021.9',
	'22' => 'NC_000022.11',
	'X'  => 'NC_000023.11',
	'Y'  => 'NC_000024.10',
	'MT' => 'NC_012920.1'
    },
    );

const my %BAM_REQUIRED => (
    'FB'    => 0,
    'MGI'   => 1,
    'RGD'   => 1,
    'SGD'   => 0,
    'WB'    => 0,
    'ZFIN'  => 1,
    'HUMAN' => 1,
    );

my ($url, $test, $logfile, $debug, $password, $cleanup, $nocheck, $overwrite, $external_human_gff, $external_mouse_fasta, $help);

my $stages = '1,2,3,4,5';
my $mods_string = 'FB,MGI,RGD,SGD,WB,ZFIN,HUMAN';

GetOptions(
    "mods|m=s"              => \$mods_string,
    "stages|s=s"            => \$stages,
    "password|p=s"          => \$password,
    "logfile|l=s"           => \$logfile,
    "debug|d=s"             => \$debug,
    "url|u=s"               => \$url,
    "test|t"                => \$test,
    "cleanup|c"             => \$cleanup,
    "overwrite|o"           => \$overwrite,
    "nocheck|n"             => \$nocheck,
    "external_human_gff|e"  => \$external_human_gff,
    "external_mouse_fasta|m" => \$external_mouse_fasta,
    "help|h"                => \$help,
    ) or print_usage();

print_usage() if $help;

make_path($BASE_DIR) unless -d $BASE_DIR;
    
my $start_time = DateTime->now->strftime('%Y%m%d%H%M%S');

my @mods = split(',', $mods_string);
$logfile = "${BASE_DIR}/submission.${start_time}.log" if !$logfile;;

my $log = Log_files->make_log($logfile, $debug);
download_from_agr(\@mods, $start_time, $url, $overwrite, $external_human_gff, $external_mouse_fasta, $log) if $stages =~ /1/;

for my $mod (@mods) {
    my ($checksums, $run_stages) = check_for_new_data($mod, $nocheck, $log);
    chdir "$BASE_DIR/$mod";
    if ($stages =~ /2/) {
	if ($run_stages->{2}) {
	    process_input_files($mod, $external_human_gff, $log);
	}
	else {
	    $log->write_to("Skipping processing of input files for $mod " .
			   'as input files unchanged and no further analyses ' .
			   "are being carried out\n");
	}
    }
    if ($stages =~ /3/) {
	if ($run_stages->{3}) {
	    calculate_pathogenicity_predictions($mod, $password, $test, $log);
	}
	else {
	    $log->write_to('Skipping pathogenicity prediction calculations ' .
			   "for $mod as input files unchanged\n\n");
	}
    }
    if ($stages =~ /4/) {
	if ($run_stages->{4}) {
	    run_vep_on_phenotypic_variations($mod, $password, $test, $log);
	    update_checksums($mod, 'VCF.vcf', $checksums, $log) if !$test;
	}
	else {
	    $log->write_to('Skipping VEP analysis of phenotypic variants ' .
			   "for $mod as input files unchanged\n\n");
	}
    }
    if ($stages =~ /5/) {
	if ($run_stages->{5}) {
	    run_vep_on_htp_variations($mod, $password, $test, $log);
	    update_checksums($mod, 'HTVCF.vcf', $checksums, $log) if !$test;
	}
	else {
	    $log->write_to('Skipping VEP analysis of HTP variants ' .
			   "for $mod as input files unchanged\n\n");
	}
    }
    if (!$test and $stages =~ /3/ and $stages =~ /4/ and $stages =~ /5/) {
	update_checksums($mod, 'FASTA.fa', $checksums, $log);
	update_checksums($mod, 'GFF.gff', $checksums, $log);
	update_checksums($mod, 'BAM.bam', $checksums, $log) if $BAM_REQUIRED{$mod};
    }
    cleanup_intermediate_files($mod, $log) if $cleanup;
}

$log->mail;
exit(0);


sub check_for_new_data {
    my ($mod, $nocheck, $log) = @_;

    my $old_checksums = get_old_checksums($mod, $log);
    my $new_checksums = get_new_checksums($mod, $log);
    
    my %run_stages;
    if ($nocheck) {
	%run_stages = map {$_ => 1} (2 .. 5);
	return ($new_checksums, \%run_stages) if $nocheck;
    }

    %run_stages = map {$_ => 0} (2 .. 5);

    if (!exists $old_checksums->{"${mod}_VCF.vcf"} or
	$old_checksums->{"${mod}_VCF.vcf"} ne $new_checksums->{"${mod}_VCF.vcf"}) {
	$run_stages{2} = 1;
	$run_stages{4} = 1 ;
    }
    if (!exists $old_checksums->{"${mod}_HTVCF.vcf"} or
	$old_checksums->{"${mod}_HTVCF.vcf"} ne $new_checksums->{"${mod}_HTVCF.vcf"}) {
	$run_stages{2} = 1;
	$run_stages{5} = 1;
    }
    if (!exists $old_checksums->{"${mod}_FASTA.fa"} or
	($BAM_REQUIRED{$mod} and !exists $old_checksums->{"${mod}_BAM.bam"}) or
	!exists $old_checksums->{"${mod}_GFF.gff"} or
	$old_checksums->{"${mod}_FASTA.fa"} ne $new_checksums->{"${mod}_FASTA.fa"} or
	($BAM_REQUIRED{$mod} and $old_checksums->{"${mod}_BAM.bam"} ne $new_checksums->{"${mod}_BAM.bam"}) or
	$old_checksums->{"${mod}_GFF.gff"} ne $new_checksums->{"${mod}_GFF.gff"}) {
	%run_stages = map {$_ => 1} (2 .. 5);
    }
    
    return ($new_checksums, \%run_stages);
}


sub get_new_checksums {
    my ($mod, $log) = @_;

    my %checksums;
    for my $suffix (@CHECKSUM_SUFFIXES) {
	my $file = "$BASE_DIR/$mod/${mod}_$suffix";
	next unless -e $file;
	open (my $fh, '<', $file) or $log->log_and_die("Cannot open $file for reading\n");
	my $md5 = Digest::MD5->new;
	$md5->addfile($fh);
	$checksums{"${mod}_${suffix}"} = $md5->hexdigest;
	close ($fh);
    }

    return \%checksums;
}


sub get_old_checksums {
    my ($mod, $log) = @_;
    my %checksums;
    run_system_cmd("touch $CHECKSUMS_FILE", "Comparing checksums of new and old $mod files", $log);
    open (CHECKSUM, '<', $CHECKSUMS_FILE) or $log->log_and_die("Couldn't open $CHECKSUMS_FILE for reading\n");
    while (<CHECKSUM>) {
	chomp;
	next unless $_ =~ /^${mod}_/;
	my ($file, $checksum) = /^([^\s]+)\s(.+)$/;
	$checksums{$file} = $checksum;
    }
    close (CHECKSUM);

    return \%checksums;
}


sub update_checksums {
    my ($mod, $suffix, $checksums, $log) = @_;

    my $file_to_update = $mod . '_' . $suffix;
    open (IN, '<', $CHECKSUMS_FILE) or $log->log_and_die("Cannot open $CHECKSUMS_FILE for reading\n");
    open (OUT, '>', $CHECKSUMS_FILE . '.tmp') or $log->log_and_die("Cannot open $CHECKSUMS_FILE.tmp for writing\n");
    while (<IN>) {
	my ($file, $checksum) = split("\s", $_);
	next if $file eq $file_to_update;
	print OUT $_;
    }
    print OUT $file_to_update . ' ' . $checksums->{$file_to_update} . "\n";
    close (IN);
    close (OUT);
    run_system_cmd("mv ${CHECKSUMS_FILE}.tmp ${CHECKSUMS_FILE}", "Updating ${mod}_${suffix} checksum", $log);

    return;
}


sub download_from_agr {
    my ($mods, $start_time, $url, $overwrite, $external_human_gff, $external_mouse_fasta, $log) = @_;
    
    my $download_urls = defined $url ? get_urls_from_snapshot($url) : get_latest_urls();
    
    my $input_files_file = "${BASE_DIR}/VEP_input_files.txt";
    
    open (FILES, '>>', $input_files_file) or $log->lod_and_die("Couldn't open $input_files_file to append data\n");
    for my $mod (@$mods) {
	my $time = localtime();
	print FILES "${mod}: $time\n";
	my $mod_dir = "$BASE_DIR/$mod";
	make_path($mod_dir) unless -d $mod_dir;
	chdir $mod_dir;
	for my $datatype (keys %{$download_urls->{$mod}}){
	    next if $external_human_gff and $mod eq 'HUMAN' and $datatype eq 'GFF';
	    if ($external_mouse_fasta and $mod eq 'MGI' and $datatype eq 'FASTA') {
		run_system_cmd("gunzip -c ${MOUSE_FILES_DIR}/Mus_musculus.GRCm39.dna.toplevel.fa.gz > MGI_FASTA.fa", "Unzipping local MGI FASTA", $log);
		next;
	    }
	    my $extension = $DATATYPE_EXTENSIONS{$datatype};
	    if (-e "${mod}_${datatype}.${extension}" and !$overwrite) {
		$log->write_to("Using previously downloaded $mod $datatype\n");
		next;
	    }
	    my ($filename) = $download_urls->{$mod}{$datatype} =~ /\/([^\/]+)$/;
	    print FILES "\t${datatype}: ${filename}\n";
	    run_system_cmd('curl -O ' . $download_urls->{$mod}{$datatype}, "Downloading $mod $datatype file", $log);
	    $filename = check_if_actually_compressed($filename, $log) if $filename !~ /\.gz/; # temporary hack to get around gzipped files in FMS without .gz extension
	    run_system_cmd("gunzip $filename", "Decompressing $filename", $log) if $filename =~ /\.gz$/; # if clause only required in interim while some FMS files not gzipped
	    $filename =~ s/\.gz$//;
	    run_system_cmd("mv $filename ${mod}_${datatype}.${extension}", "Renaming $filename", $log);
	}
	fix_rgd_headers($log) if $mod eq 'RGD'; # Temporary hack

	if ($BAM_REQUIRED{$mod}) {
	    merge_bam_files($mod, $log);
	}
	else {
	    run_system_cmd('cp ' . $RESOURCES_DIR . '/dummy.bam ' . $mod . '_BAM.bam', "Copying dummy BAM file", $log);
	}
	run_slurm_job("samtools index ${mod}_BAM.bam", "Indexing $mod BAM file", $log, '00:30:00', 4, '/dev/null', '/dev/null');
	
	unlink "${mod}_FASTA.fa.fai" if -e "${mod}_FASTA.fa.fai";
	run_slurm_job('python3 ' . $ENV{'CVS_DIR'} . "/AGR/agr_variations_json2vcf.py -j ${mod}_VARIATION.json -m $mod -g ${mod}_GFF.gff " .
		      "-f ${mod}_FASTA.fa -o ${mod}_VCF.vcf", "Converting $mod phenotypic variants JSON to VCF", $log, '00:30:00', 4, '/dev/null', '/dev/null') if -e "${mod}_VARIATION.json";
	if ($mod eq 'HUMAN') {
	    # May need to reimplement below once we move back to full set of human variants and not just RGD-submitted ClinVar variants
	    
	    # Copy across files not on FMS from local directory
	    #my @files_to_copy = ('HUMAN_HTVCF.vcf');
	    #push @files_to_copy, 'HUMAN_GFF.gff' if $external_human_gff;
	    #for my $file (@files_to_copy) {
		#run_system_cmd("cp ${HUMAN_FILES_DIR}/${file} $file", "Copying local $file file to working directory", $log);
	    #}   
	}
	elsif ($mod eq 'MGI') {

	}
    
    }
    close (FILES);

    
    return;
}


sub fix_rgd_headers {
    # Temporary hack
    my $log = shift;
    open (IN, '<', 'RGD_HTVCF.vcf') or $log->log_and_die("Couldn't open RGD_HTVCF.vcf for reading\n");
    open (OUT, '>', 'RGD_HTVCF.vcf.tmp') or $log->log_and_die("Couldn't open RGD_HTVCF.vcf.tmp for writing\n");
    while (<IN>) {
	if ($_ =~ /^##fileformat=VCF4.2/) {
	    print OUT "##fileformat=VCFv4.2\n";
	}
	elsif ($_ =~ /^##INFO=<ID=VT,Number=0,Type=Integer/) {
	    print OUT "##INFO=<ID=VT,Number=1,Type=String,Description=\"Variant type: SNP, INS, and DEL\">\n";
	}
	elsif ($_ =~ /^##FORMAT=<ID=DP,Number=2,Type=String/) {
	    print OUT "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n";
	}
	else {
	    print OUT $_;
	}
    }
    close (IN);
    close (OUT);
    run_system_cmd('mv RGD_HTVCF.vcf.tmp RGD_HTVCF.vcf', 'Replacing RGD HTVCF file with fixed header version', $log);

    return;
}


sub get_latest_urls {
    my %download_urls;
    for my $datatype (@FMS_DATATYPES) {
	my $content = get($FMS_LATEST_PREFIX . $datatype . $FMS_LATEST_SUFFIX);
	my $latest = decode_json($content);
	for my $entry (@$latest) {
	    my $mod = $entry->{dataSubType}->{name};
	    if ($datatype eq 'FASTA' or $datatype eq 'VCF') {
		next unless exists $ASSEMBLIES{$mod};
		$mod = $ASSEMBLIES{$mod};
	    }
	    $download_urls{$mod}{$datatype} = $entry->{s3Url};
	}
    }
    
    return \%download_urls;
}


sub get_urls_from_snapshot {
    my $snapshot_url = shift;
    
    my %download_urls;
    my $content = get($snapshot_url);
    my $snapshot = decode_json($content);
    for my $entry (@{$snapshot->{snapShot}{dataFiles}}) {
	my $datatype = $entry->{dataType}{name};
	next unless exists $DATATYPE_EXTENSIONS{$datatype};
	my $mod = $entry->{dataSubType}{name};
	if ($datatype eq 'FASTA' or $datatype eq 'VCF' ) {
	    next unless exists $ASSEMBLIES{$mod};
	    $mod = $ASSEMBLIES{$mod};
	}
	$download_urls{$mod}{$datatype} = $entry->{s3Url};
    }
    
    return \%download_urls;
}


sub check_if_actually_compressed {
    my ($filename, $log) = @_;
    
    my $is_gzipped = 0;
    open (FILE, "file $filename |") or $log->log_and_die($!);
    while (<FILE>) {
	chomp;
	$is_gzipped = 1 if $_ =~ /gzip compressed data/ || $_ =~ /Zip archive data/;
    }
    close (FILE);
    
    if ($is_gzipped) {
	run_system_cmd("mv $filename $filename.gz", "Adding .gz extension to $filename", $log);
	return $filename . '.gz';
    }
    
    return $filename;
}


sub process_input_files {
    my ($mod, $external_human_gff, $log) = @_;
    
    cleanup_intermediate_files($mod, $log);
    
    sort_vcf_files($mod, $log); # unless $mod eq 'HUMAN'; # No need for human as using local (sorted) file Update 25092024: (not any more - may revert in future)

    my $chr_map;
    check_chromosome_map($mod, $log);
    convert_fasta_headers($mod, $log);
    convert_vcf_chromosomes($mod, 'VCF', $log);
  
    munge_gff($mod, $external_human_gff, $log);
    run_slurm_job("bgzip -c ${mod}_FASTA.refseq.fa", "Compressing $mod FASTA", $log, '01:00:00', 4, "${mod}_FASTA.refseq.fa.gz", '/dev/null');
    run_slurm_job("sort -k1,1 -k4,4n -k5,5n -t\$'\\t' ${mod}_GFF.refseq.gff", "Sorting $mod GFF", $log, '00:30:00', 4, "${mod}_GFF.refseq.sorted.gff", '/dev/null');
    run_system_cmd("mv ${mod}_GFF.refseq.sorted.gff ${mod}_GFF.refseq.gff", "Renaming sorted GFF", $log);
    run_slurm_job("bgzip -c ${mod}_GFF.refseq.gff", "Compressing sorted GFF", $log, '01:00:00', 4, "${mod}_GFF.refseq.gff.gz", '/dev/null');
    run_slurm_job("tabix -p gff ${mod}_GFF.refseq.gff.gz", "Indexing $mod GFF", $log, '01:00:00', 4, '/dev/null', '/dev/null');
    
    return;
}


sub calculate_pathogenicity_predictions {
    my ($mod, $password, $test, $log) = @_;
    
    backup_pathogenicity_prediction_db($mod, $password, $log);
    
    my $init_cmd = "init_pipeline.pl VepProteinFunction::VepProteinFunction_conf -mod $mod" .
	" -agr_fasta ${mod}_FASTA.refseq.fa -agr_gff ${mod}_GFF.refseq.gff -agr_bam ${mod}_BAM.bam" . 
	' -hive_root_dir ' . $ENV{'HIVE_ROOT_DIR'} . ' -pipeline_base_dir ' . $ENV{'PATH_PRED_WORKING_DIR'} .
	' -pipeline_host ' . $ENV{'WORM_DBHOST'} . ' -pipeline_user ' . $ENV{'WORM_DBUSER'} .
	' -pipeline_port ' . $ENV{'WORM_DBPORT'} .
	' -sift_dir ' . $ENV{'SIFT_DIR'} . ' -pph_dir ' . $ENV{'PPH_DIR'} . ' -pph_conf_dir ' . $ENV{'PPH_CONF_DIR'} . '/' . $mod .
	' -ncbi_dir ' . $ENV{'NCBI_DIR'} . ' -blastdb ' . $ENV{'BLAST_DB'} . ' -pph_blast_db ' . $ENV{'PPH_BLAST_DB'} .
	' -uniprot_dbs ' . $ENV{'UNIPROT_DBS'} . ' -password ' . $password;
    
    run_system_cmd($init_cmd, "Initialising $mod pathogenicity prediction eHive pipeline", $log);
    
    my $ehive_url = 'mysql://' . $ENV{'WORM_DBUSER'} . ':' . $password . '@' . $ENV{'WORM_DBHOST'} . ':' . 
	$ENV{'WORM_DBPORT'} . '/agr_pathogenicity_predictions_' . lc($mod) . '_ehive';
    $ENV{EHIVE_URL} = $ehive_url;
    run_system_cmd("beekeeper.pl -url $ehive_url -loop", "Running $mod pathogenicity prediction eHive pipeline", $log);
    
    return;
}


sub run_vep_on_phenotypic_variations {
    my ($mod, $password, $test, $log) = @_;
    
    return unless -e "${mod}_VCF.refseq.vcf";
    
    my $base_vep_cmd = "vep -i ${mod}_VCF.refseq.vcf -gff ${mod}_GFF.refseq.gff.gz -fasta ${mod}_FASTA.refseq.fa.gz --force_overwrite " .
	"-hgvsg -hgvs -shift_hgvs=0 --symbol --distance 0 --plugin ProtFuncSeq,mod=$mod,pass=$password " .
	"--remove_hgvsp_version --safe --bam ${mod}_BAM.bam";
    my $gl_vep_cmd = $base_vep_cmd . " --per_gene --output_file ${mod}_VEPGENE.txt";
    my $tl_vep_cmd = $base_vep_cmd . " --output_file ${mod}_VEPTRANSCRIPT.txt";
    
    run_slurm_job($gl_vep_cmd, "Running VEP for $mod phenotypic variants at gene level", $log, '00:30:00', 4, '/dev/null', '/dev/null');
    run_slurm_job($tl_vep_cmd, "Running VEP for $mod phenotypic variants at transcript level", $log, '00:30:00', 4, '/dev/null', '/dev/null');
    
    my $reverse_map = get_reverse_chromosome_map($mod);

    for my $level ('GENE', 'TRANSCRIPT') {
	open (IN, '<', "${mod}_VEP${level}.txt") or die $log->log_and_die("Cannot open ${mod}_VEP${level}.txt for reading\n");
	open (OUT, '>', "${mod}_VEP${level}.txt.tmp") or die $log->log_and_die("Cannot open ${mod}_VEP${level}.txt.tmp for writing\n");
	while (<IN>) {
	    chomp;
	    if ($_ =~ /^#/) {
		print OUT $_ . "\n";
	    }
	    else {
		my @columns = split("\t", $_);
		my ($chr, $pos) = split(':', $columns[1]);
		

		my ($before, $hgvsg, $after) = $columns[13] =~ /(.*HGVSg=)([^;]+)(.*)/;
		if (defined $hgvsg) {
		    $log->write_to('WARNING: HGVSg in input VCF (' . $columns[0] . ") doesn't match VEP generated HGVSg ($hgvsg)\n")
			unless $columns[0] eq $hgvsg;
		    # HGVSg in extras column needs to have chromosome name not RefSeq chr ID
		    my @hgvsg_parts = split(':', $hgvsg);
		    $hgvsg_parts[0] = $reverse_map->{$chr};
		    my $old_hgvsg = join(':', @hgvsg_parts);
		    $columns[13] = $before . $old_hgvsg . $after;
		}

		#Position should have chromosome name not RefSeq chr ID
		$columns[1] = $reverse_map->{$chr} . ':' . $pos;
		print OUT join("\t", @columns) . "\n";
	    }
	}
	close (IN);
	close (OUT);
	run_system_cmd("mv ${mod}_VEP${level}.txt.tmp ${mod}_VEP${level}.txt",
		       "Replacing $mod VEP${level} file with original chromosome ID version", $log);

	my $lclevel = lc($level);
	run_slurm_job("gzip -f -9 ${mod}_VEP${level}.txt", 'Compressing ${lclevel}-level VEP results', $log, '00:05:00', 1, '/dev/null', '/dev/null');
	submit_data($mod, 'VEP' . $level, $mod . '_VEP' . $level . '.txt.gz', $log) unless $test;
    }
    
    return;
}


sub run_vep_on_htp_variations{
    my ($mod, $password, $test, $log) = @_;

    my $init_cmd = "init_pipeline.pl ModVep::ModVep_conf -mod $mod -vcf ${mod}_HTVCF.vcf -gff ${mod}_GFF.refseq.gff.gz" .
	" -fasta ${mod}_FASTA.refseq.fa.gz -bam ${mod}_BAM.bam -hive_root_dir " . $ENV{'HIVE_ROOT_DIR'} . ' -pipeline_base_dir ' .
	$ENV{'HTP_VEP_WORKING_DIR'} . ' -pipeline_host ' . $ENV{'WORM_DBHOST'} . ' -pipeline_user ' . $ENV{'WORM_DBUSER'} .
	' -pipeline_port ' . $ENV{'WORM_DBPORT'} . ' -vep_dir ' . $ENV{'VEP_DIR'} . 
	" -debug_mode 0 -password $password";
    
    run_system_cmd($init_cmd, "Initialising $mod HTP variants VEP eHive pipeline: $init_cmd", $log);
    my $ehive_url = 'mysql://' . $ENV{'WORM_DBUSER'} . ':' . $password . '@' . $ENV{'WORM_DBHOST'} . ':' . 
	$ENV{'WORM_DBPORT'} . '/agr_htp_' . lc($mod) . '_vep_ehive';
    $ENV{EHIVE_URL} = $ehive_url;
 
    run_system_cmd("beekeeper.pl -url $ehive_url -loop", "Running $mod HTP variations VEP eHive pipeline", $log);
    run_slurm_job("cp " . $ENV{'HTP_VEP_WORKING_DIR'} . "/${mod}_vep/${mod}.vep.vcf.gz .", "Copying $mod combined HTP variations VEP output", $log, '01:00:00', 1, '/dev/null', '/dev/null');
    if ($mod eq 'RGD') {
	submit_data($mod, 'HTPOSTVEPVCF', "${mod}.vep.vcf.gz", $log) unless $test;
    }
    run_system_cmd("mkdir HTPVEP", "Creating folder for $mod HTP VEP output", $log);
    run_slurm_job("mv " . $ENV{'HTP_VEP_WORKING_DIR'} . "/${mod}_vep/${mod}.* HTPVEP/",
		   "Moving $mod HTP variations VEP output", $log,'01:05:00', 1, '/dev/null', '/dev/null');
    
    return;
}


sub submit_data {
    my ($mod, $fms_datatype, $file, $log) = @_;

    my $cmd = 'curl -H "Authorization: Bearer ' . $ENV{'TOKEN'} . '" -X POST ' .
	'"https://fms.alliancegenome.org/api/data/submit" -F "' . $ENV{'AGR_RELEASE'} . '_' .
	$fms_datatype . '_' . $mod . '=@' . $file . '"';

    my $response_json = `$cmd`;
    my $response = decode_json($response_json);
    if ($response->{status} eq 'failed') {
	$log->error("Upload of $mod $fms_datatype failed:\n$response_json\n\n");
    } 
    else {
	$log->write_to("Upload of $mod ${fms_datatype} succeeded\n\n");
    }

    return;
}    


sub merge_bam_files {
    my ($mod, $log) = @_;
    
    if (-e "${mod}_MOD-GFF-BAM-KNOWN.bam") {
	if (-e "${mod}_MOD-GFF-BAM-MODEL.bam") {
	    run_slurm_job("samtools merge -f ${mod}_BAM.bam ${mod}_MOD-GFF-BAM-KNOWN.bam ${mod}_MOD-GFF-BAM-MODEL.bam",
			   "Merging $mod BAM files", $log, '00:30:00', 4, '/dev/null', '/dev/null');
	    run_system_cmd("rm ${mod}_MOD-GFF-BAM-KNOWN.bam", "Deleting $mod unmerged known transcripts BAM file", $log);
	    run_system_cmd("rm ${mod}_MOD-GFF-BAM-MODEL.bam", "Deleting $mod unmerged model transcripts BAM file", $log);
	}
	else{
	    run_system_cmd("mv ${mod}_MOD-GFF-BAM-KNOWN.bam ${mod}_BAM.bam", "Renaming $mod MOD-GFF-BAM-KNOWN file", $log);
	}
    }
    else{
	$log->log_and_die("$mod BAM files could not be found\n") unless -e "${mod}_MOD-GFF-BAM-MODEL.bam";
	run_system_cmd("mv ${mod}_MOD-GFF-BAM-MODEL.bam ${mod}_BAM.bam", "Renaming $mod MOD-GFF-BAM-MODEL file", $log); 
    }
    
    run_slurm_job("samtools sort -o ${mod}_BAM.sorted.bam -T tmp ${mod}_BAM.bam", "Sorting $mod BAM file", $log, '00:30:00', 4, '/dev/null', '/dev/null');
    run_system_cmd("mv ${mod}_BAM.sorted.bam ${mod}_BAM.bam", "Replacing $mod BAM file with sorted version", $log);
    
    return;
}


sub sort_vcf_files {
    my ($mod, $log) = @_;
    
    for my $datatype ('VCF', 'HTVCF') { # not required for phenotypic variants VCF as JSON conversion script sorts output
	next unless -e "${mod}_${datatype}.vcf";
	run_slurm_job("vcf-sort ${mod}_${datatype}.vcf", "Sorting $mod $datatype file", $log, '01:00:00', 4, "${mod}_${datatype}.sorted.vcf", '/dev/null');
	run_system_cmd("mv ${mod}_${datatype}.sorted.vcf ${mod}_${datatype}.vcf",
		       "Replacing unsorted $mod $datatype file with sorted version", $log);
    }
    
    return;
}


sub remove_mirna_primary_transcripts {
    my $log = shift;
    my (%mirna_parents, %parents);

    $log->write_to("Getting RefSeq HUMAN miRNA details from GFF\n");

    # Do first pass to get parents
    open (GFF, "grep -v '^#' HUMAN_GFF.gff |") or $log->log_and_die("Could not open HUMAN_GFF.gff for reading\n");
    while (<GFF>) {
	my $line = $_;
	next unless $line =~ /BestRefSeq/;
	chomp $line;
	
	my @columns = split("\t", $line);
	my %attr = split(/[=;]/, $columns[8]);
	$parents{$attr{'ID'}} = $attr{'Parent'} if exists $attr{'Parent'};
	if ($columns[2] eq 'miRNA') {
	    $mirna_parents{$attr{'Parent'}} = 1 ;
	}
    }
    close (GFF);

    $log->write_to("Removing miRNA primary transcripts from HUMAN GFF\n");

    open (IN, '< HUMAN_GFF.gff') or $log->log_and_die("Could not open HUMAN_GFF.gff for reading\n");
    open (OUT, '> HUMAN_GFF.tmp.gff') or $log->log_and_die("Could not open HUMAN_GFF.tmp.gff for writing\n");
    while (<IN>) {
	my $line = $_;
	if ($line =~ /BestRefSeq/) {
	    chomp $line;
	    my @columns = split("\t", $line);
	    my %attr = split(/[=;]/, $columns[8]);
	    next if exists $mirna_parents{$attr{'ID'}}; # Skip primary_transcript miRNA parents
	    if ($columns[2] eq 'miRNA') {
		$attr{'Parent'} = $parents{$attr{'Parent'}}; # Replace parent of miRNA is miRNA gene ID
		my @key_values;
		for my $key (keys %attr) {
		    push @key_values, $key . '=' . $attr{$key};
		}
		$columns[8] = join(';', @key_values);
		print OUT join("\t", @columns) . "\n";
	    }
	    else {
		print OUT $line . "\n";
	    }
	}
	else {
	    print OUT $line;
	}
    }
    close (IN);
    close (OUT);

    run_system_cmd("mv HUMAN_GFF.tmp.gff HUMAN_GFF.gff", "Replacing GFF with version without pre-miRNA lines", $log);

    return;
}

sub remove_fasta_portion {
    my $mod = shift;

    my $in_file = "${mod}_GFF.gff";
    my $out_file = "${mod}_GFF.gff.tmp";
    open(IN, "<$in_file");
    open(OUT, ">$out_file");
    while(<IN>) {
        my $line = $_;
        last if $line =~ /^#+FASTA/;
        print OUT $line;
    }
    close (IN);
    close (OUT);

    system("mv $out_file $in_file");

    return;
}

sub munge_gff {
    my ($mod, $external_human_gff, $log) = @_;
    
    run_system_cmd("rm ${mod}_GFF.refseq.gff", "Deleting old munged GFF file", $log) if -e "${mod}_GFF.refseq.gff";

    remove_fasta_portion($mod) if $mod eq 'SGD';
    
    my $reverse_map = get_reverse_chromosome_map($mod);

    my $hgnc_id_map;
    if ($mod eq 'HUMAN' and $external_human_gff) {
	remove_mirna_primary_transcripts($log);
	$hgnc_id_map = get_hgnc_id_map($log);
    }

    $log->write_to("Munging $mod GFF\n");
    open(IN, "grep -v '^#' ${mod}_GFF.gff |") or $log->log_and_die("Could not open ${mod}_GFF.gff for reading\n");
    open(OUT, "> ${mod}_GFF.refseq.gff") or $log->log_and_die("Could not open ${mod}_GFF.refseq.gff for writing\n");
    while (<IN>) {
	my $line = $_;
	chomp $line;
	next if $line eq '';
	
	my @columns = split("\t", $line);
	if (exists $REFSEQ_CHROMOSOMES{$mod}{$columns[0]}) {
	    $columns[0] = $REFSEQ_CHROMOSOMES{$mod}{$columns[0]};
	    $line = join("\t", @columns);
	}
	else {
	    $log->log_and_die("Could not map $mod chromosome in GFF " . $columns[0] . " to RefSeq ID\n")
		unless exists $reverse_map->{$columns[0]};
	}
	
	if ($mod eq 'FB') {
	    $line = change_FB_transgene_exons($line, \@columns);
	}
	elsif ($mod eq 'WB') {
	    $line = make_wb_changes($line, \@columns);
	}
	elsif ($mod eq 'ZFIN') {
	    $line = make_zfin_changes($line, \@columns);
	}
   
	if ($columns[2] eq 'guide_RNA' or $columns[2] eq 'scRNA') {
	    $columns[2] = 'ncRNA';
	    $line = join("\t", @columns);
	}
	elsif ($columns[2] eq 'protein_coding_gene') {
	    $columns[2] = 'gene';
	    $line = join("\t", @columns);
	}
	
	my @biotypes = ('lnc_RNA', 'lincRNA', 'lincRNA_gene', 'miRNA', 'miRNA_gene', 'pre_miRNA', 
			'mt_gene', 'nc_primary_transcript', 'NMD_transcript_variant', 'pseudogene',
			'processed_pseudogene', 'processed_transcript', 'RNA', 'rRNA', 'rRNA_gene',
			'snoRNA', 'snoRNA_gene', 'snRNA', 'snRNA_gene', 'transcript', 'ncRNA',
			'tRNA', 'pseudogenic_transcript'
	    );
	
	for my $biotype (@biotypes) {
	    if ($columns[2] eq $biotype and $line !~ /biotype=/) {
		if ($line =~ /;Parent/) {
		    my $replacement = ';biotype=' . $biotype . ';Parent';
		    $line =~ s/;Parent/$replacement/;
		}
		else {
		    my $replacement = 'biotype=' . $biotype . ';Parent=';
		    $line =~ s/Parent=/$replacement/;
		}
		last;
	    }
	}

	if ($mod eq 'HUMAN' and $external_human_gff) {
	    $line = convert_to_hgnc_gene_ids(\@columns, $hgnc_id_map);
	}
	
	print OUT $line . "\n";
    }
    
    close(IN);
    close(OUT);
    
    return;
}


sub backup_pathogenicity_prediction_db {
    my ($mod, $password, $log) = @_;
    
    my $date = localtime->strftime('%Y%m%d');
    my $dump_file = join('.', $mod, $date, 'sql');
    my $dump_dir = $ENV{'PATH_PRED_SQL_DUMP_DIR'};
    my $dump_cmd = 'mysqldump -h ' . $ENV{'WORM_DBHOST'} . ' -u ' . $ENV{'WORM_DBUSER'} .
	' -P ' . $ENV{'WORM_DBPORT'} . ' -p' . $password . ' ' .$ENV{'PATH_PRED_DB_PREFIX'} . 
	$mod . ' > ' . $dump_dir . '/' . $dump_file;
    run_slurm_job($dump_cmd, "Dumping $mod pathogenicity predictions database to $dump_dir", $log, '00:30:00', 4, '/dev/null', '/dev/null');
    run_slurm_job("gzip -f -9 $dump_dir/$dump_file",
		   "Compressing $mod pathogenicity predictions database dump", $log, '00:30:00', 4, '/dev/null', '/dev/null');
    
    return;
}

sub change_FB_transgene_exons {
    my ($line, $columns) = @_;
    
    # Change strand of transgene exons to prevent errors
    if ($line =~ /FBtr0084079/ or $line =~ /FBtr008408[345]/ or 
	$line =~ /FBtr0307759/ or $line =~ /FBtr0307760/ or
	$line =~ /FBtr008408[012]/) {
	$columns->[6] = '-';
	$line = join("\t", @$columns);
    }
    
    return $line;
}


sub convert_to_hgnc_gene_ids {
    my ($columns, $hgnc_id_map) = @_;

    my %attr = split /[=;]/, $columns->[8];
    my @pairs;
    my %keys_to_map = map {$_ => 1} ('ID', 'Parent', 'gene_id', 'id');
    # Dbxref=GeneID:100996442,Genbank:XR_001737578.2
    for my $key (keys %attr) {
	my $value = $attr{$key};
	if ($key eq 'Dbxref') {
	    my @dbxrefs = split ',', $value;
	    my @dbxrefs_to_keep;
	    for my $dbxref (@dbxrefs) {
		push @dbxrefs_to_keep, $dbxref unless $dbxref =~ /^GeneID:/
	    }
	    push @pairs, 'Dbxref=' . join(',', @dbxrefs_to_keep) if @dbxrefs_to_keep;
	    next;
	}
	$value =~ s/^gene[:\-]//;
	if (exists $keys_to_map{$key} and exists $hgnc_id_map->{$value}) {
	    push @pairs, $key . '=' . $hgnc_id_map->{$value};
	}
	else {
	    push @pairs, $key . '=' . $value;
	}
    }
    $columns->[8] = join(';', @pairs);

    return join("\t", @$columns);
}


sub make_wb_changes {
    my ($line, $columns) = @_;
    
    my @from = ('piRNA', 'nc_primary_transcript', 'miRNA_primary_transcript', 'pseudogenic_tRNA');
    my @to = ('ncRNA', 'ncRNA', 'miRNA', 'pseudogenic_transcript');
    for my $ix (0 .. @from - 1) {
	if ($columns->[2] eq $from[$ix]) {
	    $columns->[2] = $to[$ix];
	    $line = join("\t", @$columns);
	    last;
	}
    }
    
    return $line;
}


sub make_zfin_changes {
    my ($line, $columns) = @_;
    
    if (($line =~ /protein_coding_gene/ or $columns->[2] eq 'mRNA') and $line !~ /biotype=/) {
	$columns->[8] =~ s/\n/;biotype=protein_coding\n/;
	$line = join("\t", @$columns);
    }
    
    if ($columns->[2] eq 'gene' and $columns->[5] == 1) {
	$columns->[5] = '.';
	$line = join("\t", @$columns);
    }
    
    return $line;
}


sub print_usage {
print <<USAGE;
run_agr_vep_pipelines.pl options:
    -password PASSWORD      password for MySQL database server
    -mods MODS              comma-separated string of MODs to process
    -stages STAGES          comma-separated string of stages to carry out
                            1 = download files from AGR
			    2 = process GFF and FASTA files
			    3 = update pathogenicity prediction score databases
			    4 = run VEP on phenotypic variations
			    5 = run VEP on HTP variations
    -test                   do not upload generated files to AGR
    -logfile                filename to write log to
    -debug                  specify recipients of log email
    -url                    URL of FMS snapshot to use if latest files are not desired
    -cleanup                delete intermediate files generated by pipeline
    -overwrite              overwrite previously downloaded input files
    -nocheck                run analyses even if no new data
    -external_human_gff     do not use the RGD-supplied GFF for human data
    -help                   print this message
USAGE
    
exit 1;
}


sub run_system_cmd {
    my ($cmd, $description, $log) = @_;
    
    $log->write_to("$description\n\n");
    
    my $error = system($cmd);
    if ($error) {
	$log->log_and_die("$description failed: $cmd (Exit code: $error)\n");
    }
    
    return;
}

sub run_slurm_job {
    my ($cmd, $description, $log, $time, $gb_mem, $stdout, $stderr) = @_;

    $log->write_to("$description\n\n");
    my $mem = $gb_mem * 1000;
    $mem .= 'm';
    my $job_id = WormSlurm::submit_job_and_wait($cmd, 'production', $mem, $time, $stdout, $stderr);
    print "Running ${cmd} on Slurm with job ID ${job_id}\n";
    my $exit_code = WormSlurm::get_exit_code($job_id);
    if ($exit_code != 0) {
	$log->log_and_die("Job ${job_id} terminated with exit code ${exit_code} ($cmd)\n") unless $cmd =~ /samtools merge/;
    }
    return;
}


sub check_chromosome_map {
    my ($mod, $log) = @_;

    unless (-e "${mod}_VARIATION.json") {
	$log->write_to("WARNING: No variations file for $mod - cannot check RefSeq chromosome IDs are correct\n");
	return;
    }
    
    my $variation_json = read_file("${mod}_VARIATION.json"); 
    my $variations = decode_json $variation_json;
    
    for my $variation (@{$variations->{data}}) {
	my $refseq_chr = $variation->{sequenceOfReferenceAccessionNumber};
	next unless $refseq_chr =~ /^RefSeq:(.+)$/;
	if ($1 ne $REFSEQ_CHROMOSOMES{$mod}{$variation->{chromosome}}) {
	   $log->log_and_die("RefSeq chromosome ID does not match version submitted in variation JSON for $mod: " .
			     $variation->{chromosome} . " found $1 but should be " . 
			     $REFSEQ_CHROMOSOMES{$mod}{$variation->{chromosome}} . ' for variation ' .
			     $variation->{alleleId} . "\n");
	}
    }
	
    return;
}
    

sub convert_fasta_headers {
    my ($mod, $log) = @_;

    return if -e "${mod}_FASTA.refseq.fa";

    print "Converting $mod FASTA chromosome IDs to RefSeq\n";
    open (IN, '<', "${mod}_FASTA.fa") or $log->log_and_die("Cannot open ${mod}_FASTA.fa for reading\n");
    open (OUT, '>', "${mod}_FASTA.refseq.fa") or $log->log_and_die("Cannot open ${mod}_FASTA.refseq.fa for writing\n");
    while (<IN>) {
	if ($_ =~ /^>(\S+)/) {
	    if (exists $REFSEQ_CHROMOSOMES{$mod}{$1}) {
		print OUT '>' . $REFSEQ_CHROMOSOMES{$mod}{$1} . "\n";
		next;
	    }
	}
	print OUT $_;
    }
    close (IN);
    close (OUT);

    return;
}


sub convert_vcf_chromosomes {
    my ($mod, $type, $log) = @_;

    return unless -e "${mod}_${type}.vcf";
    return if -e "${mod}_${type}.refseq.vcf";

    my $reverse_map = get_reverse_chromosome_map($mod);
    
    print "Converting $mod $type chromosome IDs to RefSeq\n";
    open (IN, '<', "${mod}_${type}.vcf") or $log->log_and_die("Cannot open ${mod}_${type}.vcf for reading\n");
    open (OUT, '>', "${mod}_${type}.refseq.vcf") or $log->log_and_die("Cannot open ${mod}_${type}.refseq.vcf for writing\n");
    while (<IN>) {
	if ($_ !~ /^#/) {
	    my @columns = split("\t", $_);
	    if (exists $REFSEQ_CHROMOSOMES{$mod}{$columns[0]}) {
		$columns[0] = $REFSEQ_CHROMOSOMES{$mod}{$columns[0]};
		print OUT join("\t", @columns);
		next;
	    }
	    else {
		$log->log_and_die("Could not map $mod chromosome in $type " . $columns[0] . " to RefSeq ID\n")
		    unless exists $reverse_map->{$columns[0]};
	    }
	}
	print OUT $_;
    }
    close (IN);
    close (OUT);

    return;
}


sub get_reverse_chromosome_map {
    my $mod = shift;

    my %reverse_map;
    for my $chr (keys %{$REFSEQ_CHROMOSOMES{$mod}}) {
	$reverse_map{$REFSEQ_CHROMOSOMES{$mod}{$chr}} = $chr;
    }

    return \%reverse_map;
}


sub cleanup_intermediate_files {
    my ($mod, $log) = @_;

    for my $file ("${mod}_GFF.refseq.gff", "${mod}_GFF.refseq.gff.gz", "${mod}_GFF.refseq.gff.gz.tbi", 
		  "${mod}_FASTA.refseq.fa", "${mod}_FASTA.refseq.fa.gz", "${mod}_FASTA.refseq.fa.gz",
		  "${mod}_FASTA.refseq.fa.gz.fai", "${mod}_FASTA.refseq.fa.gz.gzi",
		  "${mod}_VCF.refseq.vcf") {
	run_system_cmd("rm $file", "Deleting $file", $log) if -e $file;
    }

    return;
}
 

sub get_hgnc_id_map {
    my $log = shift;

    my ($file) = $HGNC_FILE_URL =~ /\/([^\/]+)$/;
    run_system_cmd("curl -O $HGNC_FILE_URL", 'Downloading HGNC gene ID map file', $log) unless -e $file;

    my $first_line = 1;
    open (HGNC, '<', $file) or $log->log_and_die("Cannot open $file for reading\n");
    my %ids_to_store = map {$_ => 1} @IDS_TO_MAP;
    my @ix_to_store;
    my %id_map;
    while (<HGNC>) {
	chomp;
	my @columns = split("\t", $_);
	if ($first_line) {
	    for my $ix (1 .. scalar @columns - 1) {
		push @ix_to_store, $ix if exists $ids_to_store{$columns[$ix]};
	    }
	    $first_line = 0;
	    next;
	}
	for my $ix (@ix_to_store) {
	    $id_map{$columns[$ix]} = $columns[0] unless $columns[$ix] eq '';
	}
    }
    close (HGNC);
    
    return \%id_map;
}
