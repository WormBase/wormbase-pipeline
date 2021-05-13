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

const my $FMS_LATEST_PREFIX => 'https://fms.alliancegenome.org/api/datafile/by/';
const my $FMS_LATEST_SUFFIX => '?latest=true';
const my @FMS_DATATYPES => ('FASTA', 'GFF', 'VCF', 'HTVCF', 'MOD-GFF-BAM-KNOWN', 'MOD-GFF-BAM-MODEL', 'VARIATION');
const my %DATATYPE_EXTENSIONS => ('FASTA'             => 'fa',
				  'GFF'               => 'gff',
				  'VCF'               => 'vcf',
				  'HTVCF'             => 'vcf',
				  'MOD-GFF-BAM-KNOWN' => 'bam',
				  'MOD-GFF-BAM-MODEL' => 'bam',
				  'VARIATION'         => 'json',
    );
const my %ASSEMBLIES => ('GRCm38'   => 'MGI',
			 'R6'       => 'FB',
			 'R627'     => 'FB',
			 'Rnor60'   => 'RGD',
			 'SGDr64'   => 'SGD',
			 'WBcel235' => 'WB',
			 'GRCz11'   => 'ZFIN',
    );
const my $BASE_DIR => $ENV{'AGR_VEP_BASE_DIR'} . '/' . $ENV{'AGR_RELEASE'} . '/' . $ENV{'DOWNLOAD_DATE'};

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
	'1'  => 'NC_000067.6',
	'2'  => 'NC_000068.7',
	'3'  => 'NC_000069.6',
	'4'  => 'NC_000070.6',
	'5'  => 'NC_000071.6',
	'6'  => 'NC_000072.6',
	'7'  => 'NC_000073.6',
	'8'  => 'NC_000074.6',
	'9'  => 'NC_000075.6',
	'10' => 'NC_000076.6',
	'11' => 'NC_000077.6',
	'12' => 'NC_000078.6',
	'13' => 'NC_000079.6',
	'14' => 'NC_000080.6',
	'15' => 'NC_000081.6',
	'16' => 'NC_000082.6',
	'17' => 'NC_000083.6',
	'18' => 'NC_000084.6',
	'19' => 'NC_000085.6',
	'X'  => 'NC_000086.7',
	'Y'  => 'NC_000087.7',
	'MT' => 'NC_005089.1',
    },
    'RGD'  => {
	'1'  => 'NC_005100.4',
	'2'  => 'NC_005101.4',
	'3'  => 'NC_005102.4',
	'4'  => 'NC_005103.4',
	'5'  => 'NC_005104.4',
	'6'  => 'NC_005105.4',
	'7'  => 'NC_005106.4',
	'8'  => 'NC_005107.4',
	'9'  => 'NC_005108.4',
	'10' => 'NC_005109.4',
	'11' => 'NC_005110.4',
	'12' => 'NC_005111.4',
	'13' => 'NC_005112.4',
	'14' => 'NC_005113.4',
	'15' => 'NC_005114.4',
	'16' => 'NC_005115.4',
	'17' => 'NC_005116.4',
	'18' => 'NC_005117.4',
	'19' => 'NC_005118.4',
	'20' => 'NC_005119.4',
	'X'  => 'NC_005120.4',
	'Y'  => 'NC_024475.1',
	'MT' => 'NC_001665.2',
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

my ($url, $test, $password, $cleanup, $help);
my $stages = '1,2,3,4,5';
my $mods_string = 'FB,MGI,RGD,SGD,WB,ZFIN,HUMAN';

GetOptions(
    "mods|m=s"     => \$mods_string,
    "stages|s=s"   => \$stages,
    "password|p=s" => \$password,
    "url|u=s"      => \$url,
    "test|t"       => \$test,
    "cleanup|c"    => \$cleanup,
    "help|h"       => \$help,
    ) or print_usage();

print_usage() if $help;

my @mods = split(',', $mods_string);

download_from_agr(\@mods, $url) if $stages =~ /1/;

for my $mod (@mods) {
    chdir "$BASE_DIR/$mod";
    process_input_files($mod) if $stages =~ /2/;
    calculate_pathogenicity_predictions($mod, $password, $test) if $stages =~ /3/;
    run_vep_on_phenotypic_variations($mod, $password, $test) if $stages =~ /4/;
    run_vep_on_htp_variations($mod, $password, $test) if $stages =~ /5/;
    cleanup_intermediate_files($mod) if $cleanup;
}


sub download_from_agr {
    my ($mods, $url) = @_;
    
    my $download_urls = defined $url ? get_urls_from_snapshot($url) : get_latest_urls();
    
    make_path($BASE_DIR) unless -d $BASE_DIR;
    my $input_files_file = "${BASE_DIR}/VEP_input_files.txt";
    open (FILES, '>>', $input_files_file) or die $!;
    for my $mod (@$mods) {
	my $time = localtime();
	print FILES "${mod}: $time\n";
	my $mod_dir = "$BASE_DIR/$mod";
	make_path($mod_dir) unless -d $mod_dir;
	chdir $mod_dir;
	for my $datatype (keys %{$download_urls->{$mod}}){
	    my ($filename) = $download_urls->{$mod}{$datatype} =~ /\/([^\/]+)$/;
	    print FILES "\t${datatype}: ${filename}\n";
	    run_system_cmd('wget ' . $download_urls->{$mod}{$datatype}, "Downloading $mod $datatype file");
	    $filename = check_if_actually_compressed($filename) if $filename !~ /\.gz/; # temporary hack to get around gzipped files in FMS without .gz extension
	    run_system_cmd("gunzip $filename", "Decompressing $filename") if $filename =~ /\.gz$/; # if clause only required in interim while some FMS files not gzipped
	    $filename =~ s/\.gz$//;
	    my $extension = $DATATYPE_EXTENSIONS{$datatype};
	    run_system_cmd("mv $filename ${mod}_${datatype}.${extension}", "Renaming $filename");
	}
	merge_bam_files($mod);
	sort_vcf_files($mod);
    }
    close (FILES);
    
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
    my $filename = shift;
    
    my $is_gzipped = 0;
    open (FILE, "file $filename |");
    while (<FILE>) {
	chomp;
	$is_gzipped = 1 if $_ =~ /gzip compressed data/;
    }
    close (FILE);
    
    if ($is_gzipped) {
	run_system_cmd("mv $filename $filename.gz", "Adding .gz extension to $filename");
	return $filename . '.gz';
    }
    
    return $filename;
}


sub process_input_files {
    my $mod = shift;
    
    cleanup_intermediate_files($mod);

    my $chr_map;
    check_chromosome_map($mod);
    convert_fasta_headers($mod);
    convert_vcf_chromosomes($mod, 'VCF');
  
    munge_gff($mod);
    run_system_cmd("bgzip -c ${mod}_FASTA.refseq.fa > ${mod}_FASTA.refseq.fa.gz", "Compressing $mod FASTA");
    run_system_cmd("sort -k1,1 -k4,4n -k5,5n -t\$'\\t' ${mod}_GFF.refseq.gff | bgzip -c > ${mod}_GFF.refseq.gff.gz",
    	             "Sorting and compressing $mod GFF");
    run_system_cmd("tabix -p gff ${mod}_GFF.refseq.gff.gz", "Indexing $mod GFF");
    
    return;
}


sub calculate_pathogenicity_predictions {
    my ($mod, $password, $test) = @_;
    
    backup_pathogenicity_prediction_db($mod, $password);
    
    my $lsf_queue = $test ? $ENV{'LSF_TEST_QUEUE'} : $ENV{'LSF_DEFAULT_QUEUE'};
        
    my $init_cmd = "init_pipeline.pl VepProteinFunction::VepProteinFunction_conf -mod $mod" .
	" -agr_fasta ${mod}_FASTA.refseq.fa -agr_gff ${mod}_GFF.refseq.gff -agr_bam ${mod}_BAM.bam" . 
	' -hive_root_dir ' . $ENV{'HIVE_ROOT_DIR'} . ' -pipeline_base_dir ' . $ENV{'PATH_PRED_WORKING_DIR'} .
	' -pipeline_host ' . $ENV{'WORM_DBHOST'} . ' -pipeline_user ' . $ENV{'WORM_DBUSER'} .
	' -pipeline_port ' . $ENV{'WORM_DBPORT'} . ' -lsf_queue ' . $lsf_queue .
	' -sift_dir ' . $ENV{'SIFT_DIR'} . ' -pph_dir ' . $ENV{'PPH_DIR'} . ' -pph_conf_dir ' . $ENV{'PPH_CONF_DIR'} . '/' . $mod .
	' -ncbi_dir ' . $ENV{'NCBI_DIR'} . ' -blastdb ' . $ENV{'BLAST_DB'} . ' -password ' . $password;
    
    run_system_cmd($init_cmd, "Initialising $mod pathogenicity prediction eHive pipeline");
    
    my $ehive_url = 'mysql://' . $ENV{'WORM_DBUSER'} . ':' . $password . '@' . $ENV{'WORM_DBHOST'} . ':' . 
	$ENV{'WORM_DBPORT'} . '/agr_pathogenicity_predictions_' . lc($mod) . '_ehive';
    $ENV{EHIVE_URL} = $ehive_url;
    run_system_cmd("beekeeper.pl -url $ehive_url -loop", "Running $mod pathogenicity prediction eHive pipeline");
    
    return;
}


sub run_vep_on_phenotypic_variations {
    my ($mod, $password, $test) = @_;
    
    return unless -e "${mod}_VCF.refseq.vcf";
    
    my $base_vep_cmd = "vep -i ${mod}_VCF.refseq.vcf -gff ${mod}_GFF.refseq.gff.gz -fasta ${mod}_FASTA.refseq.fa.gz --force_overwrite " .
	"--bam ${mod}_BAM.bam -hgvsg -hgvs -shift_hgvs=0 --symbol --distance 0 --plugin ProtFuncSeq,mod=$mod,pass=$password";
   
    my $gl_vep_cmd = $base_vep_cmd . " --per_gene --output_file ${mod}_VEPGENE.txt";
    my $tl_vep_cmd = $base_vep_cmd . " --output_file ${mod}_VEPTRANSCRIPT.txt";
    
    run_system_cmd($gl_vep_cmd, "Running VEP for $mod phenotypic variants at gene level");
    run_system_cmd($tl_vep_cmd, "Running VEP for $mod phenotypic variants at transcript level");
    
    my $reverse_map = get_reverse_chromosome_map($mod);

    for my $level ('GENE', 'TRANSCRIPT') {
	open (IN, '<', "${mod}_VEP${level}.txt");
	open (OUT, '>', "${mod}_VEP${level}.txt.tmp");
	while (<IN>) {
	    if ($_ =~ /^#/) {
		print OUT $_;
	    }
	    else {
		my @columns = split("\t", $_);
		my ($chr, $pos) = split(':', $columns[1]);
		$columns[1] = $reverse_map->{$chr} . ':' . $pos;
		print OUT join("\t", @columns);
	    }
	}
	close (IN);
	close (OUT);
	run_system_cmd("mv ${mod}_VEP${level}.txt.tmp ${mod}_VEP${level}.txt",
		       "Replacing $mod VEP${level} file with original chromosome ID version");

	run_system_cmd("gzip -9 ${mod}_VEP${level}.txt", 'Compressing ' . lc($level) . '-level VEP results');
	my $curl_cmd = 'curl -H "Authorization: Bearer ' . $ENV{'TOKEN'} . 
	    '" -X POST "https://fms.alliancegenome.org/api/data/submit" -F "' . $ENV{'AGR_RELEASE'} .
	    '_VEP' . $level . '_' . $mod . '=@' . $mod . '_VEP' . $level . '.txt.gz"';
	run_system_cmd($curl_cmd, "Uploading $mod " . lc($level) . '-level VEP results to AGR') unless $test;
    }
    
    return;
}


sub run_vep_on_htp_variations{
    my ($mod, $password, $test) = @_;
    
    my $lsf_queue = $test ? $ENV{'LSF_TEST_QUEUE'} : $ENV{'LSF_DEFAULT_QUEUE'};
    
    my $init_cmd = "init_pipeline.pl ModVep::ModVep_conf -mod $mod -vcf ${mod}_HTVCF.vcf -gff ${mod}_GFF.refseq.gff.gz" .
	" -fasta ${mod}_FASTA.refseq.fa.gz -bam ${mod}_BAM.bam -hive_root_dir " . $ENV{'HIVE_ROOT_DIR'} . ' -pipeline_base_dir ' .
	$ENV{'HTP_VEP_WORKING_DIR'} . ' -pipeline_host ' . $ENV{'WORM_DBHOST'} . ' -pipeline_user ' . $ENV{'WORM_DBUSER'} .
	' -pipeline_port ' . $ENV{'WORM_DBPORT'} . ' -lsf_queue ' . $lsf_queue . ' -vep_dir ' . $ENV{'VEP_DIR'} . 
	" -debug_mode 0 -password $password";
    
    run_system_cmd($init_cmd, "Initialising $mod HTP variants VEP eHive pipeline: $init_cmd");
    my $ehive_url = 'mysql://' . $ENV{'WORM_DBUSER'} . ':' . $password . '@' . $ENV{'WORM_DBHOST'} . ':' . 
	$ENV{'WORM_DBPORT'} . '/agr_htp_' . lc($mod) . '_vep_ehive';
    $ENV{EHIVE_URL} = $ehive_url;
 
    run_system_cmd("beekeeper.pl -url $ehive_url -loop", "Running $mod HTP variations VEP eHive pipeline");

    my $uncompressed_file = $ENV{'HTP_VEP_WORKING_DIR'} . "/${mod}_vep/${mod}.vep.vcf";
    my $bsub_cmd = 'bsub -J ' . $mod . '_VEP_compress -o /dev/null -e /dev/null -n 20 -R "span[ptile=20]" ' .
	'pigz -9 -p 20 ' . $uncompressed_file;
    run_system_cmd($bsub_cmd, "Compressing $mod HTP variations VEP output");
    while (-e $uncompressed_file) {
	sleep(60);
    }
    
    my $compressed_file = "${uncompressed_file}.gz";
    run_system_cmd("mv ${compressed_file} ${mod}_HTPOSTVEPVCF.vcf.gz",
		   "Moving $mod HTP variations VEP output");


    my $curl_cmd = 'curl -H "Authorization: Bearer ' . $ENV{'TOKEN'} . 
	'" -X POST "https://fms.alliancegenome.org/api/data/submit" -F "' . $ENV{'AGR_RELEASE'} .
	'_HTPOSTVEPVCF' . '_' . $mod . '=@' . $mod . '_HTPOSTVEPVCF.vcf.gz"';
    run_system_cmd($curl_cmd, "Uploading $mod HTP VEP results to AGR") unless $test;
    
    return;
}


sub merge_bam_files {
    my $mod = shift;
    
    if (-e "${mod}_MOD-GFF-BAM-KNOWN.bam") {
	if (-e "${mod}_MOD-GFF-BAM-MODEL.bam") {
	    run_system_cmd("samtools merge ${mod}_BAM.bam ${mod}_MOD-GFF-BAM-KNOWN.bam ${mod}_MOD-GFF-BAM-MODEL.bam",
			   "Merging $mod BAM files");
	    run_system_cmd("rm ${mod}_MOD-GFF-BAM-KNOWN.bam", "Deleting $mod unmerged known transcripts BAM file");
	    run_system_cmd("rm ${mod}_MOD-GFF-BAM-MODEL.bam", "Deleting $mod unmerged model transcripts BAM file");
	}
	else{
	    run_system_cmd("mv ${mod}_MOD-GFF-BAM-KNOWN.bam ${mod}_BAM.bam", "Renaming $mod MOD-GFF-BAM-KNOWN file");
	}
    }
    elsif (-e "${mod}_MOD-GFF-BAM-MODEL.bam") {
	run_system_cmd("mv ${mod}_MOD-GFF-BAM-MODEL.bam ${mod}_BAM.bam", "Renaming $mod MOD-GFF-BAM-MODEL file"); 
    }
    else {
	run_system_cmd("touch ${mod}_BAM.bam", "Creating dummy BAM file");
	return;
    }
    
    run_system_cmd("samtools sort -o ${mod}_BAM.sorted.bam -T tmp ${mod}_BAM.bam", "Sorting $mod BAM file");
    run_system_cmd("mv ${mod}_BAM.sorted.bam ${mod}_BAM.bam", "Replacing $mod BAM file with sorted version");
    run_system_cmd("samtools index ${mod}_BAM.bam", "Indexing $mod BAM file");
    
    return;
}


sub sort_vcf_files {
    my $mod = shift;
    
    for my $datatype ('VCF', 'HTVCF') {
	next unless -e "${mod}_${datatype}.vcf";
	run_system_cmd("vcf-sort ${mod}_${datatype}.vcf > ${mod}_${datatype}.sorted.vcf",
		       "Sorting $mod $datatype file");
	run_system_cmd("mv ${mod}_${datatype}.sorted.vcf ${mod}_${datatype}.vcf",
		       "Replacing unsorted $mod $datatype file with sorted version");
    }
    
    return;
}


sub remove_mirna_primary_transcripts {
    my (%mirna_parents, %parents);

    print "Getting RefSeq HUMAN miRNA details from GFF\n";

    # Do first pass to get parents
    open (GFF, "grep -v '^#' HUMAN_GFF.gff |") or die "Could not open HUMAN_GFF.gff for reading\n";
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

    open (IN, '< HUMAN_GFF.gff') or die "Could not open HUMAN_GFF.gff for reading\n";
    open (OUT, '> HUMAN_GFF.tmp.gff') or die "Could not open HUMAN_GFF.tmp.gff for writing\n";
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

    run_system_cmd("mv HUMAN_GFF.tmp.gff HUMAN_GFF.gff", "Replacing GFF with version without pre-miRNA lines");

    return;
}


sub munge_gff {
    my ($mod) = @_;
    
    run_system_cmd("rm ${mod}_GFF.refseq.gff", "Deleting old munged GFF file") if -e "${mod}_GFF.refseq.gff";

    my $reverse_map = get_reverse_chromosome_map($mod);

    remove_mirna_primary_transcripts() if $mod eq 'HUMAN';

    print "Munging $mod GFF\n";
    open(IN, "grep -v '^#' ${mod}_GFF.gff |") or die "Could not open ${mod}_GFF.gff for reading\n";
    open(OUT, "> ${mod}_GFF.refseq.gff") or die "Could not open ${mod}_GFF.refseq.gff for writing\n";
    while (<IN>) {
	my $line = $_;
	next if $line =~ /^\n$/;
	
	my @columns = split("\t", $_);
	if (exists $REFSEQ_CHROMOSOMES{$mod}{$columns[0]}) {
	    $columns[0] = $REFSEQ_CHROMOSOMES{$mod}{$columns[0]};
	    $line = join("\t", @columns);
	}
	else {
	    die "Could not map $mod chromosome in GFF " . $columns[0] . " to RefSeq ID\n"
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
	
	print OUT $line;
    }
    
    close(IN);
    close(OUT);
    
    return;
}


sub backup_pathogenicity_prediction_db {
    my ($mod, $password) = @_;
    
    my $date = localtime->strftime('%Y%m%d');
    my $dump_file = join('.', $mod, $date, 'sql');
    my $dump_dir = $ENV{'PATH_PRED_SQL_DUMP_DIR'};
    my $dump_cmd = 'mysqldump -h ' . $ENV{'WORM_DBHOST'} . ' -u ' . $ENV{'WORM_DBUSER'} .
	' -P ' . $ENV{'WORM_DBPORT'} . ' -p' . $password . ' ' .$ENV{'PATH_PRED_DB_PREFIX'} . 
	$mod . ' > ' . $dump_dir . '/' . $dump_file;
    run_system_cmd($dump_cmd, "Dumping $mod pathogenicity predictions database to $dump_dir");
    run_system_cmd("gzip -9 $dump_dir/$dump_file",
		   "Compressing $mod pathogenicity predictions database dump");
    
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
    -cleanup                delete intermediate files generated by pipeline
    -help                   print this message
USAGE
    
exit 1;
}


sub run_system_cmd {
    my ($cmd, $description) = @_;
    
    system("echo $description\n\n");
    
    my $error = system($cmd);
    if ($error) {
	die "$description failed: $cmd (Exit code: $error)\n";
    }
    
    return;
}


sub check_chromosome_map {
    my ($mod) = shift;

    unless (-e "${mod}_VARIATION.json") {
	warn "No variations file for $mod - cannot check RefSeq chromosome IDs are correct\n";
	return;
    }
    
    my $variation_json = read_file("${mod}_VARIATION.json"); 
    my $variations = decode_json $variation_json;
    
    for my $variation (@{$variations->{data}}) {
	my $refseq_chr = $variation->{sequenceOfReferenceAccessionNumber};
	next unless $refseq_chr =~ /^RefSeq:(.+)$/;
	if ($1 ne $REFSEQ_CHROMOSOMES{$mod}{$variation->{chromosome}}) {
	    die "RefSeq chromosome ID does not match version submitted in variation JSON for $mod: " .
		$variation->{chromosome} . " should be $1 but found " . 
		$REFSEQ_CHROMOSOMES{$mod}{$variation->{chromosome}} . "\n";
	}
    }
	
    return;
}
    

sub convert_fasta_headers {
    my ($mod) = @_;

    return if -e "${mod}_FASTA.refseq.fa";

    print "Converting $mod FASTA chromosome IDs to RefSeq\n";
    open (IN, '<', "${mod}_FASTA.fa");
    open (OUT, '>', "${mod}_FASTA.refseq.fa");
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
    my ($mod, $type) = @_;

    return unless -e "${mod}_${type}.vcf";
    return if -e "${mod}_${type}.refseq.vcf";

    my $reverse_map = get_reverse_chromosome_map($mod);
    
    print "Converting $mod $type chromosome IDs to RefSeq\n";
    open (IN, '<', "${mod}_${type}.vcf");
    open (OUT, '>', "${mod}_${type}.refseq.vcf");
    while (<IN>) {
	if ($_ !~ /^#/) {
	    my @columns = split("\t", $_);
	    if (exists $REFSEQ_CHROMOSOMES{$mod}{$columns[0]}) {
		$columns[0] = $REFSEQ_CHROMOSOMES{$mod}{$columns[0]};
		print OUT join("\t", @columns);
		next;
	    }
	    else {
		die "Could not map $mod chromosome in $type " . $columns[0] . " to RefSeq ID\n"
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


sub cleanup_intermediate_files{
    my $mod = shift;

    for my $file ("${mod}_GFF.refseq.gff", "${mod}_GFF.refseq.gff.gz", "${mod}_GFF.refseq.gff.gz.tbi", 
		  "${mod}_FASTA.refseq.fa", "${mod}_FASTA.refseq.fa.gz", "${mod}_FASTA.refseq.fa.gz",
		  "${mod}_FASTA.refseq.fa.gz.fai", "${mod}_FASTA.refseq.fa.gz.gzi", "${mod}_VCF.refseq.vcf",
		  "${mod}_HTVCF.refseq.vcf") {
	run_system_cmd("rm $file", "Deleting $file") if -e $file;
    }

    return;
}
 
