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

my ($url, $test, $password, $help);
my $stages = '1,2,3,4,5';
my $mods_string = 'FB,MGI,RGD,SGD,WB,ZFIN,HUMAN';

GetOptions(
    "mods|m=s"     => \$mods_string,
    "stages|s=s"   => \$stages,
    "password|p=s" => \$password,
    "url|u=s"      => \$url,
    "test|t"       => \$test,
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
}


sub download_from_agr {
    my ($mods, $url) = @_;

    my $download_urls = defined $url ? get_urls_from_snapshot($url) : get_latest_urls();
    
    for my $mod (@$mods) {
	my $mod_dir = "$BASE_DIR/$mod";
	make_path($mod_dir) unless -d $mod_dir;
	chdir $mod_dir;
	for my $datatype (keys %{$download_urls->{$mod}}){
	    my ($filename) = $download_urls->{$mod}{$datatype} =~ /\/([^\/]+)$/;
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

    my $chr_map;
    unless ($mod eq 'HUMAN') {
	$chr_map = create_chromosome_map($mod);
	convert_fasta_headers($mod, $chr_map);
	convert_vcf_chromosomes($mod, $chr_map, 'VCF');
	convert_vcf_chromosomes($mod, $chr_map, 'HTVCF') if -e "${mod}_HTVCF.vcf";
    }
    munge_gff($mod, $chr_map);
    run_system_cmd("bgzip -c ${mod}_FASTA.fa > ${mod}_FASTA.fa.gz", "Compressing $mod FASTA");
    run_system_cmd("sort -k1,1 -k4,4n -k5,5n -t\$'\\t' ${mod}_GFF.gff | bgzip -c > ${mod}_GFF.gff.gz",
		       "Sorting and compressing $mod GFF");
    run_system_cmd("tabix -p gff ${mod}_GFF.gff.gz", "Indexing $mod GFF");
 
    return;
}


sub calculate_pathogenicity_predictions {
    my ($mod, $password, $test) = @_;

    backup_pathogenicity_prediction_db($mod, $password);
    
    my $lsf_queue = $test ? $ENV{'LSF_TEST_QUEUE'} : $ENV{'LSF_DEFAULT_QUEUE'};

    my $bam = 0;
    $bam = "${mod}_BAM.bam" if -e "${mod}_BAM.bam";
    
    my $init_cmd = "init_pipeline.pl VepProteinFunction::VepProteinFunction_conf -mod $mod" .
	" -agr_fasta ${mod}_FASTA.fa -agr_gff ${mod}_GFF.gff -agr_bam ${mod}_BAM.bam" . 
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

    return unless -e "${mod}_VCF.vcf";

    my $base_vep_cmd = "vep -i ${mod}_VCF.vcf -gff ${mod}_GFF.gff.gz -fasta ${mod}_FASTA.fa.gz --force_overwrite " .
	"-hgvsg -hgvs -shift_hgvs=0 --symbol --distance 0 --plugin ProtFuncSeq,mod=$mod,pass=$password";
    $base_vep_cmd .= " --bam ${mod}_BAM.bam" if -e "${mod}_BAM.bam";
    my $gl_vep_cmd = $base_vep_cmd . " --per_gene --output_file ${mod}_VEPGENE.txt";
    my $tl_vep_cmd = $base_vep_cmd . " --output_file ${mod}_VEPTRANSCRIPT.txt";

    run_system_cmd($gl_vep_cmd, "Running VEP for $mod phenotypic variants at gene level");
    run_system_cmd($tl_vep_cmd, "Running VEP for $mod phenotypic variants at transcript level");

    for my $level ('GENE', 'TRANSCRIPT') {
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

    return unless -e "${mod}_HTVCF.vcf";

    my $bam = -e "${mod}_BAM.bam" ? "${mod}_BAM.bam" : 0;
    my $lsf_queue = $test ? $ENV{'LSF_TEST_QUEUE'} : $ENV{'LSF_DEFAULT_QUEUE'};

    my $init_cmd = "init_pipeline.pl ModVep::ModVep_conf -mod $mod -vcf ${mod}_HTVCF.vcf -gff ${mod}_GFF.gff.gz" .
	" -fasta ${mod}_FASTA.fa.gz -bam $bam -hive_root_dir " . $ENV{'HIVE_ROOT_DIR'} . ' -pipeline_base_dir ' .
	$ENV{'HTP_VEP_WORKING_DIR'} . ' -pipeline_host ' . $ENV{'WORM_DBHOST'} . ' -pipeline_user ' . $ENV{'WORM_DBUSER'} .
	' -pipeline_port ' . $ENV{'WORM_DBPORT'} . ' -lsf_queue ' . $lsf_queue . ' -vep_dir ' . $ENV{'VEP_DIR'} . 
	" -debug_mode 0 -password $password";

    run_system_cmd($init_cmd, "Initialising $mod HTP variants VEP eHive pipeline");
    my $ehive_url = 'mysql://' . $ENV{'WORM_DBUSER'} . ':' . $password . '@' . $ENV{'WORM_DBHOST'} . ':' . 
	$ENV{'WORM_DBPORT'} . '/agr_htp_' . lc($mod) . '_vep_ehive';
    $ENV{EHIVE_URL} = $ehive_url;
    run_system_cmd("beekeeper.pl -url $ehive_url -loop", "Running $mod HTP variations VEP eHive pipeline");
    run_system_cmd('mv ' . $ENV{'HTP_VEP_WORKING_DIR'} . "/${mod}_vep/${mod}.vep.vcf ${mod}_HTPOSTVEPVCF.vcf",
		   "Moving $mod HTP variations VEP output");
    run_system_cmd("gzip -9 ${mod}_HTPOSTVEPVCF.vcf", "Compressing $mod HTP variations VEP output");
    my $curl_cmd = 'curl -H "Authorization: Bearer ' . $ENV{'TOKEN'} . 
	'" -X POST "https://fms.alliancegenome.org/api/data/submit" -F "' . $ENV{'AGR_RELEASE'} .
	'_HTPOSTVEPVCF' . '_' . $mod . '=@' . $mod . '_HTPOSTVEPVCF.txt.gz"';
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
	return;
    }

    run_system_cmd("samtools sort -o ${mod}_BAM.sorted.bam -T tmp ${mod}_BAM.bam", "Sorting $mod BAM file");
    run_system_cmd("mv ${mod}_BAM.sorted.bam ${mod}_BAM.bam", "Replacing $mod BAM file with sorted version");

    return;
}


sub sort_vcf_files {
    my $mod = shift;
    
    for my $datatype ('VCF', 'HTVCF') {
	next unless -e "${mod}_${datatype}.vcf";
	run_system_cmd("vcf-sort ${mod}_${datatype}.vcf > ${mod}_${datatype}.sorted.vcf",
		       "Sorting $mod $datatype file");
	run_system_cmd("mv ${mod}_${datatype}.sorted.vcf ${mod}_${datatype}.vcf",
		       "Replacing sorted $mod $datatype file with sorted version");
    }

    return;
}


sub munge_gff {
    my ($mod, $chr_map) = @_;

    my $gff = "${mod}_GFF.gff";
    open(IN, "grep -v '^#' $gff |") or die "Could not open $gff for reading\n";
    open(OUT, "> $gff.tmp") or die "Could not open $gff.tmp for writing\n";

    while (<IN>) {
	my $line = $_;
	next if $line =~ /^\n$/;

	my @columns = split("\t", $_);
	if ($mod ne 'HUMAN') {
	    if (exists $chr_map->{$columns[0]}) {
		$columns[0] = $chr_map->{$columns[0]};
		$line = join("\t", @columns);
	    }
	    else {
		die "Could not map $mod chromosome in GFF " . $columns[0] . "to RefSeq ID\n";
	    }
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
    
    run_system_cmd("mv $gff.tmp $gff");

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


sub create_chromosome_map {
    my ($mod) = shift;

    my %chromosome_map;
    my $variation_json = read_file("${mod}_VARIATION.json"); 
    my $variations = decode_json $variation_json;

    for my $variation (@{$variations->{data}}) {
	my $refseq_chr = $variation->{sequenceOfReferenceAccessionNumber};
	next unless $refseq_chr =~ /^RefSeq:(.+)$/;
	$chromosome_map{$variation->{chromosome}} = $1;
    }

    open (MAP, '>', "${mod}_chromosome_map.txt");
    for my $chr (sort values %chromosome_map) {
	print MAP $chr . "\t" . $chromosome_map{$chr} . "\n";
    }
    close (MAP);

    return \%chromosome_map;
}
    

sub convert_fasta_headers {
    my ($mod, $chr_map) = @_;

    open (IN, '<', "${mod}_FASTA.fa");
    open (OUT, '>', "${mod}_FASTA.fa.tmp");
    while (<IN>) {
	if ($_ =~ /^>(\S+)/) {
	    if (exists $chr_map->{$1}) {
		print OUT '>' . $chr_map->{$1} . "\n";
		next;
	    }
	}
	print OUT $_;
    }
    close (IN);
    close (OUT);

    run_system_cmd("mv ${mod}_FASTA.fa.tmp ${mod}_FASTA.fa",
		   'Replacing FASTA file with RefSeq chr ID version');

    return;
}


sub convert_vcf_chromosomes {
    my ($mod, $chr_map, $type) = @_;

    open (IN, '<', "${mod}_${type}.vcf");
    open (OUT, '>', "${mod}_${type}.vcf.tmp");
    while (<IN>) {
	if ($_ !~ /^#/) {
	    my @columns = split("\t", $_);
	    if (exists $chr_map->{$columns[0]}) {
		$columns[0] = $chr_map->{$columns[0]};
		print OUT join("\t", @columns);
		next;
	    }
	    else {
		die "Could not map $mod chromosome in VCF " . $columns[0] . "to RefSeq ID\n";
	    }
	}
	print OUT $_;
    }
    close (IN);
    close (OUT);

    run_system_cmd("mv ${mod}_${type}.vcf.tmp ${mod}_${type}.vcf",
		   "Replacing $type file with RefSeq chr ID version");

    return;
}
